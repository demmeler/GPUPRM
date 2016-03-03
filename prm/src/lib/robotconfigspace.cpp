#include "config.h"

#include "robotconfigspace.h"

#include "cuda_head.h"

#include "collision4.h"
#include "polytope4data.h"
#include "kinematics.h"

#include <vector>
#include <cmath>


//! robot_: Robot Object with Denavit-Hartenberg data
//! polydata_: Object containing compressed Polytope data of robot arm and surrounding
//! mins[ndof]: maximal q values
//! maxs[ndof]: minimal   -"-
//! dq:         edge resolution (stepsize)
template<int ndof>
RobotConfigspace<ndof>::RobotConfigspace(const Robot<ndof>* robot_,
                                         const collision4::polytope4* polys_,  const int* sys_, const int N_,
                                         const float* mins_, const float* maxs_, const float dq_,
                                         const int nbuf_, const int numthreadsmax_):
    kin(robot_)
{
  robot=robot_;
  for(int i=0;i<ndof;++i){
    mins[i]=mins_[i];
    maxs[i]=maxs_[i];
  }
  dq=dq_;

  polylist.polys=polys_;
  polylist.sys=sys_;
  polylist.N=N_;

  polylist.from=0x0;
  polylist.to=0x0;
  polylist.M=0x0;

  nbufres=nbuf_;
  nbufqto=ndof*nbuf_;
  nbufqfrom=ndof*nbuf_;
  nbuftest=nbuf_;

  numthreadsmax=numthreadsmax_;

  devloaded=false;

}

template<int ndof>
RobotConfigspace<ndof>::RobotConfigspace(const Robot<ndof>* robot_,
                                         const collision4::polytope4* polys_,  const int* sys_, const int N_,
                                         const int *from_, const int *to_, const int M_,
                                         const float* mins_, const float* maxs_, const float dq_,
                                         const int nbuf_, const int numthreadsmax_):
    RobotConfigspace(robot_, polys_, sys_, N_, mins_, maxs_, dq_, nbuf_, numthreadsmax_)
{
    polylist.from=from_;
    polylist.to=to_;
    polylist.M=M_;
}



//!initialization function copy polytope and robot data to gpu etc..
template<int ndof>
int RobotConfigspace<ndof>::init()
{
  if(devloaded){
    clear();
  }

  polydata=new collision4::polytope4data;
  polydata->build(polylist.polys,polylist.sys,polylist.N,ndof,polylist.from, polylist.to, polylist.M);

#ifdef CUDA_IMPLEMENTATION

  polydatadev_hostref=new collision4::polytope4data;

  collision4::copy_host_to_device(*polydatadev_hostref,*polydata);
  cudaMalloc((void**)&polydatadev, sizeof(collision4::polytope4data));
  cudaMemcpy((void*)polydatadev, (void*)polydatadev_hostref, sizeof(collision4::polytope4data), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&robotdev, sizeof(Robot<ndof>));
  cudaMemcpy((void*)robotdev,(void*)robot, sizeof(Robot<ndof>), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&qdevbufferfrom, nbufqfrom*sizeof(float));
  cudaMalloc((void**)&testnumdev, nbuftest*sizeof(int));
  cudaMalloc((void**)&testposdev, nbuftest*sizeof(int));
  cudaMalloc((void**)&qdevbufferto, nbufqto*sizeof(float));
  cudaMalloc((void**)&resdevbuffer, nbufres*sizeof(int));
  cudaMalloc((void**)&resdevbufferext, numthreadsmax*sizeof(int));

#else

  resbufferext=new int[numthreadsmax];

#endif

  devloaded=true;

  return 0; //TODO: errorhandling?

}

template<int ndof>
int RobotConfigspace<ndof>::clear()
{
  // TODO
    msg("error: RobotConfigspace<ndof>::clear() not implemented");
    return -1;
}

//!
//! indicator functions on CPU
//!

//! indicator function of obstacles
//! q: length d*N, array of structures: q[N*k+i]= k-th component of i-th q-vector
//! res: length N
template<int ndof>
int RobotConfigspace<ndof>::indicator(const float* q, int* res, const int N, const int offset)
{
  msg("error: not implemented");
  return -1;
}

//!case N=1
template<int ndof>
int RobotConfigspace<ndof>::indicator(const float* q)
{
  if(check_boundaries(q)==1) return 1;
  kin.calculate(q,1);

#if 0
  for(int dof0=0;dof0<=ndof;++dof0) for(int dof1=dof0+1;dof1<=ndof;++dof1){
    int numsys0=polydata->get_numsys(dof0), numsys1=polydata->get_numsys(dof1);
    for(int k0=0;k0<numsys0;++k0) for(int k1=0;k1<numsys1;++k1){
      collision4::polytope4 poly0, poly1;
      polydata->get_polytope(poly0, dof0,k0);
      polydata->get_polytope(poly1, dof1,k1);
      int res=seperating_vector_algorithm(poly0,poly1,kin.trafos[dof0],kin.trafos[dof1]);
      if(res!=0){
        return res;
      }
    }
  }
#else

        //new version, only testing specific pairs


    //! collision algorithm

    for(int k0=0;k0<polydata->N;++k0){
      collision4::polytope4 poly0;
      polydata->get_polytope(poly0, k0);
      int *dest;
      int destnum;
      polydata->get_collision_list(k0,dest,destnum);
      for(int l=0;l<destnum;++l){
          int k1=dest[l];
          collision4::polytope4 poly1;
          polydata->get_polytope(poly1, k1);

          dprintarr(dest,destnum);
          dprintvard(polydata->sys[k0]);
          dprintvard(polydata->sys[k1]);
          dt4print(kin.trafos[polydata->sys[k0]]);
          dt4print(kin.trafos[polydata->sys[k1]]);

          int result=collision4::seperating_vector_algorithm(poly0,poly1,kin.trafos[polydata->sys[k0]],kin.trafos[polydata->sys[k1]]);
          if(result!=0){
            return result;
          }
      }
    }

#endif
  return 0;
}









//!
//! indicator functions on GPU
//!


//! ************************
//! *                      *
//! *   collision kernel   *
//! *                      *
//! ************************

template<int ndof>
#ifdef CUDA_IMPLEMENTATION
__global__ void kernel_indicator2(Robot<ndof>* robot,
                                  collision4::polytope4data* polydata,
                                  float* qs, int offsets,
                                  float* qe, int offsete,
                                  int* res, int* resext,
                                  int* testpos, int* testnum,
                                  int N, int numthreads){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i<numthreads){
#else
void kernel_indicator2(const Robot<ndof>* robot,
                       const collision4::polytope4data* polydata,
                       const float* qs, int offsets,
                       const float* qe, int offsete,
                       int* res, int* resext,
                       const int* testpos, const int* testnum,
                       int N, int numthreads){
  for(int i=0;i<numthreads;++i){
#endif

    //! determine the line in which thread is involved
    int k; //! index of the line
           //! line k is handles by threads testpos[k], ...., testpos[k]+numpos[k]
    for(k=N-1;testpos[k]>i;--k); //binaersuche machen?

    //! calculate q (convex combination)
    float q[ndof];
    float c1,c2;
    if(testnum[k]>0){
      c1=(float)(i-testpos[k])/(float)(testnum[k]-1);
      c2=1.0-c1;
    }else{
      c1=1.0;
      c2=0.0;
    }
    for(int j=0;j<ndof;++j){
      q[j]=c1*qs[k+offsets*j]+c2*qe[k+offsete*j];
    }

    Kinematics<ndof> kin(robot);
    resext[i]=0;
    kin.calculate(&q[0],1);

#if 0

    //! collision algorithm
    for(int dof0=0;dof0<=ndof;++dof0) for(int dof1=dof0+1;dof1<=ndof;++dof1){
      int numsys0=polydata->get_numsys(dof0), numsys1=polydata->get_numsys(dof1);
      for(int k0=0;k0<numsys0;++k0) for(int k1=0;k1<numsys1;++k1){
        collision4::polytope4 poly0, poly1;
        polydata->get_polytope(poly0, dof0, k0);
        polydata->get_polytope(poly1, dof1, k1);
        int result=collision4::seperating_vector_algorithm(poly0,poly1,kin.trafos[dof0],kin.trafos[dof1]);
        if(result!=0){
          resext[i]=result;
        }
      }
    }

#else

        //new version, only testing specific pairs


    //! collision algorithm

    for(int k0=0;k0<polydata->N;++k0){
      collision4::polytope4 poly0;
      polydata->get_polytope(poly0, k0);
      int *dest;
      int destnum;
      polydata->get_collision_list(k0,dest,destnum);
      for(int l=0;l<destnum;++l){
          int k1=dest[l];
          collision4::polytope4 poly1;
          polydata->get_polytope(poly1, k1);
          int result=collision4::seperating_vector_algorithm(poly0,poly1,kin.trafos[polydata->sys[k0]],kin.trafos[polydata->sys[k1]]);
          if(result!=0){
            resext[i]=result;
          }
      }
    }

#endif

    //!reduce resext to res
#ifdef CUDA_IMPLEMENTATION
    //TODO
#else
    if(resext[i]!=0 || i==testpos[k]) res[k]=resext[i];
#endif

  }//if/for
}


//! checks if indicator function = 1 somewhere on the line between (qs[i],qs[i+N]) and (qe[i],qe[i+N]) for all i=1,...,N-1
//! res is return value
//! number of pairs
template<int ndof>
int RobotConfigspace<ndof>::indicator2(const float* qs, const float* qe, int *res, const int N, const int offset)
{
  //! calculate number of threads needed
  std::vector<int> testnum(N,0);
  std::vector<int> testpos(N,0);
  int numthreads=0;
  for(int l=0;l<N;++l){
    testpos[l]=numthreads;
    float norm=0.0;
    for(int i=0;i<ndof;++i){
      float diff=qe[l+offset*i]-qs[l+offset*i];
      norm+=diff*diff;
    }
    norm=sqrt(norm);
    testnum[l]=(int)(norm/dq)+1; //start and end points included
    numthreads+=testnum[l];
  }

  if(numthreads>numthreadsmax){
    msg("error, thread limit reached")
    return -1;
  }


#ifdef CUDA_IMPLEMENTATION
  int BLOCK = 256, GRID = (numthreads + BLOCK - 1)/BLOCK;

  for(int i=0;i<ndof;++i){
      //pointer inkrement in cuda??
    cudaMemcpy((void*)(qdevbufferfrom+nbufqfrom*i),(void*)&(qs[offset*i]), N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)(qdevbufferto+nbufqto*i),(void*)&(qe[offset*i]), N*sizeof(float), cudaMemcpyHostToDevice);
  }

  cudaMemcpy((void*)testposdev,(void*)testpos.data(), N*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy((void*)testnumdev,(void*)testnum.data(), N*sizeof(int), cudaMemcpyHostToDevice);


  kernel_indicator2<ndof><<<GRID,BLOCK>>>(robotdev,polydatadev,qdevbufferfrom,nbufqfrom,qdevbufferto,nbufqto,resdevbuffer,resdevbufferext,testposdev,testnumdev,N, numthreads);


  cudaMemcpy((void*)res,(void*)resdevbuffer,N*sizeof(int), cudaMemcpyDeviceToHost);
#else
  kernel_indicator2<ndof>(robot,polydata,qs,offset,qe,offset,res,resbufferext,testpos.data(),testnum.data(),N, numthreads);
#endif

}










//! same paircheck as above, but with compressed storage:
//! checks pairs: (qs[i],...) ->  (qe(posqe[i]),...) , ...., (qe[posqe[i]+numqe[i]-1],...) for i=0,...,M-1
template<int ndof>
int RobotConfigspace<ndof>::indicator2(const float* qs, const int M, const float* qe, int *res, const int *posqe, const int *numqe, const int offset)
{

  msg("this indicator2 is not implemented!");
  return -1;

}


//! boundary check on CPU

//! structure like indicator function
//! returns if lies in boundaries
template<int ndof>
int RobotConfigspace<ndof>::check_boundaries(const float* q, int* res, const int N, const int offset)
{
  for(int ix=0;ix<N;++ix){
    int resix=0;
    for(int i=0;i<ndof;++i){
      float qi=q[ix+i*offset];
      if(qi<mins[i] || qi>maxs[i]){resix=1; break;}
    }
    res[ix]=resix;
  }
  return 0;
}

//! for N=1
template<int ndof>
int RobotConfigspace<ndof>::check_boundaries(const float* q)
{
  for(int i=0;i<ndof;++i){
    float qi=q[i];
    if(qi<mins[i] || qi>maxs[i]){return 1;}
  }
  return 0;
}




template class RobotConfigspace<2>;











