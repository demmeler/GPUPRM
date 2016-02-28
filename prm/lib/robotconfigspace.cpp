#include "robotconfigspace.h"

#include "cuda_head.h"
#include <vector>
#include <cmath>

//! robot_: Robot Object with Denavit-Hartenberg data
//! polydata_: Object containing compressed Polytope data of robot arm and surrounding
//! mins[ndof]: maximal q values
//! maxs[ndof]: minimal   -"-
//! dq:         edge resolution (stepsize)
template<int ndof>
RobotConfigspace<ndof>::RobotConfigspace(const Robot<ndof>* robot_, const collision4::polytope4* polys_, const int* sys_, const int N_, const float* mins_, const float* maxs_, const float dq_, const int nbuf_):
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

  nbufres=nbuf_;
  nbufqto=ndof*nbuf_;
  nbufqfrom=ndof*nbuf_;
  nbuftest=nbuf_;

  devloaded=false;

}


//!initialization function copy polytope and robot data to gpu etc..
template<int ndof>
int RobotConfigspace<ndof>::init()
{
  if(devloaded){
    clear();
  }

  polydata=new collision4::polytope4data;
  polydata->build(polylist.polys,polylist.sys,polylist.N);

  polydatadev_hostref=new collision4::polytope4data;

  collision4::copy_host_to_device(*polydatadev_hostref,*polydata);
  cudaMalloc((void**)&polydatadev, sizeof(polytope4data));
  cudaMemcpy((void*)polydatadev, (void*)polydatadev_hostref, sizeof(polytope4data), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&robotdev, sizeof(Robot<ndof>));
  cudaMemcpy((void*)robotdev,(void*)robot, sizeof(Robot<ndof>), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&qdevbufferfrom, nbufqfrom*sizeof(float));
  cudaMalloc((void**)&testnumdev, nbuftest*sizeof(int));
  cudaMalloc((void**)&testposdev, nbuftest*sizeof(int));
  cudaMalloc((void**)&qdevbufferto, nbufqto*sizeof(float));
  cudaMalloc((void**)&resdevbuffer, nbufres*sizeof(int));
  cudaMalloc((void**)&resdevbufferext, numthreadsmax*sizeof(int));



  devloaded=true;

}

template<int ndof>
int RobotConfigspace<ndof>::clear()
{
  // TODO
}

//!
//! indicator functions on CPU
//!

//! indicator function of obstacles
//! q: length d*N, array of structures: q[N*k+i]= k-th component of i-th q-vector
//! res: length N
template<int ndof>
int RobotConfigspace<ndof>::indicator(const float* q, int* res, int N)
{
  msg("error: not implemented");
  return -1;
}

//!case N=1
template<int ndof>
int RobotConfigspace<ndof>::indicator(const float* q)
{
  kin.calculate(q,1);
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
  return 0;
}

//!
//! indicator functions on GPU
//!



template<int ndof>
__global__ void kernel_indicator2(Robot<ndof>* robot, collision4::polytope4data* polydata, float* qs, int offsets, float* qe, int offsete, int* res, int* resext, int* testpos, int* testnum, int N, int numthreads){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i<numthreads){
    coll[i]=0;
    Kinematics<ndof> kin(robot);

    //!calculate q
    int k; //! index of line in which thread is involved
           //! line k is handles by threads testpos[k], ...., testpos[k]+numpos[k]
    for(k=N-1;testpos[k]>i;--k);

    //! calculate q (convex combination)
    float q[ndof];
    float c1,c2;
    if(testnum[k]>0){
      c1=(float)(j-testpos[k])/(float)(testnum[k]-1);
      c2=1.0-c2;
    }else{
      c1=1.0;
      c2=0.0;
    }
    for(int j=0;j<ndof;++j){
      q[j]=c1*qs[k+offsets*j]+c2*qe[k+offsete*j];
    }


    kin.calculate(&q[0],1);
    for(int dof0=0;dof0<=ndof;++dof0)
    for(int dof1=dof0+1;dof1<=ndof;++dof1){
      int numsys0=polydata->get_numsys(dof0), numsys1=polydata->get_numsys(dof1);
      for(int k0=0;k0<numsys0;++k0)
      for(int k1=0;k1<numsys1;++k1){
        polytope4 poly0, poly1;
        polydata->get_polytope(poly0, dof0,k0);
        polydata->get_polytope(poly1, dof1,k1);
        int res=seperating_vector_algorithm(poly0,poly1,kin.trafos[dof0],kin.trafos[dof1]);
        if(res!=0){
          coll[i]=res;
        }
      }
    }
  }//if
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

  int BLOCK = 256, GRID = (numthreads + BLOCK - 1)/BLOCK;

  for(int i=0;i<ndof;++i){
      //pointer inkrement in cuda??
    cudaMemcpy((void*)(qdevbufferfrom+nbufqfrom*i),(void*)&(qs[offset*i]), N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)(qdevbufferto+nbufqto*i),(void*)&(qe[offset*i]), N*sizeof(float), cudaMemcpyHostToDevice);
  }

  cudaMemcpy((void*)testposdev,(void*)testpos.data(), M*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy((void*)testnumdev,(void*)testnum.data(), M*sizeof(int), cudaMemcpyHostToDevice);


  kernel_indicator2<ndof><<<GRID,BLOCK>>>(robotdev,polydatadev,qdevbufferfrom,nbufqfrom,qdevbufferto,nbufqto,resdevbuffer,resdevbufferext,testposdev,testnumdev,N, numthreads);


  cudaMemcpy((void*)res,(void*)resdevbuffer,N*sizeof(int), cudaMemcpyDeviceToHost);

}

//! same paircheck as above, but with compressed storage:
//! checks pairs: (qs[i],...) ->  (qe(posqe[i]),...) , ...., (qe[posqe[i]+numqe[i]-1],...) for i=0,...,M-1
template<int ndof>
int RobotConfigspace<ndof>::indicator2(const float* qs, const int M, const float* qe, int *res, const int *posqe, const int *numqe, const int N, const int offset)
{

  msg("this indicator2 is not implemented!");
  return -1;

}


//! boundary check on CPU

//! structure like indicator function
//! returns if lies in boundaries
template<int ndof>
int RobotConfigspace<ndof>::check_boundaries(const float* q, int* res, int offset, int N)
{
  for(int ix=0;ix<N;++ix,++iy){
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
    if(qi>=mins[i] || qi<maxs[i]){return 1;}
  }
  return 0;
}

















