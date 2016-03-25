
#include "config.h"

#include "robotconfigspace.h"


#include "cuda_head.h"
#include "util.hpp"

#include "collision4.h"
#include "polytope4data.h"
#include "kinematics.h"

#include <vector>
#include <cmath>


template<int ndof>
void RobotConfigspace<ndof>::construct(const Robot<ndof>* robot_,
                 const polytope *polys_,
                 const int* sys_,
                 const int N_,
                 const int *from_, const int *to_, const int M_,
                 const float* mins_, const float* maxs_, const float dq_,
                 const int nbuf_){//, const int numthreadsmax_){
  robot=robot_;
  kin=new Kinematics<ndof>(robot);
  for(int i=0;i<ndof;++i){
    mins[i]=mins_[i];
    maxs[i]=maxs_[i];
  }
  dq=dq_;

  polylist.polys=polys_;
  polylist.sys=sys_;
  polylist.N=N_;

  polylist.from=from_;
  polylist.to=to_;
  polylist.M=M_;

  nbufres=nbuf_;
  nbufqto=ndof*nbuf_;
  nbufqfrom=ndof*nbuf_;
  nbuftest=nbuf_;


  devloaded=false;

  requeststack_id=0;
}

//! robot_: Robot Object with Denavit-Hartenberg data
//! polydata_: Object containing compressed Polytope data of robot arm and surrounding
//! mins[ndof]: maximal q values
//! maxs[ndof]: minimal   -"-
//! dq:         edge resolution (stepsize)
template<int ndof>
RobotConfigspace<ndof>::RobotConfigspace(const Robot<ndof>* robot_,
                                         const polytope* polys_,  const int* sys_, const int N_,
                                         const int *from_, const int *to_, const int M_,
                                         const float* mins_, const float* maxs_, const float dq_,
                                         const int nbuf_)
{
    construct(robot_, polys_, sys_, N_, from_, to_, M_, mins_, maxs_, dq_, nbuf_);
}

template<int ndof>
RobotConfigspace<ndof>::RobotConfigspace(std::string path, const float dq_, const int nbuf_)
{
    polytope *polys_;
    int *sys_, *from_, *to_;
    int N_,M_;
    Robot<ndof>* robot_;
    float *mins_, *maxs_;
    load_config(path,robot_,polys_,sys_,N_,from_,to_,M_, mins_, maxs_,false);

    construct(robot_, polys_, sys_, N_, from_, to_, M_, mins_, maxs_, dq_, nbuf_);

    delete mins_, maxs_;
}

template<int ndof>
int RobotConfigspace<ndof>::load_config(std::string path, Robot<ndof>* &robot, polytope* &polys,
                int* &sys, int &N, int* &from, int* &to, int& M, float *&mins, float *&maxs, bool printmsg) const{

    //! DH params

  int ndof_loaded=0;
  read_file(path+"/ndof.bin",&ndof_loaded,1);
  if(ndof_loaded!=ndof){
    msg("error: wrong ndof");
    return -1;
  }

  robot=new Robot<ndof>;

  read_file(path+"/dh/a.bin",&(robot->a[0]),ndof);
  read_file(path+"/dh/alpha.bin",&(robot->alpha[0]),ndof);
  read_file(path+"/dh/q.bin",&(robot->q[0]),ndof);
  read_file(path+"/dh/d.bin",&(robot->d[0]),ndof);
  read_file(path+"/dh/types.bin",(int*)&(robot->types[0]),ndof);

  if(printmsg){
      printvar(ndof);
      printarr(robot->a,ndof);
      printarr(robot->alpha,ndof);
      printarr(robot->q,ndof);
      printarr(robot->d,ndof);
      printarr(robot->types,ndof);
  }



  mins=new float[ndof];
  maxs=new float[ndof];
  /*for(int i=0;i<ndof;++i){
    mins[i]=-pi; maxs[i]=1.5*pi;
  }*/
  read_file(path+"/mins.bin",mins,ndof);
  read_file(path+"/maxs.bin",maxs,ndof);

  printvar(mins[0]+pi);
  printvar(maxs[0]-1.5*pi);

    //! Polytopes

  read_file(path+"/polys/N.bin",&N,1);
  if(N<=0){
    msg("error: N<=0");
    return -1;
  }

  sys=new int[N];
  read_file(path+"/polys/sys.bin",sys,N);

  read_file(path+"/pairs/M.bin",&M, 1);   //--> check machen ob file existiert
  check(M>0);
  from=new int[M];
  to=new int[M];
  read_file(path+"/pairs/from.bin", from, M);
  read_file(path+"/pairs/to.bin", to, M);


  if(printmsg){
    printvar(M);
    printarr(from, M);
    printarr(to,M);
    printarr(sys,N);
  }

  polys=new polytope[N];
  for(int i=0;i<N;++i){
    std::stringstream polypathstream;
    polypathstream<<path<<"/polys/poly"<<i;
    std::string polypath=polypathstream.str();
    int size[2];
    read_file(polypath+"/size.bin",&(size[0]),2);
    int n,m;
    n=polys[i].n=size[0];
    m=polys[i].m=size[1];
    polys[i].dsp=new int[n];
    polys[i].cnt=new int[n];
    polys[i].dest=new int[m];
    polys[i].vertices=new polytope::vec4[4*n];
    float* vertices=new float[3*n];
    read_file(polypath+"/vertices.bin",vertices,3*n);
    for(int j=0;j<n;++j){
        polys[i].vertices[j].x=vertices[3*j];
        polys[i].vertices[j].y=vertices[3*j+1];
        polys[i].vertices[j].z=vertices[3*j+2];
        polys[i].vertices[j].w=1.0;
    }
    read_file(polypath+"/dsp.bin",polys[i].dsp,n);
    read_file(polypath+"/cnt.bin",polys[i].cnt,n);
    read_file(polypath+"/dest.bin",polys[i].dest,m);


    if(printmsg){
      msg("-----");
      printvar(i);
      //geo4::trafo4 t(0.0,0.0,0.0,0.0);
      //p4print(polys[i],t);

      printvar(n);
      printvar(m);
      printarr(polys[i].dsp,n);
      printarr(polys[i].cnt,n);
      printarr(polys[i].dest,m);
    }

  }
  return 0; //TODO error handling?
}


//!initialization function copy polytope and robot data to gpu etc..
template<int ndof>
int RobotConfigspace<ndof>::init(const int ressource_rank, const int ressource_size)
{
  if(devloaded){
    clear();
  }

  polydata=new collision4::polytope4data;
  polydata->build(polylist.polys,polylist.sys,polylist.N,ndof,polylist.from, polylist.to, polylist.M);

#ifdef CUDA_IMPLEMENTATION
  int devcount;
  cudaGetDeviceCount(&devcount);

  int device = (ressource_rank*devcount)/ressource_size;
  cudaSetDevice(device);

  cudaDeviceProp p;
  cudaGetDeviceProperties(&p, device);

  std::stringstream stream;
  stream<<"ressource_rank="<<ressource_rank<<"/"<<ressource_size<<" device="<<device<<"/"<<devcount<<std::endl;

  stream << "Device: " << p.name << std::endl;
  stream << "MP: " << p.multiProcessorCount << std::endl;
  stream << "Compute: " << p.major << "." << p.minor << std::endl;
  std::cout<<stream.str();

  polydatadev_hostref=new collision4::polytope4data;

  assert(0==collision4::copy_host_to_device(*polydatadev_hostref,*polydata,true));
  cudaassert(cudaMalloc((void**)&polydatadev, sizeof(collision4::polytope4data)));
  cudaassert(cudaMemcpy((void*)polydatadev, (void*)polydatadev_hostref, sizeof(collision4::polytope4data), cudaMemcpyHostToDevice));

  cudaassert(cudaMalloc((void**)&robotdev, sizeof(Robot<ndof>)));
  cudaassert(cudaMemcpy((void*)robotdev,(void*)robot, sizeof(Robot<ndof>), cudaMemcpyHostToDevice));

  cudaassert(cudaMalloc((void**)&qdevbufferfrom, ndof*nbufqfrom*sizeof(float)));
  cudaassert(cudaMalloc((void**)&qdevbufferto, ndof*nbufqto*sizeof(float)));
  cudaassert(cudaMalloc((void**)&testnumdev, nbuftest*sizeof(int)));
  cudaassert(cudaMalloc((void**)&testposdev, nbuftest*sizeof(int)));
  cudaassert(cudaMalloc((void**)&resdevbuffer, nbufres*sizeof(int)));

#else


#endif

  devloaded=true;

  return 0; //TODO: errorhandling?

}

template<int ndof>
int RobotConfigspace<ndof>::clear()
{
  // TODO
    msg("error: RobotConfigspace<ndof>::clear() not yet implemented");
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
  kin->calculate(q,1);

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

          /*dprintarr(dest,destnum);
          dprintvard(polydata->sys[k0]);
          dprintvard(polydata->sys[k1]);
          dt4print(kin->trafos[polydata->sys[k0]]);
          dt4print(kin->trafos[polydata->sys[k1]]);
          */

          int result=collision4::seperating_vector_algorithm(poly0,poly1,kin->trafos[polydata->sys[k0]],kin->trafos[polydata->sys[k1]]);
          if(result!=0){
            return result;
          }
      }
    }

  return 0;
}









//!
//! indicator functions on GPU
//!

//!help kernels

#ifdef CUDA_IMPLEMENTATION
template<class T>
__global__ void set_kernel(T *array, T val, int n){
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i<n){
        array[i]=val;
    }
}
#endif

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
                                  int* res,
                                  int* testpos, int* testnum,
                                  int N, int numthreads){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i<numthreads){
#else
void kernel_indicator2(const Robot<ndof>* robot,
                       const collision4::polytope4data* polydata,
                       const float* qs, int offsets,
                       const float* qe, int offsete,
                       int* res,
                       const int* testpos, const int* testnum,
                       int N, int numthreads){
  for(int i=0;i<numthreads;++i){
#endif

    //! determine the line in which the thread is involved
    int k; //! index of the line
           //! line k is handled by threads testpos[k], ...., testpos[k]+numpos[k]
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
    int resext=0;
    kin.calculate(&q[0],1);


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
            resext=result;
          }
      }
    }

    //printf("resext=%d\n",resext);

    //! reduce resext to res with ||
    if(resext!=0) res[k]=resext;

  }//if/for
}


//! checks if indicator function = 1 somewhere on the line between (qs[i],qs[i+N]) and (qe[i],qe[i+N]) for all i=1,...,N-1
//! res is return value
//! number of pairs
template<int ndof>
int RobotConfigspace<ndof>::indicator2_async(const float* qs, const float* qe, int *res, const int N, const int offset, int &request){
    assert(N<=nbuftest);
    assert(N<=nbufres);

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

  #ifdef CUDA_IMPLEMENTATION
    int BLOCK = 256, GRID = (numthreads + BLOCK - 1)/BLOCK;

    for(int i=0;i<ndof;++i){
        //pointer inkrement in cuda??
      cudaassert(cudaMemcpy((void*)(qdevbufferfrom+nbufqfrom*i),(void*)&(qs[offset*i]), N*sizeof(float), cudaMemcpyHostToDevice));
      cudaassert(cudaMemcpy((void*)(qdevbufferto+nbufqto*i),(void*)&(qe[offset*i]), N*sizeof(float), cudaMemcpyHostToDevice));
    }

    cudaassert(cudaMemcpy((void*)testposdev,(void*)testpos.data(), N*sizeof(int), cudaMemcpyHostToDevice));
    cudaassert(cudaMemcpy((void*)testnumdev,(void*)testnum.data(), N*sizeof(int), cudaMemcpyHostToDevice));


    //cudaassert(cudaMemcpy((void*)resdevbuffer,(void*)res,N*sizeof(int), cudaMemcpyHostToDevice));
    int GRIDN=(N + BLOCK - 1)/BLOCK;
    set_kernel<int><<<GRIDN,BLOCK>>>(resdevbuffer,0,N);

    kernel_indicator2<ndof><<<GRID,BLOCK>>>(robotdev,polydatadev,qdevbufferfrom,nbufqfrom,qdevbufferto,nbufqto,resdevbuffer,testposdev,testnumdev,N, numthreads);

  #else
    for(int k=0;k<N;++k){
      res[k]=0; //-> kernel?
    }

    kernel_indicator2<ndof>(robot,polydata,qs,offset,qe,offset,res,testpos.data(),testnum.data(),N, numthreads);

  #endif

    ++requeststack_id;
    request_data data;
    data.res=res;
    data.N=N;
    requeststack[requeststack_id]=data;
    request=requeststack_id;

    //printarr(res,N);
    //printvar(numthreads);

    return 0; //TODO error handling?
}

template<int ndof>
int RobotConfigspace<ndof>::indicator2_async_wait(int request){
    typename std::map<int,request_data>::iterator it;
    it=requeststack.find(request);
    assert(it!=requeststack.end());
    request_data data=requeststack[request];
    requeststack.erase(it);
#ifdef CUDA_IMPLEMENTATION
    cudaassert(cudaMemcpy((void*)data.res,(void*)resdevbuffer,data.N*sizeof(int), cudaMemcpyDeviceToHost));
#else
#endif

    return 0;
}

template<int ndof>
int RobotConfigspace<ndof>::indicator2(const float* qs, const float* qe, int *res, const int N, const int offset)
{
  int request;
  indicator2_async(qs, qe, res, N, offset, request);
  indicator2_async_wait(request);

  return 0; //TODO error handling?
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




//template class RobotConfigspace<2>;
//template class RobotConfigspace<3>;
template class RobotConfigspace<4>;











