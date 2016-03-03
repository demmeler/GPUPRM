#include <iostream>

#define SILENT
#define CUDA_IMPLEMENTATION

#include "lib/cuda_head.h"

#include "lib/geo4.h"
#include "lib/collision4.h"
#include "lib/polytope4data.h"
#include "lib/kinematics.h"

using namespace collision4;

const float pi=3.14159265358;

#ifdef CUDA_IMPLEMENTATION
template<int ndof>
__global__ void kernel(Robot<ndof>* robot, polytope4data* polydata, float* q, int* coll, int n){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i<n){
    dmsg("Hello");
    coll[i]=0;
    Kinematics<ndof> kin(robot);
    kin.calculate(&q[i],n);
    dmsg("Hello");
    for(int dof0=0;dof0<=ndof;++dof0)
    for(int dof1=dof0+1;dof1<=ndof;++dof1){
      dmsg("Hello");
      dprintvard(dof0);
      dprintvard(dof1);
      int numsys0=polydata->get_numsys(dof0), numsys1=polydata->get_numsys(dof1);
      dprintvard(numsys0);
      dprintvard(numsys1);
      dmsg("Hello4");
      dmsg("Hello4");
      dmsg("Hello4");
      dmsg("Hello4");
      for(int k0=0;k0<numsys0;++k0)
      for(int k1=0;k1<numsys1;++k1){
        dmsg("Hello5");
        polytope4 poly0, poly1;
        polydata->get_polytope(poly0, dof0,k0);
        polydata->get_polytope(poly1, dof1,k1);
#if 1
        dp4print(poly0);
        dp4print(poly1);
#endif
        int res=seperating_vector_algorithm(poly0,poly1,kin.trafos[dof0],kin.trafos[dof1]);
        if(res!=0){
          coll[i]=res;
        }
      }
    }
  }//if
}
#endif


template<int ndof>
void kernel_(Robot<ndof>* robot, polytope4data* polydata, float* q, int* coll, int n){
  for(int i=0;i<n;++i){
    if(i<n){
      coll[i]=0;
      Kinematics<ndof> kin(robot);
      kin.calculate(&q[i],n);
      for(int dof0=0;dof0<=ndof;++dof0) for(int dof1=dof0+1;dof1<=ndof;++dof1){
        int numsys0=polydata->get_numsys(dof0), numsys1=polydata->get_numsys(dof1);
        for(int k0=0;k0<numsys0;++k0) for(int k1=0;k1<numsys1;++k1){
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
}


void build_example1(Robot<1>* &robot, polytope4data* &polydata){
  robot=new Robot<1>;

  robot->a[0]=0.9;
  robot->alpha[0]=0.0;
  robot->q[0]=0.0;
  robot->d[0]=0.3;
  robot->types[0]=rotational;



  polytope4 P[3];
  int sys[3];

  trafo4 t0(0.0, 0.0, 0.0, 0.0);
  generate_simplex(P[0], 1.0, 1.0, 1.0);
  transform(P[0],t0);
  sys[0]=0;

  trafo4 t1(0.0, 0.0, 0.0, 0.0);
  generate_simplex(P[1], 1.0, 1.0, 1.0);
  transform(P[1],t1);
  sys[1]=1;

  //trafo4 t2(0.0, -pi/2.0, 0.0, 0.0);
  //generate_quader(P[2], 1.0, 1.0, 1.0);
  //transform(P[2],t2);



  polydata=new polytope4data;
  polydata->build(&P[0],&sys[0],2,1);

}


int main()
{

  const int ndof=1;

  Robot<ndof> *robot;
  polytope4data *polydata;

  build_example1(robot, polydata);


  polytope4 poly;
  polydata->get_polytope(poly,1,0);
  polytope4 poly0;
  polydata->get_polytope(poly0,0,0);

  Kinematics<ndof> kin(robot);

  float qtest[ndof];
  qtest[0]=pi/2.0*47/48;
  kin.calculate(&qtest[0],1);

#if 0
  p4print(poly,kin.trafos[1]);
  t4print(kin.trafos[1]);

  printvar(poly.n);
  printvar(poly.m);

  printarr(poly.dsp, poly.n);
  printarr(poly.cnt, poly.n);
  printarr(poly.dest,poly.m);

  printvar(seperating_vector_algorithm(poly0,poly,kin.trafos[0],kin.trafos[1]));
#endif

#if 1
  int n=1024*1024;
  float *q;
  q=new float[n];
  float qmin=-0.0,qmax=pi/2.0;
  for(int i=0;i<n;++i){
    q[i]=qmin+i*(qmax-qmin)/n;
  }

  int *coll;
  coll=new int[n];

#if 1
  int BLOCK = 256, GRID = (n + BLOCK - 1)/BLOCK;
  float *qdev;
  int *colldev;
  cudaMalloc((void**)&qdev,n*sizeof(float));
  cudaMalloc((void**)&colldev,n*sizeof(int));

  cudaMemcpy((void*)qdev, (void*)q, n*sizeof(float), cudaMemcpyHostToDevice);
  Robot<ndof> *robotdev;
  cudaMalloc((void**)&robotdev, sizeof(Robot<ndof>));
  cudaMemcpy((void*)robotdev,(void*)robot, sizeof(Robot<ndof>), cudaMemcpyHostToDevice);

  polytope4data *polydatadev_hostref=new polytope4data;
  copy_host_to_device(*polydatadev_hostref,*polydata);
  polytope4data *polydatadev;
  cudaMalloc((void**)&polydatadev, sizeof(polytope4data));
  cudaMemcpy((void*)polydatadev, (void*)polydatadev_hostref, sizeof(polytope4data), cudaMemcpyHostToDevice);

  kernel<ndof><<<GRID, BLOCK>>>(robotdev, polydatadev, qdev, colldev, n);

  cudaMemcpy((void*)coll, (void*)colldev, n*sizeof(int), cudaMemcpyDeviceToHost);
#else
  kernel_(robot,polydata, q, coll, n);
#endif

  //printarr(q,n);
  //printarr(coll,n);
#endif


  msg("finished");
  return 0;
}

