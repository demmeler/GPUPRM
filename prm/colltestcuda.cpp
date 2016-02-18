#include <iostream>

#define CUDA_IMPLEMENTATION
#include "lib/cuda_head.h"

#include "lib/geo4.h"
#include "lib/collision4.h"

using namespace std;
using namespace geo4;
using namespace collision4;

const float pi=3.14159265358;

__global__ void kernel(polytope4* Pdev, polytope4* Qdev, float* q, int* coll, int n){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
#if 1
  dprintvard(i);
#endif
  if(i<n){
    //dprintvard(Pdev->n);

    trafo4 tp(0.0, 0.0, 0.0, 0.0);
    trafo4 tq(-0.8, 0.0, q[i], 0.0);
    coll[i]=seperating_vector_algorithm(*Pdev,*Qdev,tp,tq);
  }
}


void kernel_(polytope4* Pdev, polytope4* Qdev, float* q, int* coll, int n){
  for(int i=0;i<n;++i){
    if(i<n){
      trafo4 tp(0.0, 0.0, 0.0, 0.0);
      trafo4 tq(-0.8, 0.0, q[i], 0.0);
      coll[i]=seperating_vector_algorithm(*Pdev,*Qdev,tp,tq);
    }
  }
}

int main()
{

  trafo4 tp(0.0, 0.0, 0.0, 0.0);
  trafo4 tq(1.0, 0.0, 0.0, -0.01);
  trafo4 tr(-0.9, 0.0, 0.0005, 0.0);


  polytope4 P, *Pdev;
  generate_simplex(P, 1.0, 1.0, 1.0);
  cuda_init_and_copy(&P,&Pdev);

  polytope4 Q, *Qdev;
  generate_simplex(Q, 1.0, 1.0, 1.0);
  cuda_init_and_copy(&Q,&Qdev);

  polytope4 R, *Rdev;
  generate_quader(R, 1.0, 1.0, 1.0);
  cuda_init_and_copy(&R,&Rdev);

#if 0
  polytope4 P1;
  generate_simplex(P1, 1.0, 1.0, 1.0);
  p4print(P1,tp);
  cuda_copy_to_host(Pdev,P1);
  p4print(P1,tp);
#endif


#if 1
  int n=24, BLOCK = 4, GRID = (n + BLOCK - 1)/BLOCK;
  float *q, *qdev;
  q=new float[n];
  float qmin=-pi/2.0,qmax=pi/2.0;
  for(int i=0;i<n;++i){
    q[i]=qmin+i*(qmax-qmin)/n;
  }

  int *coll, *colldev;
  coll=new int[n];

#if 1
  cudaMalloc((void**)&qdev,n*sizeof(float));
  cudaMalloc((void**)&colldev,n*sizeof(int));

  cudaMemcpy((void*)qdev, (void*)q, n*sizeof(float), cudaMemcpyHostToDevice);
  kernel<<<GRID, BLOCK>>>(Pdev, Qdev, qdev, colldev, n);
  cudaMemcpy((void*)coll, (void*)colldev, n*sizeof(int), cudaMemcpyDeviceToHost);
#else
  kernel_(&P,&Q, q, coll, n);
#endif

  printarr(q,n);
  printarr(coll,n);
#endif

  cuda_free(Pdev);
  cuda_free(Qdev);
  cuda_free(Rdev);

  msg("finished");
  return 0;
}

