#include <iostream>

#define CUDA_IMPLEMENTATION
#include <cuda.h>


#include "lib/geo4.h"
#include "lib/kinematics.h"
#include "lib/collision4.h"

#include "lib/util.hpp"

using namespace std;
using namespace geo4;

const float pi=3.14159265358;


__global__ void kernel(float* q, float* res, int n){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i<n){
      float4 u=make_float4(1.0,0.0,0.0);
      trafo4 T(0.0,0.0,q[i],0.0);
      T.apply(u);
      res[i]=u.y;
  }
}

void handle(cudaError_t res){
  if(res!=0){
    cout<<"cudaError_t:  "<<cudaGetErrorString(res)<<endl;
  }
}

int main()
{

  int device = 0;
  cudaError_t res=cudaSetDevice(device);
  if(res!=0){
    cout<<"Error setting device: "<<cudaGetErrorString(res)<<endl;
    return -1;
  }

  cudaDeviceProp p;
  cudaGetDeviceProperties(&p, device);
  cout << "Device: " << p.name << endl;
  cout << "MP: " << p.multiProcessorCount << endl;
  cout << "Compute: " << p.major << "." << p.minor << endl;
  cout << "maxThreadsPerBlock: " << p.maxThreadsPerBlock << endl;
  cout << "maxThreadsDim: " << p.maxThreadsDim <<endl;
  cout << "maxGridSize: " << p.maxGridSize << endl;
  cout << "warpSize: " << p.warpSize << endl;
  cout << "totalGlobalMem: " << p.totalGlobalMem/1024/1024 << "MB" << endl;

  int n = 32, BLOCK = 32, GRID = (n + BLOCK - 1)/BLOCK;

  printvar(n);
  printvar(GRID);
  printvar(BLOCK);


  float *hst_q=new float[n]();
  float *dev_q;
  float *hst_res=new float[n]();
  float *dev_res;

  for(int i=0;i<n;++i){
    hst_q[i]=pi/4.0+i*pi/4.0/n;
  }

  res=cudaMalloc((void**)&dev_q, n * sizeof(float));
  handle(res);
  res=cudaMalloc((void**)&dev_res, n * sizeof(float));
  handle(res);

  cout<<"allocated"<<endl;

  res=cudaMemcpy(dev_q, hst_q, n * sizeof(float), cudaMemcpyHostToDevice);
  handle(res);

  kernel<<<GRID, BLOCK>>>(dev_q,dev_res,n);

  res=cudaMemcpy(hst_res, dev_res, n * sizeof(float), cudaMemcpyDeviceToHost);
  handle(res);

  cudaFree(dev_q);
  cudaFree(dev_res);

  printarr(hst_res,n);

  return 0;
}

