#include <iostream>

#include <cuda.h>

#define CUDA_IMPLEMENTATION
#include "lib/geo4.h"


using namespace std;
using namespace geo4;

const float pi=3.14159265358;


__global__ void kernel(float* q, float* res, int n){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i<n){
      res[i]=q[i];
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

  int n = 256 * 1024, BLOCK = 256, GRID = (n + BLOCK - 1)/BLOCK;


  float *hst_q=new float[n]();
  float *dev_q;
  float *hst_res=new float[n]();
  float *dev_res;

  for(int i=0;i<n;++i){
    hst_q[i]=pi/2;
  }

  cudaMalloc((void**)&dev_q, GRID * sizeof(float));
  cudaMalloc((void**)&dev_res, GRID * sizeof(float));

  cudaMemcpy(dev_q, hst_q, n * sizeof(float), cudaMemcpyHostToDevice);

  kernel<<<GRID, BLOCK>>>(dev_q,dev_res,n);

  cudaMemcpy(dev_res, hst_res, n * sizeof(float), cudaMemcpyHostToDevice);

  cout<<hst_res[100]<<endl;

  cout<<"test.cpp"<<endl;
  return 0;
}

