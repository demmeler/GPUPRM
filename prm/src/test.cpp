#include <iostream>
#include <time.h>

using namespace std;

//for parsing parameters
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

const int N=1000;
__constant__ float carray[N];


__global__ void set_kernel(float *array, const float *src, int n){
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i<n){
        array[i]=src[i];
    }
}



int main(int argc, char** argv)
{


  //! Parameters

  int num=32;

  //! Parse inputs
  static struct option long_options[] = {
      {"num",      required_argument,        0,  'n' },
  };

  char opt= 0;
  int long_index =0;
  while ((opt = getopt_long(argc, argv,"",
                            long_options, &long_index )) != -1) {
      switch (opt) {
      case 'n' : num = atoi(optarg);
      default:
          cout<<"bad argument: "<< opt <<endl;
          exit(EXIT_FAILURE);
      }
  }

  int n=100;
  float *dev_array;
  float *host_array=new float[n];
  float src_array[N];

  cudaMalloc((void**)&dev_array, n*sizeof(float));

  float *carray_ptr;
  cudaGetSymbolAddress((void**)&carray_ptr, carray[0]);

  for(int i=0;i<50;++i){
    src_array[i]=0.9;
  }
  cudaMemcpyToSymbol(carray, (void*)src_array, 65*sizeof(float), 0, cudaMemcpyHostToDevice);
  cout << cudaGetErrorString(cudaGetLastError())<<endl;

  int BLOCK=256, GRID=(n + BLOCK - 1)/BLOCK;

  set_kernel<<<GRID,BLOCK>>>(dev_array,carray_ptr,n);


  cudaMemcpy((void*)host_array, (void*)dev_array, n*sizeof(float), cudaMemcpyDeviceToHost);

  cout<<"num: "<<num<<endl;
  for(int i=0;i<n;++i){
  cout<<i<<" "<<host_array[i]<<" "<<src_array[i]<<endl;
  }
}
