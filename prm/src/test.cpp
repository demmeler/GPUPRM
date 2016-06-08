#include <iostream>
#include <time.h>

using namespace std;

//for parsing parameters
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>



__global__ void set_kernel(float *array, float val, int n){
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i<n){
        array[i]=val;
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

  cudaMalloc((void**)&dev_array, n*sizeof(float));
  cudaMemcpy((void*)dev_array, (void*)host_array, n*sizeof(float), cudaMemcpyHostToDevice);


  int BLOCK=256, GRID=(n + BLOCK - 1)/BLOCK;

  set_kernel<<<GRID,BLOCK>>>(dev_array,0,n);


  cudaMemcpy((void*)host_array, (void*)dev_array, n*sizeof(float), cudaMemcpyDeviceToHost);


  cout<<host_array[n/2]<<endl;

}
