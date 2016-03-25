#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


const float pi=3.1415926535897;


///   **************************
///   *     general print      *
///   *    implementations     *
///   **************************

#ifndef MPI_CODE
    #define check(x) if(!(x)){ std::cout<<"check "<<#x<<" failed"<<std::endl;}
    #define printvar(x) std::cout<<#x<<"="<<x<<std::endl;
    #define printarr(x,n) std::cout<<#x<<"=";for(int index___=0;index___<n;++index___){std::cout<<x[index___]<<" ";}std::cout<<std::endl;
    #define msg(x) std::cout<<x<<std::endl;
#else
    #include <mpi.h>
    #define check(x) if(!(x)){int rank=0; MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
                              std::stringstream stream__; stream__<<rank<<": check "<<#x<<" failed"<<std::endl;  \
                              std::cout<<stream__.str(); }
    #define printvar(x) {int rank=0; MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
                         std::stringstream stream__; stream__<<rank<<": "<<#x<<"="<<x<<std::endl; \
                         std::cout<<stream__.str(); }
    #define printarr(x,n) {int rank=0; MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
                           std::stringstream stream__; \
                           stream__<<rank<<": "<<#x<<"="; \
                           for(int index___=0;index___<n;++index___){stream__<<x[index___]<<" ";} \
                           stream__<<std::endl; \
                           std::cout<<stream__.str();}
    #define msg(x) {int rank=0; MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
                    std::stringstream stream__; stream__<<rank<<": "<<x<<std::endl; \
                    std::cout<<stream__.str();}
#endif

#if 1//def CUDA_IMPLEMENTATION
  #include <assert.h>
  #define cudaassert(x) {if(x!=cudaSuccess){printvar(x); printvar(cudaGetErrorString(x));}  assert(x==cudaSuccess);}
#endif


///   **************************
///   *    device printing     *
///   *    implementations     *
///   **************************

#ifndef SILENT
#ifdef CUDA_IMPLEMENTATION
#if 0
    #define dprintsetup std::stringstream str; str<<blockIdx.x<<" "<<threadIdx.x;
    #define dcheck(x) if(!(x)) printf( "check %s failed.\n", #x );
    #define dprintvar(x) { dprintsetup str<<#x<<"="<<x<<"\n"; printf("%s",str.str()); }
    #define dprintarr(x,n) { dprintsetup str<<#x<<"=";                         \
                             for(int i=0;i<n;++i){str<<x[i]<<" ";} \
                             str<<"\n"; printf("%s",str.str()); }
    #define dmsg(x) { dprintsetup printf("%s %s\n",str.str(),x ); }
#endif
#if 0
    #define dcheck(v) if(!(v)) printf( "%d %d: check %s failed.\n", blockIdx.x, threadIdx.x, #v);
    #define dprintvarf(v) printf( "%d %d: %s=%f\n", blockIdx.x, threadIdx.x, #v, v);
    #define dprintvard(v) printf( "%d %d: %s=%d\n", blockIdx.x, threadIdx.x, #v, v);
    #define dprintarr(x,n) { printf("%s=",#x);                              \
                             for(int i=0;i<n;++i){printf("%d ",x[i]);}      \
                             printf("\n");
    #define dmsg(s) printf("%d %d: %s\n", blockIdx.x, threadIdx.x, s);
#endif
#if 1
    #define dcheck(v) if(!(v)) printf( "check %s failed.\n",  #v);
    #define dprintvarf(v) printf( "%s=%f\n",  #v, v);
    #define dprintvard(v) printf( "%s=%d\n",  #v, v);
    #define dprintarr(x,n) { printf("%s=",#x);                              \
                             for(int i=0;i<n;++i){printf("%d ",x[i]);}      \
                             printf("\n"); }
    #define dmsg(x) printf("%s\n",  x);
#endif

#else
    #define dcheck(x) if(!(x)){ std::cout<<"check "<<#x<<" failed"<<std::endl; }
    #define dprintvard(x) std::cout<<#x<<"="<<x<<std::endl;
    #define dprintvarf(x) std::cout<<#x<<"="<<x<<std::endl;
    #define dprintarr(x,n) std::cout<<#x<<"=";for(int i=0;i<n;++i){std::cout<<x[i]<<" ";}std::cout<<std::endl;
    #define dmsg(x) std::cout<<x<<std::endl;
#endif
#else
    #define dcheck(x)
    #define dprintvard(x)
    #define dprintvarf(x)
    #define dprintarr(x,n)
    #define dmsg(x)
#endif





///   **************************
///   *      file writing      *
///   *    implementations     *
///   **************************


inline void write_file(std::string path, const double* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(double));
  outfile.close();
}

inline void write_file(std::string path, const float* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(float));
  outfile.close();
}

inline void write_file(std::string path, const int* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(int));
  outfile.close();
}

inline void read_file(std::string path, double* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(double));
  infile.close();
}

inline void read_file(std::string path, float* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(float));
  infile.close();
}

inline void read_file(std::string path, int* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(int));
  infile.close();
}






#endif // UTIL_HPP
