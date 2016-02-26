#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


///   **************************
///   *     general print      *
///   *    implementations     *
///   **************************

#define check(x) if(!(x)){ std::cout<<"check "<<#x<<" failed"<<std::endl;}
#define printvar(x) std::cout<<#x<<"="<<x<<std::endl;
#define printarr(x,n) std::cout<<#x<<"=";for(int i=0;i<n;++i){std::cout<<x[i]<<" ";}std::cout<<std::endl;
#define msg(x) std::cout<<x<<std::endl;


#define tick(t1) std::chrono::high_resolution_clock::time_point t1=std::chrono::high_resolution_clock::now();
#define tock(t1) {std::chrono::high_resolution_clock::time_point endtime___=std::chrono::high_resolution_clock::now(); \
                 int time=std::chrono::duration_cast<std::chrono::milliseconds>(endtime___-t1).count();                \
                 std::cout << "Since " << #t1 << ": \t" << time << "ms" << endl; }
#define timedelta(t1,t2) {int time=std::chrono::duration_cast<std::chrono::milliseconds>(endtime___-t1).count(); \
                          std::cout << "Since " << #t1 << ": \t" << time << "ms" << endl;}



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
                             printf("\n");
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


void write_file(std::string path, double* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(double));
  outfile.close();
}

void write_file(std::string path, float* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(float));
  outfile.close();
}

void write_file(std::string path, int* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(int));
  outfile.close();
}

void read_file(std::string path, double* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(double));
  infile.close();
}

void read_file(std::string path, float* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(float));
  infile.close();
}

void read_file(std::string path, int* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(int));
  infile.close();
}






#endif // UTIL_HPP
