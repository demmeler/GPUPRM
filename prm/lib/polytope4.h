#ifndef POLYTOPE4_H
#define POLYTOPE4_H

#include "cuda_head.h"


#include "geo4.h"
#include <sstream>

using namespace geo4;

namespace collision4{

  struct polytope4{
    float4* vertices;
    int n; //Reihenfolge?
    //!edges saved in crs format
    int* dsp;
    int* cnt;
    int* dest;
    int m;
  };


  ///   **************************
  ///   *       polytope4         *
  ///   *         memory         *
  ///   *    implementations     *
  ///   **************************

  //!allocate devpoly according to hostpoly
  int cuda_alloc(polytope4& devpoly, int n, int m){
    hostonly(
      return -1;
    )
    cudaonly(
      devpoly.n=n;
      devpoly.m=m;
      int res[4];
      res[0]=cudaMalloc((void**)&devpoly.vertices, devpoly.n * sizeof(float4));
      res[1]=cudaMalloc((void**)&devpoly.dsp, devpoly.n * sizeof(int));
      res[2]=cudaMalloc((void**)&devpoly.cnt, devpoly.n * sizeof(int));
      res[3]=cudaMalloc((void**)&devpoly.dest, devpoly.m * sizeof(int));
      for(int i=0;i<4;++i)if(res[i]!=0)return res[i];
      return 0;
    )
  }

  //!copy hostpoly to (initialized) devpoly
  int cuda_copy_to_device(const polytope4& hostpoly, polytope4& devpoly){
    hostonly(
      return -1;
    )
    cudaonly(
      int res[4];
      res[0]=cudaMemcpy((void*)devpoly.vertices, (void*)hostpoly.vertices, hostpoly.n * sizeof(float4), cudaMemcpyHostToDevice);
      res[1]=cudaMemcpy((void*)devpoly.dsp, (void*)hostpoly.dsp, hostpoly.n * sizeof(int), cudaMemcpyHostToDevice);
      res[2]=cudaMemcpy((void*)devpoly.cnt, (void*)hostpoly.cnt, hostpoly.n * sizeof(int), cudaMemcpyHostToDevice);
      res[3]=cudaMemcpy((void*)devpoly.dest, (void*)hostpoly.dest, hostpoly.m * sizeof(int), cudaMemcpyHostToDevice);
      for(int i=0;i<4;++i)if(res[i]!=0)return res[i];
      return 0;
    )
  }

  //!alloc (uninitialized) devpoly and copy from host
  int cuda_init_and_copy(const polytope4& hostpoly, polytope4& devpoly){
    int res[2];
    res[0]=cuda_alloc(devpoly, hostpoly.n, hostpoly.m);
    res[1]=cuda_copy_to_device(hostpoly, devpoly);
    for(int i=0;i<2;++i)if(res[i]!=0)return res[i];
    return 0;
  }

  //!copy devpoly to (initialized) hostpoly
  int cuda_copy_to_host(polytope4& devpoly, const polytope4& hostpoly){
    hostonly(
      return -1;
    )
    cudaonly(
      int res[4];
      res[0]=cudaMemcpy((void*)hostpoly.vertices, (void*)devpoly.vertices, devpoly.n * sizeof(float4), cudaMemcpyDeviceToHost);
      res[1]=cudaMemcpy((void*)hostpoly.dsp, (void*)devpoly.dsp, devpoly.n * sizeof(int), cudaMemcpyDeviceToHost);
      res[2]=cudaMemcpy((void*)hostpoly.cnt, (void*)devpoly.cnt, devpoly.n * sizeof(int), cudaMemcpyDeviceToHost);
      res[3]=cudaMemcpy((void*)hostpoly.dest, (void*)devpoly.dest, devpoly.m * sizeof(int), cudaMemcpyDeviceToHost);
      for(int i=0;i<4;++i)if(res[i]!=0)return res[i];
      return 0;
    )
  }

  //!cudaFree devpoly arrays
  int cuda_free(polytope4& devpoly){
    hostonly(
      return -1;
    )
    cudaonly(
      int res[4];
      res[0]=cudaFree(devpoly.vertices);
      res[1]=cudaFree(devpoly.dsp);
      res[2]=cudaFree(devpoly.cnt);
      res[3]=cudaFree(devpoly.dest);
      for(int i=0;i<4;++i)if(res[i]!=0)return res[i];
      return 0;
    )
  }

  //!delete hostpoly arrays
  int host_free(polytope4& hostpoly){
    delete hostpoly.vertices;
    delete hostpoly.dsp;
    delete hostpoly.cnt;
    delete hostpoly.dest;
    return 0;
  }


  ///   **************************
  ///   *       polytope4        *
  ///   *        example         *
  ///   *    implementations     *
  ///   **************************

  //!create simpley (for polytope on host)
  void generate_simplex(polytope4& P, float lx, float ly, float lz){
    P.n=4;
    P.vertices=new float4[P.n];
    P.vertices[0]=make_float4(0.0, 0.0, 0.0);
    P.vertices[1]=make_float4(lx, 0.0, 0.0);
    P.vertices[2]=make_float4(0.0, ly, 0.0);
    P.vertices[3]=make_float4(0.0, 0.0, lz);


    P.m=12;
    P.dsp=new int[P.n];
    P.cnt=new int[P.n];
    P.dest=new int[P.m];

    for(int i=0;i<P.n;++i){
      P.dsp[i]=3*i;
      P.cnt[i]=3;
    }

    P.dest[0]=1;
    P.dest[1]=2;
    P.dest[2]=3;
    P.dest[3]=0;
    P.dest[4]=2;
    P.dest[5]=3;
    P.dest[6]=0;
    P.dest[7]=1;
    P.dest[8]=3;
    P.dest[9]=0;
    P.dest[10]=1;
    P.dest[11]=2;
  }

  //!create quader (for polytope on host)
  void generate_quader(polytope4& P, float lx, float ly, float lz){
    P.n=8;
    P.vertices=new float4[P.n];
    P.vertices[0]=make_float4(0.0, 0.0, 0.0);
    P.vertices[1]=make_float4(lx, 0.0, 0.0);
    P.vertices[2]=make_float4(0.0, ly, 0.0);
    P.vertices[3]=make_float4(lx, ly, 0.0);
    P.vertices[4]=make_float4(0.0, 0.0, lz);
    P.vertices[5]=make_float4(lx, 0.0, lz);
    P.vertices[6]=make_float4(0.0, ly, lz);
    P.vertices[7]=make_float4(lx, ly, lz);

    P.m=24;
    P.dsp=new int[P.n];
    P.cnt=new int[P.n];
    P.dest=new int[P.m];

    for(int i=0;i<P.n;++i){
      P.dsp[i]=3*i;
      P.cnt[i]=3;
    }

    P.dest[0]=1;
    P.dest[1]=3;
    P.dest[2]=4;
    P.dest[3]=0;
    P.dest[4]=2;
    P.dest[5]=5;
    P.dest[6]=1;
    P.dest[7]=3;
    P.dest[8]=6;
    P.dest[9]=0;
    P.dest[10]=2;
    P.dest[11]=7;
    P.dest[12]=5;
    P.dest[13]=7;
    P.dest[14]=0;
    P.dest[15]=4;
    P.dest[16]=6;
    P.dest[17]=1;
    P.dest[18]=5;
    P.dest[19]=7;
    P.dest[20]=2;
    P.dest[21]=4;
    P.dest[22]=6;
    P.dest[23]=3;
  }


  ///   **************************
  ///   *        output          *
  ///   *    implementations     *
  ///   **************************

  void print(const polytope4& p, const geo4::trafo4& t, std::ostream& out, const std::string& name=""){
    out<<"polytope4: "<<name<<std::endl;
    for(int i=0;i<p.n;++i){
      float4 x;
      t.apply(p.vertices[i],x);
      std::stringstream namei;
      namei<<name<<".vertices["<<i<<"]";
      geo4::print(x,out,namei.str());
    }
  }

#define p4print(p,t) print(p,t,std::cout,#p);
}

#endif // POLYTOPE4_H
