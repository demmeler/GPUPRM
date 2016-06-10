#ifndef POLYTOPE4_H
#define POLYTOPE4_H

#include "cuda_head.h"
#include "util.hpp"

#include "geo4.h"
#include <sstream>

//using namespace geo4;

namespace collision4{

  struct polytope4{
    float4* restrict vertices;
    int n;
    //!edges saved in crs format
    int* restrict dsp;
    int* restrict cnt;
    int* restrict dest;
    int m;
  };


  ///   **************************
  ///   *       polytope4        *
  ///   *         memory         *
  ///   *    implementations     *
  ///   **************************
#if 0
  //!allocate polytope4 struct on device
  int cuda_alloc(polytope4** pdevpoly, int n, int m){
    hostonly(
      return -1;
    )
    cudaonly(
      polytope4 polytemp;
      polytemp.n=n;
      polytemp.m=m;
      int res[6];
      res[0]=cudaMalloc((void**)&(polytemp.vertices), n * sizeof(float4));
      res[1]=cudaMalloc((void**)&(polytemp.dsp), n * sizeof(int));
      res[2]=cudaMalloc((void**)&(polytemp.cnt), n * sizeof(int));
      res[3]=cudaMalloc((void**)&(polytemp.dest), m * sizeof(int));
      res[4]=cudaMalloc((void**)pdevpoly,sizeof(polytope4));
      res[5]=cudaMemcpy((void*)*pdevpoly,&polytemp, sizeof(polytope4), cudaMemcpyHostToDevice);
      for(int i=0;i<6;++i)if(res[i]!=0)return res[i];
      return 0;
    )
  }

  //!copy hostpoly to (initialized) devpoly
  int cuda_copy_to_device(const polytope4* hostpoly, polytope4* devpoly){
    hostonly(
      return -1;
    )
    cudaonly(
      int res[5];
      polytope4 polytemp;
      res[0]=cudaMemcpy((void*)&polytemp,(void*)devpoly, sizeof(polytope4), cudaMemcpyDeviceToHost);
      res[1]=cudaMemcpy((void*)polytemp.vertices, (void*)hostpoly->vertices, hostpoly->n * sizeof(float4), cudaMemcpyHostToDevice);
      res[2]=cudaMemcpy((void*)polytemp.dsp, (void*)hostpoly->dsp, hostpoly->n * sizeof(int), cudaMemcpyHostToDevice);
      res[3]=cudaMemcpy((void*)polytemp.cnt, (void*)hostpoly->cnt, hostpoly->n * sizeof(int), cudaMemcpyHostToDevice);
      res[4]=cudaMemcpy((void*)polytemp.dest, (void*)hostpoly->dest, hostpoly->m * sizeof(int), cudaMemcpyHostToDevice);
      for(int i=0;i<5;++i)if(res[i]!=0)return res[i];
      return 0;
    )
  }

  //!alloc (uninitialized) devpoly and copy from host
  int cuda_init_and_copy(const polytope4* hostpoly, polytope4** pdevpoly){
    int res[2];
    res[0]=cuda_alloc(pdevpoly, hostpoly->n, hostpoly->m);
    res[1]=cuda_copy_to_device(hostpoly, *pdevpoly);
    for(int i=0;i<2;++i)if(res[i]!=0)return res[i];
    return 0;
  }

  //!copy devpoly to (initialized) hostpoly
  int cuda_copy_to_host(polytope4* devpoly, const polytope4* hostpoly){
    hostonly(
      return -1;
    )
    cudaonly(
      int res[5];
      polytope4 polytemp;
      res[0]=cudaMemcpy(&polytemp,(void*)devpoly, sizeof(polytope4), cudaMemcpyDeviceToHost);
      res[1]=cudaMemcpy((void*)hostpoly->vertices, (void*)polytemp.vertices, polytemp.n * sizeof(float4), cudaMemcpyDeviceToHost);
      res[2]=cudaMemcpy((void*)hostpoly->dsp, (void*)polytemp.dsp, polytemp.n * sizeof(int), cudaMemcpyDeviceToHost);
      res[3]=cudaMemcpy((void*)hostpoly->cnt, (void*)polytemp.cnt, polytemp.n * sizeof(int), cudaMemcpyDeviceToHost);
      res[4]=cudaMemcpy((void*)hostpoly->dest, (void*)polytemp.dest, polytemp.m * sizeof(int), cudaMemcpyDeviceToHost);
      for(int i=0;i<5;++i)if(res[i]!=0)return res[i];
      return 0;
    )
  }

  //!cudaFree devpoly arrays
  int cuda_free(polytope4* devpoly){
    hostonly(
      return -1;
    )
    cudaonly(
      int res[6];
      polytope4 polytemp;
      res[0]=cudaMemcpy(&polytemp,(void*)devpoly, sizeof(polytope4), cudaMemcpyDeviceToHost);
      res[1]=cudaFree(polytemp.vertices);
      res[2]=cudaFree(polytemp.dsp);
      res[3]=cudaFree(polytemp.cnt);
      res[4]=cudaFree(polytemp.dest);
      res[5]=cudaFree(devpoly);
      for(int i=0;i<6;++i)if(res[i]!=0)return res[i];
      return 0;
    )
  }
#endif

  //!delete hostpoly arrays
  inline int host_free(polytope4& hostpoly){
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

  inline void transform(polytope4& P, geo4::trafo4& t){
    for(int i=0;i<P.n;++i){
      t.apply(P.vertices[i]);
    }
  }

#if 0
  //!create simpley (for polytope on host)
  inline void generate_simplex(polytope4& P, float lx, float ly, float lz){
    P.n=4;
    P.vertices=new float4[P.n];
    P.vertices[0]=geo4::make_float4(0.0, 0.0, 0.0);
    P.vertices[1]=geo4::make_float4(lx, 0.0, 0.0);
    P.vertices[2]=geo4::make_float4(0.0, ly, 0.0);
    P.vertices[3]=geo4::make_float4(0.0, 0.0, lz);


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
  inline void generate_quader(polytope4& P, float lx, float ly, float lz){
    P.n=8;
    P.vertices=new float4[P.n];
    P.vertices[0]=geo4::make_float4(0.0, 0.0, 0.0);
    P.vertices[1]=geo4::make_float4(lx, 0.0, 0.0);
    P.vertices[2]=geo4::make_float4(0.0, ly, 0.0);
    P.vertices[3]=geo4::make_float4(lx, ly, 0.0);
    P.vertices[4]=geo4::make_float4(0.0, 0.0, lz);
    P.vertices[5]=geo4::make_float4(lx, 0.0, lz);
    P.vertices[6]=geo4::make_float4(0.0, ly, lz);
    P.vertices[7]=geo4::make_float4(lx, ly, lz);

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
#endif

  ///   **************************
  ///   *        output          *
  ///   *    implementations     *
  ///   **************************

  inline void print(const polytope4& p, const geo4::trafo4& t, std::ostream& out, const std::string& name=""){
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
#ifdef SILENT
    #define dp4print(P)
#else
    #define dp4print(P) {dmsg(#P); for(int i=0;i<P.n;++i) df4print(P.vertices[i]);}
#endif
}

#endif // POLYTOPE4_H
