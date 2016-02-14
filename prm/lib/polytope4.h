#ifndef POLYTOPE4_H
#define POLYTOPE4_H

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
