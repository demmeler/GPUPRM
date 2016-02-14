#ifndef COLLISION4_H
#define COLLISION4_H

#undef qualifier
#ifdef CUDA_IMPLEMENTATION
  #include <cuda.h>
  #define qualifier __device__
  #define cudaonly(x) x
  #define devonly(x)
#else
  #define qualifier inline
  #define cudaonly(x)
  #define devonly(x) x
#endif

#include "geo4.h"
#include "robot.h"
#include "polytope4.h"

#include "util.hpp"
#include <vector>
#include <math.h>

namespace collision4{

  const int L=10;


  class support_vertex_searcher{
  public:
    qualifier support_vertex_searcher(int* vmarks_buffer){
      counter=0;
      vmarks=vmarks_buffer;
    }

    qualifier int search_support_vertex(const polytope4& P, int p, float4& Sp){
      int newp=p;
      float max=sprod(P.vertices[p],Sp);
      ++counter;
      vmarks[p]=counter;
      //msg("start supp vert");
      //printarr(Sp,3);
      //printvar(max);
      do{
        //printvar(p);
        p=newp;
        int imax=P.dsp[p]+P.cnt[p];
        for(int i=P.dsp[p];i<imax;++i){
          //printvar(i);
          int v=P.dest[i];
          if(vmarks[v]!=counter){
            //printvar(v);
            float dp=sprod(P.vertices[v],Sp);
            //printvar(dp);
            if(dp>max){
              max=dp;
              newp=v;
            }
            vmarks[v]=counter;
          }
        }
      }while(newp!=p);
      //printvar(p);
      //msg("end supp vert");
      return p;
    }

  private:
    int* vmarks;
    int counter;
  };


  class orthtrafo23{
  public:
    float M[6];
    qualifier orthtrafo23(float4 rk){
      float x=rk.x;
      float y=rk.y;
      float z=rk.z;
      float d2=x*x+y*y;
      float d=sqrt(d2);
      M[0]=x*z/d; M[1]=y*z/d; M[2]=-d;
      M[3]=-y/d;  M[4]=x/d;   M[5]=0.0;
      //M[6]=x;     M[7]=y;     M[8]=z;
    }
    qualifier void mult23(const float4& x, float2& res){
      res.x=M[0]*x.x+M[1]*x.y+M[2]*x.z;
      res.y=M[3]*x.x+M[4]*x.y+M[5]*x.z;
    }

    qualifier void mult23T(const float2& x, float4& res){
      res.x=M[0]*x.x+M[3]*x.y;
      res.y=M[1]*x.x+M[4]*x.y;
      res.z=M[2]*x.x+M[5]*x.y;
    }
  };

  qualifier bool find_half_plane(const float4* R, int k, float4& w, const float4& rk){
    double dp=sprod(w,rk);
    if(dp>0)return true;
    orthtrafo23 M(rk);
    float2 ra,rb,T;

    M.mult23(R[0],ra);
    M.mult23(R[1],rb);
    if(cross2(ra,rb)<0){
      float2 temp=ra; //optimieren!
      ra=rb;
      rb=temp;
    }
    for(int i=2;i<=k;++i){
      M.mult23(R[i],T);
      if(cross2(ra,T)>0){
          if(cross2(rb,T)>0)
            rb=T;
      }else{
          if(cross2(rb,T)<0){
            ra=T;
          }else{
            return false;
          }
      }
    }
    add(ra,rb,T);
    normalize(T);
    M.mult23T(T,w);

    return true;
  }


  //!0 = no collision, 1 = collision, 2 = max iterations reached
  qualifier int seperating_vector_algorithm(const polytope4& P, const polytope4& Q, const trafo4& tp, const trafo4& tq){
    int p=0;
    int q=0;
    float4 S;
    sub(tq.translation(),tp.translation(),S);
    normalize(S);
    //TODO: get better starting values

    float4 Sp,Sq,p0,q0,w;
    float4* rk;
    float4 R[L];
    rk=&R[0];

    int* vmarks_buffer_P=new int[P.n];
    int* vmarks_buffer_Q=new int[Q.n];

    support_vertex_searcher psearcher(&vmarks_buffer_P[0]);
    support_vertex_searcher qsearcher(&vmarks_buffer_Q[0]);

    for(int k=0;k<L;++k){
        /*devonly(printvar(k);)*/

        tp.apply_rot_inv(S,Sp);
        S*=-1.0;
        tq.apply_rot_inv(S,Sq);
        S*=-1.0;
        p=psearcher.search_support_vertex(P,p,Sp);
        q=qsearcher.search_support_vertex(Q,q,Sq);
        tp.apply(P.vertices[p],p0);
        tq.apply(Q.vertices[q],q0);
        sub(q0,p0,*rk); normalize(*rk);

        float dp=sprod(S,*rk);

        /*devonly(
          f4print(S);
          printvar(p);
          printvar(q);
          f4print(p0);
          f4print(q0);
          f4print(*rk);
          printvar(dp);
        )*/

        if(dp>=0.0){
          //save S, p, q
          return 0;
        }
        if(false){ //wenn (p,q) bereits aufgetreten
          S=w;
        }else{
          if(k==1){
            add(R[0],*rk,w); //w=<r0+r1>
            normalize(w);
          }
          if(find_half_plane(R,k,w,*rk)==false){
            //devonly(msg("no half plane");)
            return 1;
          }
        }

        lin(S,-2.0*dp,*rk,S);

        ++rk;
    }
    //devonly(msg("maximal iterations reached!"));
    return 2;
  }






}


#endif // COLLISION4_H
