#ifndef COLLISION4_H
#define COLLISION4_H

#include "cuda_head.h"

#include "geo4.h"
#include "polytope4.h"

#include "util.hpp"


namespace collision4{

  using namespace geo4;


  #define max_for_loop 10
  #define max_for_loop_whole 10
  #define max_vertices_number 100
  #define collision_eps 0.000001


  class support_vertex_searcher{
  public:
    qualifierd support_vertex_searcher(const polytope4& P, int* vmarks_buffer){
      counter=0;
      vmarks=vmarks_buffer;
      //dprintvard(P.n);
      for(int i=0;i<P.n;++i){
          vmarks[i]=0;
      }
    }

    qualifierd int search_support_vertex(const polytope4& P, int p, float4& Sp){
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
    qualifierd orthtrafo23(float4 rk){
      float x=rk.x;
      float y=rk.y;
      float z=rk.z;
      float d2=x*x+y*y;
      float d=sqrt_(d2);
      M[0]=x*z/d; M[1]=y*z/d; M[2]=-d;
      M[3]=-y/d;  M[4]=x/d;   M[5]=0.0;
      //M[6]=x;     M[7]=y;     M[8]=z;
    }
    qualifierd void mult23(const float4& x, float2& res){
      res.x=M[0]*x.x+M[1]*x.y+M[2]*x.z;
      res.y=M[3]*x.x+M[4]*x.y+M[5]*x.z;
    }

    qualifierd void mult23T(const float2& x, float4& res){
      res.x=M[0]*x.x+M[3]*x.y;
      res.y=M[1]*x.x+M[4]*x.y;
      res.z=M[2]*x.x+M[5]*x.y;
    }
  };

  qualifierd bool find_half_plane(const float4* R, int k, float4& w, const float4& rk){
    dmsg("find_half_plane: -----------");
    df4print(rk);
    df4print(w);

    float dp=sprod(w,rk);

    dprintvarf(dp);

    if(dp>0){
        dmsg("nothing to do");
        return true;
    }
    orthtrafo23 M(rk);
    float2 ra,rb,T;

    df4print(R[0]);
    df4print(R[1]);

    M.mult23(R[0],ra);
    M.mult23(R[1],rb);

    if(cross2(ra,rb)<0){
      float2 temp=ra; //optimieren!
      ra=rb;
      rb=temp;
      dmsg("swapping ra, rb");
    }

    df2print(ra);
    df2print(rb);

    float4 R0test;
    M.mult23T(ra,R0test);
    df4print(R0test);

    for(int i=2;i<k;++i){
      M.mult23(R[i],T);
      dprintvard(i);
      df4print(R[i]);
      df2print(T);
      if(cross2(ra,T)>0){
          if(cross2(rb,T)>0){
            rb=T;
            df2print(rb);
          }else{
            dmsg("T in < ra,rb");
          }
      }else{
          if(cross2(rb,T)<0){
            ra=T;
            df2print(ra);
          }else{
            dmsg("exit");
            return false;
          }
      }
    }


    dmsg("result:");

    normalize(ra);
    normalize(rb);
    add(ra,rb,T);

    df2print(ra);
    df2print(rb);
    df2print(T);

    normalize(T);  //Achtung: kann nan rauskommen
    M.mult23T(T,w);

    df2print(T);
    df4print(w);
    dprintvarf(sprod(w,w));

    for(int i=0;i<=k;++i){
        dprintvarf(sprod(w,R[i]));
    }

    dmsg("end");

    return true;
  }


  //!0 = no collision, 1 = collision, 2 = max iterations reached
  qualifierd int seperating_vector_algorithm(const polytope4& P, const polytope4& Q, const trafo4& tp, const trafo4& tq){
    int p=0;
    int q=0;
    float4 S;
    //sub(tq.translation(),tp.translation(),S);
    float4 Pq, Pp;
    tq.apply(Q.vertices[0],Pq);
    tp.apply(P.vertices[0],Pp);
    sub(Pq,Pp,S);

    df4print(tp.translation());
    df4print(tq.translation());
    df4print(S);

    float norm=normalize(S);
    if(norm<collision_eps){
      S.x=1.0;
    }

    float4 Sp,Sq,p0,q0,w;
    float4 R[max_for_loop];
    int combsp[max_for_loop];
    int combsq[max_for_loop];
#if 0
    int* vmarks_buffer_P=new int[P.n];
    int* vmarks_buffer_Q=new int[Q.n];
#else
    int vmarks_buffer_P[max_vertices_number];
    int vmarks_buffer_Q[max_vertices_number];
#endif

    support_vertex_searcher psearcher(P,&vmarks_buffer_P[0]);
    support_vertex_searcher qsearcher(Q,&vmarks_buffer_Q[0]);


    for(int k=0, l=0; k<max_for_loop && l<max_for_loop_whole; ++l ){
        float4 *rk=&R[k];

        /*hostonly(printvar(k);)*/
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

#if 1
          dmsg("+++++++++++++++++++++++++++++++++++++++++");
          dprintvard(k);
          df4print(S);
          dprintvarf(sprod(S,S));
          dprintvard(p);
          dprintvard(q);
          df4print(p0);
          df4print(q0);
          df4print(R[k]);
          dprintvarf(dp);
          dmsg("");
#endif

        if(dp>=-collision_eps){
          //save S, p, q
          return 0;
        }
        combsp[k]=p;
        combsq[k]=q;
        bool repetition=false;
        for(int j=0;j<k;++j)if(combsp[j]==combsp[k] && combsq[j]==combsq[k]){
            dmsg("repetition");
            repetition=true;
            break;
        }
        if(repetition){ //wenn (p,q) bereits aufgetreten
          S=w;
        }else{
          if(k==1){
            add(R[0],*rk,w); //w=<r0+r1>
            normalize(w);
          }
          if(k>=2 && find_half_plane(&R[0],k,w,*rk)==false){
            dmsg("no half plane");
            return l+1;
          }
          lin(S,-2.0*dp,*rk,S);
          ++k;
        }

        df4print(w);

    }
//dmsg("Hello");
    //hostonly(msg("maximal iterations reached!"));
    return -1;
  }






}


#endif // COLLISION4_H
