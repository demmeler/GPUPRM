#ifndef COLLISION_HPP
#define COLLISION_HPP

#include "util.hpp"
#include "trafo.hpp"
#include <vector>
#include <math.h>

namespace collision{

  const int L=10;

  struct Polytope{
    double* vertices;
    int n;
    //!edges saved in crs format
    int* dsp;
    int* cnt;
    int* dest;
    int m;
  };

  inline void normalize(double* x){
    double sum=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    if(sum>0){
      x[0]/=sum;
      x[1]/=sum;
      x[2]/=sum;
    }
  }

  inline void normalize2(double* x){
    double sum=sqrt(x[0]*x[0]+x[1]*x[1]);
    if(sum>0){
      x[0]/=sum;
      x[1]/=sum;
    }
  }

  inline void add(double* x, double* y, double* res){
    res[0]=x[0]+y[0];
    res[1]=x[1]+y[1];
    res[2]=x[2]+y[2];
  }

  inline void add2(double* x, double* y, double* res){
    res[0]=x[0]+y[0];
    res[1]=x[1]+y[1];
  }

  inline void sub(double* x, double* y, double* res){
    res[0]=x[0]-y[0];
    res[1]=x[1]-y[1];
    res[2]=x[2]-y[2];
  }

  inline double dot_prod(double* x, double* y){
    return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
  }

  inline void copy(double* from, double* to){
    to[0]=from[0];
    to[1]=from[1];
    to[2]=from[2];
  }

  inline void copy2(double* from, double* to){
    to[0]=from[0];
    to[1]=from[1];
    to[2]=from[2];
  }

  inline void scale(double* x, double f){
    x[0]*=f; x[1]*=f; x[2]*=f;
  }

  inline void mult33(double* M, double* x, double* res){
    double x0=x[0];
    double x1=x[1];
    double x2=x[2];
    res[0]=M[0]*x0+M[1]*x1+M[2]*x2;
    res[1]=M[3]*x0+M[4]*x1+M[5]*x2;
    res[2]=M[6]*x0+M[7]*x1+M[8]*x2;
  }

  inline void mult23(double* M, double* x, double* res){
    double x0=x[0];
    double x1=x[1];
    double x2=x[2];
    res[0]=M[0]*x0+M[1]*x1+M[2]*x2;
    res[1]=M[3]*x0+M[4]*x1+M[5]*x2;
  }

  inline void mult23T(double* M, double* x, double* res){
    double x0=x[0];
    double x1=x[1];
    res[0]=M[0]*x0+M[3]*x1;
    res[1]=M[1]*x0+M[4]*x1;
    res[2]=M[2]*x0+M[5]*x1;
  }

  inline double cross2(double* x, double* y){
    return x[0]*y[1]-x[1]*y[0];
  }

  class support_vertex_searcher{
  public:
    support_vertex_searcher(const Polytope& P){
      counter=0;
      vmarks=new int[P.n]();
    }

    inline int search_support_vertex(const Polytope& P, int p, double* Sp){
      int newp=p;
      double max=dot_prod(P.vertices+(3*p),Sp);
      counter++;
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
            double dp=dot_prod(P.vertices+(3*v),Sp);
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

  inline void compute_M(double* M, double* rk){
    double x=rk[0];
    double y=rk[1];
    double z=rk[2];
    double d2=x*x+y*y;
    double d=sqrt(d2);
    M[0]=x*z/d; M[1]=y*z/d; M[2]=-d;
    M[3]=-y/d;  M[4]=x/d;   M[5]=0.0;
    M[6]=x;     M[7]=y;     M[8]=z;
  }

  inline bool find_half_plane(double* R,int k,double* w, double* rk, double* M){
    double dp=dot_prod(w,rk);
    if(dp>0)return true;
    compute_M(M,rk);
    double* ra=new double[2];
    double* rb=new double[2];
    double* T=new double[2];

#define delall delete ra,rb,T;

    mult23(M,R,ra);
    mult23(M,R+3,rb);
    if(cross2(ra,rb)<0){
      double* temp=ra;
      ra=rb;
      rb=temp;
    }
    for(int i=6;i<=3*k;i+=3){
      mult23(M,R+i,T);
      if(cross2(ra,T)>0){
          if(cross2(rb,T)>0)
            copy2(T,rb);
      }else{
          if(cross2(rb,T)<0){
            copy2(T,ra);
          }else{
            delall;
            return false;
          }
      }
    }
    add2(ra,rb,T);
    normalize2(T);
    mult23T(M,T,w);

    delall;
    return true;
#undef delall
  }


  inline bool seperating_vector_algorithm(const Polytope& P, const Polytope& Q, const trafo& tp, const trafo& tq){
    int p=0;
    int q=0;
    double* S=new double[3];
    const double* Tp=tp.get_translation();
    const double* Tq=tq.get_translation();
    for(int i=0;i<3;++i)S[i]=Tq[i]-Tp[i];
    normalize(S);
    //TODO: get better starting values

    double* Sp=new double[3];
    double* Sq=new double[3];
    double* R =new double[3*L];
    double* rk=R;
    double* p0=new double[3];
    double* q0=new double[3];
    double* w=new double[3];

    double* M=new double[9];

#define delall delete S,Sp,Sq,R,p0,q0;

    support_vertex_searcher psearcher(P);
    support_vertex_searcher qsearcher(Q);

    msg("start");
    for(int k=0;k<L;++k){
        printvar(k);
        tp.apply_rot_inv(S,Sp);
        scale(S,-1.0);
        tq.apply_rot_inv(S,Sq);
        scale(S,-1.0);
        p=psearcher.search_support_vertex(P,p,Sp);
        q=qsearcher.search_support_vertex(Q,q,Sq);
        tp.apply(P.vertices+(3*p),p0);
        tq.apply(Q.vertices+(3*q),q0);
        sub(q0,p0,rk); normalize(rk);

        printarr(S,3);
        printvar(p);
        printvar(q);

        printarr(p0,3);
        printarr(q0,3);
        printarr(rk,3);

        double dp=dot_prod(S,rk);

        printvar(dp);

        if(dp>=0.0){
          //save S, p, q
          delall;
          return false;
        }
        if(false){ //wenn (p,q) bereits aufgetreten
          copy(w,S);
        }else{
          if(k==1){
            add(R,rk,w); //w=<r0+r1>
            normalize(w);
          }
          if(find_half_plane(R,k,w,rk,M)==false){
            msg("no half plane");
            delall;
            return true;
          }
        }

        for(int i=0;i<3;++i){
          S[i]=S[i]-2*dp*rk[i];
        }

        rk+=3;
    }
    delall;
    return true;
#undef delall
  }






}

#endif // COLLISION_HPP
