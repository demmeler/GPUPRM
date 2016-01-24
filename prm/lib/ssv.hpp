#ifndef SSV_ELEMENT_HPP
#define SSV_ELEMENT_HPP

#include <math.h>
#include <float.h>
#include "trafo.hpp"
#include "matvec.hpp"
#include "util.hpp"
#include <vector>

namespace ssv{

  using namespace matvec;

  class ssvele;
  class ssvref;
  class ssvvec;

  class ssvele
  {
  public:
    ssvele(){}

    void set_from(const trafo& t,ssvele ref){
      l=ref.l;
      r2=ref.r2;
      t.apply(&ref.p,&p);
      t.apply_rot(&ref.d,&d);
    }

  private:
    double p[3];
    double d[3];
    double l;
    double r;
  };

  inline double dist(ssvele ele1, ssvele ele2);



  class ssvref{
    friend class ssvvec;
  public:
    ssvref(int type_, double* P_, double r_):type(type_),r(r_){
      int n=3*type;
      P=new double[n];
      for(int i=0;i<n;++i)P[i]=P_[i];
    }
    ssvref(const std::vector<double>& Pvec_, double r_):r(r_){
      type=Pvec_.size();
      check(0<=type && type<=3);
      int n=3*type;
      P=new double[n];
      double* P_=Pvec_.data();
      for(int i=0;i<n;++i)P[i]=P_[i];
    }
    ~ssvref(){delete P;}

  private:
    //!type = 1,2 or 3 = number of vertices
    int type;
    //!vertex array
    double* P;
    //!radius
    double r;
  };

  class ssvvec{
    ssvvec(int type_, double r_, int n_):type(type_),r(r_),n(n_){
      P=new double[3*type*n]();
    }
    //!set k-th ssv element from reference with trafo
    inline void set_from(const trafo& t, const ssvref& ref, const int k){
      int ind=3*k*type;
      t.apply(ref.P+ind,P+ind,type);
    }
    //!set k-th ssv element, ind==3*k*type
    inline void set_from_ind(const trafo& t, const ssvref& ref, const int ind){
      t.apply(ref.P+ind,P+ind,type);
    }

  private:
    //!type = 1,2 or 3 = number of vertices
    int type;
    //!vertex array
    double* P;
    //!radius
    double r;
    //!number of configurations
    int n;
  };

  //!set res[k]==true if the two ssvs collide in k-th configuration, else do nothing
  inline void test_collision(const ssvvec& v1, const ssvvec& v2, bool* res){
    if(v1.type==1){
        if(v2.type==1){
            test_collision11(v1,v2,res);
        }else if(v2.type==2){
            test_collision12(v1,v2,res);
        }else{//v2.type==3
            test_collision13(v1,v2,res);
        }
    }else if(v1.type==2){
        if(v2.type==1){
            test_collision12(v2,v1,res);
        }else if(v2.type==2){
            test_collision22(v1,v2,res);
        }else{//v2.type==3
            test_collision23(v1,v2,res);
        }
    }else{//v1.type==3
        if(v2.type==1){
            test_collision13(v2,v1,res);
        }else if(v2.type==2){
            test_collision23(v2,v1,res);
        }else{//v2.type==3
            test_collision33(v1,v2,res);
        }
    }
  }

  inline void test_collision(const ssvref& v1, const ssvvec& v2, bool* res){
    if(v1.type==1){
        if(v2.type==1){
            test_collision11(v1,v2,res);
        }else if(v2.type==2){
            test_collision12(v1,v2,res);
        }else{//v2.type==3
            test_collision13(v1,v2,res);
        }
    }else if(v1.type==2){
        if(v2.type==1){
            test_collision12(v2,v1,res);
        }else if(v2.type==2){
            test_collision22(v1,v2,res);
        }else{//v2.type==3
            test_collision23(v1,v2,res);
        }
    }else{//v1.type==3
        if(v2.type==1){
            test_collision13(v2,v1,res);
        }else if(v2.type==2){
            test_collision23(v2,v1,res);
        }else{//v2.type==3
            test_collision33(v1,v2,res);
        }
    }
  }

  inline void test_collision(const ssvvec& v1, const ssvref& v2, bool* res){
    test_collision(v2,v1,res);
  }



  //!concrete implementations of ssv collision tests
  //!***********************************************
  /*
   * Assumptions:
   *    v1.n==v2.n for vec-vec methods
   */

  //!point-point
  inline void test_collision11(const ssvvec& v1, const ssvvec& v2, bool* res){
    int l=0;
    for(int k=0;k<v1.n;++k){
        if(dist3(v1.P+l,v2.P+l)<v1.r+v2.r)res[k]==true;
        l+=3;
    }
  }

  inline void test_collision11(const ssvref& v1, const ssvvec& v2, bool* res){
    int l=0;
    for(int k=0;k<v1.n;++k){
        if(dist3(v1.P,v2.P+l)<v1.r+v2.r)res[k]==true;
        l+=3;
    }
  }

  inline void test_collision11(const ssvvec& v1, const ssvref& v2, bool* res){
    test_collision11(v2,v1,res);
  }

  //!point-line
  inline void test_collision12(const ssvvec& v1, const ssvvec& v2, bool* res){
    int l=0,m=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* bc=new double[3];
    for(int k=0;k<v1.n;++k){
      double* a=v2.P+m;
      double* b=v2.P+m+3;
      double* c=v1.P+l;
      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(c,b,bc);
      double e=dot_prod3(ac,ab);
      double f=dot_prod3(ab,ab);
      double dsq;
      if(e<0.0)dsq=normsq3(ac);
      else if(e>f)dsq=normsq3(bc);
      else dsq=normsq3(ac)-e*e/f;
      if(sqrt(dsq)<v1.r+v2.r)res[k]=true;
      l+=3;
      m+=6;
    }
    delete ab,ac,bc;
  }

  inline void test_collision12(const ssvref& v1, const ssvvec& v2, bool* res){
    int m=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* bc=new double[3];
    for(int k=0;k<v1.n;++k){
      double* a=v2.P+m;
      double* b=v2.P+m+3;
      double* c=v1.P;
      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(c,b,bc);
      double e=dot_prod3(ac,ab);
      double f=dot_prod3(ab,ab);
      double dsq;
      if(e<0.0)dsq=normsq3(ac);
      else if(e>f)dsq=normsq3(bc);
      else dsq=normsq3(ac)-e*e/f;
      if(sqrt(dsq)<v1.r+v2.r)res[k]=true;
      m+=6;
    }
    delete ab,ac,bc;
  }

  inline void test_collision12(const ssvvec& v1, const ssvref& v2, bool* res){
    int l=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* bc=new double[3];
    for(int k=0;k<v1.n;++k){
      double* a=v2.P;
      double* b=v2.P+3;
      double* c=v1.P+l;
      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(c,b,bc);
      double e=dot_prod3(ac,ab);
      double f=dot_prod3(ab,ab);
      double dsq;
      if(e<0.0)dsq=normsq3(ac);
      else if(e>f)dsq=normsq3(bc);
      else dsq=normsq3(ac)-e*e/f;
      if(sqrt(dsq)<v1.r+v2.r)res[k]=true;
      l+=3;
    }
    delete ab,ac,bc;
  }

  //!point triangle
  inline void test_collision13(const ssvvec& v1, const ssvvec& v2, bool* res){
    int l=0,m=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* ap=new double[3];
    double* bp=new double[3];
    double* cp=new double[3];
    double* q=new double[3]; //closest point
    for(int k=0;k<v1.n;++k){
      double* p=v1.P+l;
      double* a=v2.P+m;
      double* b=a+3;
      double* c=a+6;
      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(p,a,ap);
      sub3(p,b,bp);
      sub3(p,c,cp);
      double d1=dot_prod3(ab,ap);
      double d2=dot_prod3(ac,ap);
      double d3 = dot_prod3(ab, bp);
      double d4 = dot_prod3(ac, bp);
      double d5 = dot_prod3(ab, cp);
      double d6 = dot_prod3(ac, cp);

      double vc = d1*d4 - d3*d2;
      double vb = d5*d2 - d1*d6;
      double va = d3*d6 - d5*d4;
      if (d1 <= 0.0f && d2 <= 0.0f){
          //return a;
          copy3(a,q);
      }
      else if (d3 >= 0.0f && d4 <= d3){
          //return b;
          copy3(b,q);
      }
      else if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) { //edge AB
        double v = d1 / (d1 - d3);
        //return a + v * ab; // barycentric coordinates (1-v,v,0)
        mult3(v,ab,q);
        add3(a,q,q);
      }
      else if (d6 >= 0.0f && d5 <= d6){
          //return c;
          copy3(c,q);
      }
      else if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) { //edge AC
        double w = d2 / (d2 - d6);
        //return a + w * ac; // barycentric coordinates (1-w,0,w)
        mult3(w,ac,q);
        add(a,q,q);
      }
      else if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) { //edge BC
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b); // barycentric coordinates (0,1-w,w)
        sub(c,b,q);
        mult3(w,q,q);
        add(b,q,q);
      }
      else{ //inside face ABC
        double denom = 1.0f / (va + vb + vc);
        double v = vb * denom;
        double w = vc * denom;
        return a+ab*v+ac*w; //=u*a+v*b+w*c,u=va*denom=1.0f-v-w
        lincomb3(v,ab,w,ac,q);
        add(a,q,q);
      }

      if(dist(p,q)<v1.r+v2.r)res[k]=true;
      l+=3;
      m+=9;
    }
    delete ab,ac,ap,bp,cp,q;
  }

  inline void test_collision13(const ssvref& v1, const ssvvec& v2, bool* res){
    int m=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* ap=new double[3];
    double* bp=new double[3];
    double* cp=new double[3];
    double* q=new double[3]; //closest point
    for(int k=0;k<v2.n;++k){
      double* p=v1.P;
      double* a=v2.P+m;
      double* b=a+3;
      double* c=a+6;
      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(p,a,ap);
      sub3(p,b,bp);
      sub3(p,c,cp);
      double d1=dot_prod3(ab,ap);
      double d2=dot_prod3(ac,ap);
      double d3 = dot_prod3(ab, bp);
      double d4 = dot_prod3(ac, bp);
      double d5 = dot_prod3(ab, cp);
      double d6 = dot_prod3(ac, cp);

      double vc = d1*d4 - d3*d2;
      double vb = d5*d2 - d1*d6;
      double va = d3*d6 - d5*d4;
      if (d1 <= 0.0f && d2 <= 0.0f){
          //return a;
          copy3(a,q);
      }
      else if (d3 >= 0.0f && d4 <= d3){
          //return b;
          copy3(b,q);
      }
      else if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) { //edge AB
        double v = d1 / (d1 - d3);
        //return a + v * ab; // barycentric coordinates (1-v,v,0)
        mult3(v,ab,q);
        add3(a,q,q);
      }
      else if (d6 >= 0.0f && d5 <= d6){
          //return c;
          copy3(c,q);
      }
      else if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) { //edge AC
        double w = d2 / (d2 - d6);
        //return a + w * ac; // barycentric coordinates (1-w,0,w)
        mult3(w,ac,q);
        add(a,q,q);
      }
      else if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) { //edge BC
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b); // barycentric coordinates (0,1-w,w)
        sub(c,b,q);
        mult3(w,q,q);
        add(b,q,q);
      }
      else{ //inside face ABC
        double denom = 1.0f / (va + vb + vc);
        double v = vb * denom;
        double w = vc * denom;
        return a+ab*v+ac*w; //=u*a+v*b+w*c,u=va*denom=1.0f-v-w
        lincomb3(v,ab,w,ac,q);
        add(a,q,q);
      }

      if(dist(p,q)<v1.r+v2.r)res[k]=true;
      m+=9;
    }
    delete ab,ac,ap,bp,cp,q;
  }

  inline void test_collision13(const ssvvec& v1, const ssvref& v2, bool* res){
    int l=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* ap=new double[3];
    double* bp=new double[3];
    double* cp=new double[3];
    double* q=new double[3]; //closest point
    for(int k=0;k<v1.n;++k){
      double* p=v1.P+l;
      double* a=v2.P;
      double* b=a+3;
      double* c=a+6;
      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(p,a,ap);
      sub3(p,b,bp);
      sub3(p,c,cp);
      double d1=dot_prod3(ab,ap);
      double d2=dot_prod3(ac,ap);
      double d3 = dot_prod3(ab, bp);
      double d4 = dot_prod3(ac, bp);
      double d5 = dot_prod3(ab, cp);
      double d6 = dot_prod3(ac, cp);

      double vc = d1*d4 - d3*d2;
      double vb = d5*d2 - d1*d6;
      double va = d3*d6 - d5*d4;
      if (d1 <= 0.0f && d2 <= 0.0f){
          //return a;
          copy3(a,q);
      }
      else if (d3 >= 0.0f && d4 <= d3){
          //return b;
          copy3(b,q);
      }
      else if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) { //edge AB
        double v = d1 / (d1 - d3);
        //return a + v * ab; // barycentric coordinates (1-v,v,0)
        mult3(v,ab,q);
        add3(a,q,q);
      }
      else if (d6 >= 0.0f && d5 <= d6){
          //return c;
          copy3(c,q);
      }
      else if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) { //edge AC
        double w = d2 / (d2 - d6);
        //return a + w * ac; // barycentric coordinates (1-w,0,w)
        mult3(w,ac,q);
        add(a,q,q);
      }
      else if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) { //edge BC
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b); // barycentric coordinates (0,1-w,w)
        sub(c,b,q);
        mult3(w,q,q);
        add(b,q,q);
      }
      else{ //inside face ABC
        double denom = 1.0f / (va + vb + vc);
        double v = vb * denom;
        double w = vc * denom;
        return a+ab*v+ac*w; //=u*a+v*b+w*c,u=va*denom=1.0f-v-w
        lincomb3(v,ab,w,ac,q);
        add(a,q,q);
      }

      if(dist(p,q)<v1.r+v2.r)res[k]=true;
      l+=3;
    }
    delete ab,ac,ap,bp,cp,q;
  }

  inline double distance13(double* p, double* a, double* b, double* c, double* ab, double* ac, double* ap, double* bp, double* cp){
    //double* ab=new double[3];
    //double* ac=new double[3];
    //double* ap=new double[3];
    //double* bp=new double[3];
    //double* cp=new double[3];
    double* q=new double[3];

    //sub3(b,a,ab);
    //sub3(c,a,ac);
    //sub3(p,a,ap);
    //sub3(p,b,bp);
    //sub3(p,c,cp);

    double d1 = dot_prod3(ab,ap);
    double d2 = dot_prod3(ac,ap);
    double d3 = dot_prod3(ab, bp);
    double d4 = dot_prod3(ac, bp);
    double d5 = dot_prod3(ab, cp);
    double d6 = dot_prod3(ac, cp);

    double vc = d1*d4 - d3*d2;
    double vb = d5*d2 - d1*d6;
    double va = d3*d6 - d5*d4;
    if (d1 <= 0.0f && d2 <= 0.0f){
        //return a;
        copy3(a,q);
    }
    else if (d3 >= 0.0f && d4 <= d3){
        //return b;
        copy3(b,q);
    }
    else if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) { //edge AB
      double v = d1 / (d1 - d3);
      //return a + v * ab; // barycentric coordinates (1-v,v,0)
      mult3(v,ab,q);
      add3(a,q,q);
    }
    else if (d6 >= 0.0f && d5 <= d6){
        //return c;
        copy3(c,q);
    }
    else if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) { //edge AC
      double w = d2 / (d2 - d6);
      //return a + w * ac; // barycentric coordinates (1-w,0,w)
      mult3(w,ac,q);
      add(a,q,q);
    }
    else if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) { //edge BC
      double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      return b + w * (c - b); // barycentric coordinates (0,1-w,w)
      sub(c,b,q);
      mult3(w,q,q);
      add(b,q,q);
    }
    else{ //inside face ABC
      double denom = 1.0f / (va + vb + vc);
      double v = vb * denom;
      double w = vc * denom;
      return a+ab*v+ac*w; //=u*a+v*b+w*c,u=va*denom=1.0f-v-w
      lincomb3(v,ab,w,ac,q);
      add(a,q,q);
    }

    //delete ab,ac,ap,bp,cp,q;
    delete q;
    return dist(p,q);

  }

  inline double distance13mod(double* p, double* a, double* b, double* c, double* ab, double* ac, double* pa, double* pb, double* pc){
    //double* ab=new double[3];
    //double* ac=new double[3];
    //double* ap=new double[3];
    //double* bp=new double[3];
    //double* cp=new double[3];
    double* q=new double[3];

    //sub3(b,a,ab);
    //sub3(c,a,ac);
    //sub3(p,a,ap);
    //sub3(p,b,bp);
    //sub3(p,c,cp);

    double d1 = -dot_prod3(ab, pa);
    double d2 = -dot_prod3(ac, pa);
    double d3 = -dot_prod3(ab, pb);
    double d4 = -dot_prod3(ac, pb);
    double d5 = -dot_prod3(ab, pc);
    double d6 = -dot_prod3(ac, pc);

    double vc = d1*d4 - d3*d2;
    double vb = d5*d2 - d1*d6;
    double va = d3*d6 - d5*d4;
    if (d1 <= 0.0f && d2 <= 0.0f){
        //return a;
        copy3(a,q);
    }
    else if (d3 >= 0.0f && d4 <= d3){
        //return b;
        copy3(b,q);
    }
    else if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) { //edge AB
      double v = d1 / (d1 - d3);
      //return a + v * ab; // barycentric coordinates (1-v,v,0)
      mult3(v,ab,q);
      add3(a,q,q);
    }
    else if (d6 >= 0.0f && d5 <= d6){
        //return c;
        copy3(c,q);
    }
    else if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) { //edge AC
      double w = d2 / (d2 - d6);
      //return a + w * ac; // barycentric coordinates (1-w,0,w)
      mult3(w,ac,q);
      add(a,q,q);
    }
    else if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) { //edge BC
      double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      return b + w * (c - b); // barycentric coordinates (0,1-w,w)
      sub(c,b,q);
      mult3(w,q,q);
      add(b,q,q);
    }
    else{ //inside face ABC
      double denom = 1.0f / (va + vb + vc);
      double v = vb * denom;
      double w = vc * denom;
      return a+ab*v+ac*w; //=u*a+v*b+w*c,u=va*denom=1.0f-v-w
      lincomb3(v,ab,w,ac,q);
      add(a,q,q);
    }

    //delete ab,ac,ap,bp,cp,q;
    delete q;
    return dist(p,q);

  }

  inline void test_collision22(const ssvvec& v1, const ssvvec& v2, bool* res){
    int l=0,m=0;
    double* d1=new double[3];// = q1 - p1; // Direction vector of segment S1
    double* d2=new double[3];// = q2 - p2; // Direction vector of segment S2
    double* r=new double[3];// = p1 - p2;
    for(int k=0;k<v1.n;++k){
      double* p1=v1.P+l;
      double* q1=p1+3;
      double* p2=v2.P+m;
      double* q2=p2+3;
      sub3(q1,p1,d1);
      sub3(q2,p2,d2);
      sub3(p1,p2,r);
      double a = dot_prod3(d1, d1); // Squared length of segment S1, always nonnegative
      double e = dot_prod3(d2, d2); // Squared length of segment S2, always nonnegative
      double f = dot_prod3(d2, r);
      double b = dot_prod3(d1, d2);
      double denom = a*e-b*b; // Always nonnegative
      // If segments not parallel, compute closest point on L1 to L2 and // clamp to segment S1. Else pick arbitrary s (here 0)
      double t,s;
      if (denom != 0.0) {
        s = clamp((b*f - c*e) / denom, 0.0, 1.0);
      } else s = 0.0;
      // Compute point on L2 closest to S1(s) using
      // t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e t = (b*s + f) / e;
      // If t in [0,1] done. Else clamp t, recompute s for the new value
      // of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a // and clamp s to [0, 1]
      if (t < 0.0) {
        t = 0.0;
        s = clamp(-c / a, 0.0, 1.0);
      } else if (t > 1.0) {
        t = 1.0;
        s = clamp((b - c) / a, 0.0, 1.0);
      }

      //delta=(r+s*d1-t*d2)(r+s*d1-t*d2)
      double dsq=normsq3(r)+s*s*a+t*t*e+2.0*(s*dot_prod3(d1,r)+t*f+s*t*b);
      if(sqrt(dsq)<v1.r+v2.r)res[k]=true;

      l+=6;
      m+=6;
    }

    delete d1,d2,r;
  }

  inline void test_collision22(const ssvref& v1, const ssvvec& v2, bool* res){
    int m=0;
    double* d1=new double[3];// = q1 - p1; // Direction vector of segment S1
    double* d2=new double[3];// = q2 - p2; // Direction vector of segment S2
    double* r=new double[3];// = p1 - p2;
    for(int k=0;k<v1.n;++k){
      double* p1=v1.P;
      double* q1=p1+3;
      double* p2=v2.P+m;
      double* q2=p2+3;
      sub3(q1,p1,d1);
      sub3(q2,p2,d2);
      sub3(p1,p2,r);
      double a = dot_prod3(d1, d1); // Squared length of segment S1, always nonnegative
      double e = dot_prod3(d2, d2); // Squared length of segment S2, always nonnegative
      double f = dot_prod3(d2, r);
      double b = dot_prod3(d1, d2);
      double denom = a*e-b*b; // Always nonnegative
      // If segments not parallel, compute closest point on L1 to L2 and // clamp to segment S1. Else pick arbitrary s (here 0)
      double t,s;
      if (denom != 0.0) {
        s = clamp((b*f - c*e) / denom, 0.0, 1.0);
      } else s = 0.0;
      // Compute point on L2 closest to S1(s) using
      // t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e t = (b*s + f) / e;
      // If t in [0,1] done. Else clamp t, recompute s for the new value
      // of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a // and clamp s to [0, 1]
      if (t < 0.0) {
        t = 0.0;
        s = clamp(-c / a, 0.0, 1.0);
      } else if (t > 1.0) {
        t = 1.0;
        s = clamp((b - c) / a, 0.0, 1.0);
      }

      //delta=(r+s*d1-t*d2)(r+s*d1-t*d2)
      double dsq=normsq3(r)+s*s*a+t*t*e+2.0*(s*dot_prod3(d1,r)+t*f+s*t*b);
      if(sqrt(dsq)<v1.r+v2.r)res[k]=true;

      m+=6;
    }

    delete d1,d2,r;
  }

  inline void test_collision22(const ssvvec& v1, const ssvref& v2, bool* res){
    test_collision22(v2,v1,res);
  }


  inline double distance22(double* p1, double* q1, double* p2, double* q2, double* d1, double* d2, double* r){
    /*double* d1=new double[3];
    double* d2=new double[3];
    double* r=new double[3];
    sub3(q1,p1,d1);
    sub3(q2,p2,d2);
    sub3(p1,p2,r);*/
    double a = dot_prod3(d1, d1); // Squared length of segment S1, always nonnegative
    double e = dot_prod3(d2, d2); // Squared length of segment S2, always nonnegative
    double f = dot_prod3(d2, r);
    double b = dot_prod3(d1, d2);
    double denom = a*e-b*b; // Always nonnegative
    // If segments not parallel, compute closest point on L1 to L2 and // clamp to segment S1. Else pick arbitrary s (here 0)
    double t,s;
    if (denom != 0.0) {
      s = clamp((b*f - c*e) / denom, 0.0, 1.0);
    } else s = 0.0;
    // Compute point on L2 closest to S1(s) using
    // t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e t = (b*s + f) / e;
    // If t in [0,1] done. Else clamp t, recompute s for the new value
    // of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a // and clamp s to [0, 1]
    if (t < 0.0) {
      t = 0.0;
      s = clamp(-c / a, 0.0, 1.0);
    } else if (t > 1.0) {
      t = 1.0;
      s = clamp((b - c) / a, 0.0, 1.0);
    }

    //delta=(r+s*d1-t*d2)(r+s*d1-t*d2)
    double dsq=normsq3(r)+s*s*a+t*t*e+2.0*(s*dot_prod3(d1,r)+t*f+s*t*b);
    //delete d1,d2,r;
    return sqrt(dsq);
  }


  inline void test_collision23(const ssvvec& v1, const ssvvec& v2, bool* res){
    int l=0,m=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* bc=new double[3];
    double* pq=new double[3];
    double* xa=new double[3];
    double* xb=new double[3];
    double* xc=new double[3];
    double* n=new double[3]; //normal vec

    for(int k=0;k<v1.n;++k){
      double* p=v1.P+l;
      double* q=a+3;
      double* a=v2.P+m;
      double* b=a+3;
      double* c=a+6;

      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(c,b,bc);
      sub3(q,p,pq);
      sub3(a,p,xa); ///x=p
      sub3(b,p,xb);
      sub3(c,p,xc);
      cross3(ab,ac,n);

      double l1=distance22(p,q,a,b,pq,ab,xa);
      double l2=distance22(p,q,a,c,pq,ac,xa);
      double l3=distance22(p,q,b,c,pq,bc,xb);

      double d1 = dot_prod3(ab, xa);
      double d2 = dot_prod3(ac, xa);
      double d3 = dot_prod3(ab, xb);
      double d4 = dot_prod3(ac, xb);
      double d5 = dot_prod3(ab, xc);
      double d6 = dot_prod3(ac, xc);

      double va = d5*d4 - d3*d6;
      double vb = d1*d6 - d5*d2;
      double vc = d3*d2 - d1*d4;

      bool pdrin=va>0.0 && vb>0.0 && vc>0.0;

      double l4=DBL_MAX;
      double dotp=dot_prod3(xa,n);
      if(pdrin){
        l4=abs(dotp)/sqrt(normsq3(n));
      }

      sub3(a,q,xa); ///x=q
      sub3(b,q,xb);
      sub3(c,q,xc);

      double d7 = dot_prod3(ab, xa);
      double d8 = dot_prod3(ac, xa);
      double d9 = dot_prod3(ab, xb);
      double d10 = dot_prod3(ac, xb);
      double d11 = dot_prod3(ab, xc);
      double d12 = dot_prod3(ac, xc);

      double wa = d11*d10 - d9*d12;
      double wb = d7*d12 - d11*d8;
      double wc = d9*d8 - d7*d10;

      bool qdrin=wa>0.0 && wb>0.0 && wc>0.0;

      double l5=DBL_MAX;
      double dotq=dot_prod3(xa,n);
      if(qdrin){
        l5=abs(dotq)/sqrt(normsq3(n));
      }

      if(pdrin && qdrin && ((dotp>0.0 && dotq<0.0) || (dotp<0.0 && dotq>0.0))){
          l4=0.0;
      }

      double dist=l1;
      l1=min(l1,l2);
      l1=min(l1,l3);
      l1=min(l1,l4);
      l1=min(l1,l5);

      if(dist<v1.r+v2.r)res[k]=true;

      l+=6;
      m+=9;
    }
    delete ab,ac,bc,pq,xa,xb,xc,n;
  }

  inline void test_collision23(const ssvref& v1, const ssvvec& v2, bool* res){
    int m=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* bc=new double[3];
    double* pq=new double[3];
    double* xa=new double[3];
    double* xb=new double[3];
    double* xc=new double[3];
    double* n=new double[3]; //normal vec

    for(int k=0;k<v1.n;++k){
      double* p=v1.P;
      double* q=a+3;
      double* a=v2.P+m;
      double* b=a+3;
      double* c=a+6;

      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(c,b,bc);
      sub3(q,p,pq);
      sub3(a,p,xa); ///x=p
      sub3(b,p,xb);
      sub3(c,p,xc);
      cross3(ab,ac,n);

      double l1=distance22(p,q,a,b,pq,ab,xa);
      double l2=distance22(p,q,a,c,pq,ac,xa);
      double l3=distance22(p,q,b,c,pq,bc,xb);

      double d1 = dot_prod3(ab, xa);
      double d2 = dot_prod3(ac, xa);
      double d3 = dot_prod3(ab, xb);
      double d4 = dot_prod3(ac, xb);
      double d5 = dot_prod3(ab, xc);
      double d6 = dot_prod3(ac, xc);

      double va = d5*d4 - d3*d6;
      double vb = d1*d6 - d5*d2;
      double vc = d3*d2 - d1*d4;

      bool pdrin=va>0.0 && vb>0.0 && vc>0.0;

      double l4=DBL_MAX;
      double dotp=dot_prod3(xa,n);
      if(pdrin){
        l4=abs(dotp)/sqrt(normsq3(n));
      }

      sub3(a,q,xa); ///x=q
      sub3(b,q,xb);
      sub3(c,q,xc);

      double d7 = dot_prod3(ab, xa);
      double d8 = dot_prod3(ac, xa);
      double d9 = dot_prod3(ab, xb);
      double d10 = dot_prod3(ac, xb);
      double d11 = dot_prod3(ab, xc);
      double d12 = dot_prod3(ac, xc);

      double wa = d11*d10 - d9*d12;
      double wb = d7*d12 - d11*d8;
      double wc = d9*d8 - d7*d10;

      bool qdrin=wa>0.0 && wb>0.0 && wc>0.0;

      double l5=DBL_MAX;
      double dotq=dot_prod3(xa,n);
      if(qdrin){
        l5=abs(dotq)/sqrt(normsq3(n));
      }

      if(pdrin && qdrin && ((dotp>0.0 && dotq<0.0) || (dotp<0.0 && dotq>0.0))){
          l4=0.0;
      }

      double dist=l1;
      l1=min(l1,l2);
      l1=min(l1,l3);
      l1=min(l1,l4);
      l1=min(l1,l5);

      if(dist<v1.r+v2.r)res[k]=true;

      m+=9;
    }
    delete ab,ac,bc,pq,xa,xb,xc,n;
  }

  inline void test_collision23(const ssvvec& v1, const ssvref& v2, bool* res){
    int l=0;
    double* ab=new double[3];
    double* ac=new double[3];
    double* bc=new double[3];
    double* pq=new double[3];
    double* xa=new double[3];
    double* xb=new double[3];
    double* xc=new double[3];
    double* n=new double[3]; //normal vec

    for(int k=0;k<v1.n;++k){
      double* p=v1.P+l;
      double* q=a+3;
      double* a=v2.P;
      double* b=a+3;
      double* c=a+6;

      sub3(b,a,ab);
      sub3(c,a,ac);
      sub3(c,b,bc);
      sub3(q,p,pq);
      sub3(a,p,xa); ///x=p
      sub3(b,p,xb);
      sub3(c,p,xc);
      cross3(ab,ac,n);

      double l1=distance22(p,q,a,b,pq,ab,xa);
      double l2=distance22(p,q,a,c,pq,ac,xa);
      double l3=distance22(p,q,b,c,pq,bc,xb);

      double d1 = dot_prod3(ab, xa);
      double d2 = dot_prod3(ac, xa);
      double d3 = dot_prod3(ab, xb);
      double d4 = dot_prod3(ac, xb);
      double d5 = dot_prod3(ab, xc);
      double d6 = dot_prod3(ac, xc);

      double va = d5*d4 - d3*d6;
      double vb = d1*d6 - d5*d2;
      double vc = d3*d2 - d1*d4;

      bool pdrin=va>0.0 && vb>0.0 && vc>0.0;

      double l4=DBL_MAX;
      double dotp=dot_prod3(xa,n);
      if(pdrin){
        l4=abs(dotp)/sqrt(normsq3(n));
      }

      sub3(a,q,xa); ///x=q
      sub3(b,q,xb);
      sub3(c,q,xc);

      double d7 = dot_prod3(ab, xa);
      double d8 = dot_prod3(ac, xa);
      double d9 = dot_prod3(ab, xb);
      double d10 = dot_prod3(ac, xb);
      double d11 = dot_prod3(ab, xc);
      double d12 = dot_prod3(ac, xc);

      double wa = d11*d10 - d9*d12;
      double wb = d7*d12 - d11*d8;
      double wc = d9*d8 - d7*d10;

      bool qdrin=wa>0.0 && wb>0.0 && wc>0.0;

      double l5=DBL_MAX;
      double dotq=dot_prod3(xa,n);
      if(qdrin){
        l5=abs(dotq)/sqrt(normsq3(n));
      }

      if(pdrin && qdrin && ((dotp>0.0 && dotq<0.0) || (dotp<0.0 && dotq>0.0))){
          l4=0.0;
      }

      double dist=l1;
      if(l2<l1)l1=l2;
      if(l3<l1)l1=l3;
      if(l4<l1)l1=l4;
      if(l5<l1)l1=l5;

      if(dist<v1.r+v2.r)res[k]=true;

      l+=6;
    }
    delete ab,ac,bc,pq,xa,xb,xc,n;
  }

  inline void test_collision33(const ssvvec& v1, const ssvvec& v2, bool* res){
    int l=0;
    double* ab1=new double[3];
    double* ac1=new double[3];
    double* bc1=new double[3];
    double* ab2=new double[3];
    double* ac2=new double[3];
    double* bc2=new double[3];
    double* a1a2=new double[3];
    double* a1b2=new double[3];
    double* a1c2=new double[3];
    double* b1a2=new double[3];
    double* b1b2=new double[3];
    double* b1c2=new double[3];
    double* c1a2=new double[3];
    double* c1b2=new double[3];
    double* c1c2=new double[3];

    for(int k=0;k<v1.n;++k){
      double* a1=v1.P+l;
      double* b1=a1+3;
      double* c1=a1+6;
      double* a2=v2.P+l;
      double* b2=a2+3;
      double* c2=a2+6;

      sub3(b1,a1,ab1);
      sub3(c1,a1,ac1);
      sub3(c1,b1,bc1);
      sub3(b2,a2,ab2);
      sub3(c2,a2,ac2);
      sub3(c2,b2,bc2);
      sub3(a2,a1,a1a2);
      sub3(b2,a1,a1b2);
      sub3(c2,a1,a1c2);
      sub3(a2,b1,b1a2);
      sub3(b2,b1,b1b2);
      sub3(c2,b1,b1c2);
      sub3(a2,c1,c1a2);
      sub3(b2,c1,c1b2);
      sub3(c2,c1,c1c2);

      double l[15];
      l[0]=distance22(a2,b2,a1,b1,ab2,ab1,a1a2);
      l[1]=distance22(a2,c2,a1,b1,ac2,ab1,a1a2);
      l[2]=distance22(b2,c2,a1,b1,bc2,ab1,a1b2);
      l[3]=distance22(a2,b2,a1,c1,ab2,ac1,a1a2);
      l[4]=distance22(a2,c2,a1,c1,ac2,ac1,a1a2);
      l[5]=distance22(b2,c2,a1,c1,bc2,ac1,a1b2);
      l[6]=distance22(a2,b2,b1,c1,ab2,bc1,b1a2);
      l[7]=distance22(a2,c2,b1,c1,ac2,bc1,b1a2);
      l[8]=distance22(b2,c2,b1,c1,bc2,bc1,b1b2);

      l[9]=distance13(a2,a1,b1,c1,ab1,ac1,a1a2,b1a2,c1a2);
      l[10]=distance13(b2,a1,b1,c1,ab1,ac1,a1b2,b1b2,c1b2);
      l[11]=distance13(c2,a1,b1,c1,ab1,ac1,a1c2,b1c2,c1c2);

      l[12]=distance13mod(a1,a2,b2,c2,ab2,ac2,a1a2,a1b2,a1c2);
      l[13]=distance13mod(b1,a2,b2,c2,ab2,ac2,b1a2,b1b2,b1c2);
      l[14]=distance13mod(c1,a2,b2,c2,ab2,ac2,c1a2,c1b2,c1c2);

      double rsum=v1.r+v2.r;
      for(int i=0;i<15;++i)if(l[i]<rsum){res[k]==true;break;}

      l+=9;
    }
    delete ab1,ac1,bc1,ab2,ac2,bc2,a1a2,a1b2,a1c2,b1a2,b1b2,b1c2,c1a2,c1b2,c1c2;
  }

  inline void test_collision33(const ssvref& v1, const ssvvec& v2, bool* res){
    int l=0;
    double* ab1=new double[3];
    double* ac1=new double[3];
    double* bc1=new double[3];
    double* ab2=new double[3];
    double* ac2=new double[3];
    double* bc2=new double[3];
    double* a1a2=new double[3];
    double* a1b2=new double[3];
    double* a1c2=new double[3];
    double* b1a2=new double[3];
    double* b1b2=new double[3];
    double* b1c2=new double[3];
    double* c1a2=new double[3];
    double* c1b2=new double[3];
    double* c1c2=new double[3];

    for(int k=0;k<v1.n;++k){
      double* a1=v1.P;
      double* b1=a1+3;
      double* c1=a1+6;
      double* a2=v2.P+l;
      double* b2=a2+3;
      double* c2=a2+6;

      sub3(b1,a1,ab1);
      sub3(c1,a1,ac1);
      sub3(c1,b1,bc1);
      sub3(b2,a2,ab2);
      sub3(c2,a2,ac2);
      sub3(c2,b2,bc2);
      sub3(a2,a1,a1a2);
      sub3(b2,a1,a1b2);
      sub3(c2,a1,a1c2);
      sub3(a2,b1,b1a2);
      sub3(b2,b1,b1b2);
      sub3(c2,b1,b1c2);
      sub3(a2,c1,c1a2);
      sub3(b2,c1,c1b2);
      sub3(c2,c1,c1c2);

      double l[15];
      l[0]=distance22(a2,b2,a1,b1,ab2,ab1,a1a2);
      l[1]=distance22(a2,c2,a1,b1,ac2,ab1,a1a2);
      l[2]=distance22(b2,c2,a1,b1,bc2,ab1,a1b2);
      l[3]=distance22(a2,b2,a1,c1,ab2,ac1,a1a2);
      l[4]=distance22(a2,c2,a1,c1,ac2,ac1,a1a2);
      l[5]=distance22(b2,c2,a1,c1,bc2,ac1,a1b2);
      l[6]=distance22(a2,b2,b1,c1,ab2,bc1,b1a2);
      l[7]=distance22(a2,c2,b1,c1,ac2,bc1,b1a2);
      l[8]=distance22(b2,c2,b1,c1,bc2,bc1,b1b2);

      l[9]=distance13(a2,a1,b1,c1,ab1,ac1,a1a2,b1a2,c1a2);
      l[10]=distance13(b2,a1,b1,c1,ab1,ac1,a1b2,b1b2,c1b2);
      l[11]=distance13(c2,a1,b1,c1,ab1,ac1,a1c2,b1c2,c1c2);

      l[12]=distance13mod(a1,a2,b2,c2,ab2,ac2,a1a2,a1b2,a1c2);
      l[13]=distance13mod(b1,a2,b2,c2,ab2,ac2,b1a2,b1b2,b1c2);
      l[14]=distance13mod(c1,a2,b2,c2,ab2,ac2,c1a2,c1b2,c1c2);

      double rsum=v1.r+v2.r;
      for(int i=0;i<15;++i)if(l[i]<rsum){res[k]==true;break;}

      l+=9;
    }
    delete ab1,ac1,bc1,ab2,ac2,bc2,a1a2,a1b2,a1c2,b1a2,b1b2,b1c2,c1a2,c1b2,c1c2;
  }

  inline void test_collision33(const ssvvec& v1, const ssvref& v2, bool* res){
    test_collision33(v2,v1,res);
  }


}


#endif // SSV_ELEMENT_H
