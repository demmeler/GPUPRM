#ifndef MATVEC_HPP
#define MATVEC_HPP

#include <math.h>

namespace matvec{

#define clamp(x,min,max) (x<min?min:(x>max?max:x))
#define sign(x) (x>0.0 ? 1.0 : -1.0)
#define minimum(x,y) (x<y?x:y)

  inline void normalize3(double* x){
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

  inline void add3(double* x, double* y, double* res){
    res[0]=x[0]+y[0];
    res[1]=x[1]+y[1];
    res[2]=x[2]+y[2];
  }

  inline void add2(double* x, double* y, double* res){
    res[0]=x[0]+y[0];
    res[1]=x[1]+y[1];
  }

  inline void sub3(double* x, double* y, double* res){
    res[0]=x[0]-y[0];
    res[1]=x[1]-y[1];
    res[2]=x[2]-y[2];
  }

  inline void mult3(double c, double* x, double* res){
    res[0]=c*x[0];
    res[1]=c*x[1];
    res[2]=c*x[2];
  }

  inline void lincomb3(double a, double* x, double b, double* y, double* res){
    res[0]=a*x[0]+b*y[0];
    res[1]=a*x[1]+b*y[1];
    res[2]=a*x[2]+b*y[2];
  }

  inline double dot_prod3(double* x, double* y){
    return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
  }

  inline void copy3(double* from, double* to){
    to[0]=from[0];
    to[1]=from[1];
    to[2]=from[2];
  }

  inline void copy2(double* from, double* to){
    to[0]=from[0];
    to[1]=from[1];
  }

  inline void scale3(double* x, double f){
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

  inline void cross3(double* x, double* y, double* z){
    z[0]=x[1]*y[2]-x[2]*y[1];
    z[1]=x[2]*y[0]-x[0]*y[2];
    z[2]=x[0]*y[1]-x[1]*y[0];
  }


  inline double dist3(double* x, double* y){
    double d0=x[0]-y[0];
    double d1=x[1]-y[1];
    double d2=x[2]-y[2];
    return sqrt(d0*d0+d1*d1+d2*d2);
  }

  inline double normsq3(double* x){
    //TODO: optimieren?
    return dot_prod3(x,x);
  }

  inline double det3(double* x, double* y, double* z){
    double a=x[0]*(y[1]*z[2]-y[2]*z[1]);
    double b=x[1]*(y[2]*z[0]-y[0]*z[2]);
    double c=x[2]*(y[0]*z[1]-y[1]*z[0]);
    return a+b+c;
  }

}

#endif // MATVEC_HPP
