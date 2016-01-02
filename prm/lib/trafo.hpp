#ifndef TRAFO_H
#define TRAFO_H

#include <math.h>

class trafo{
public:
  //!Trafo nach Denavit Hartenberg
  trafo(double q, double alpha, double a, double d){
    double cq=cos(q);
    double sq=sin(q);
    double ca=cos(alpha);
    double sa=sin(alpha);

    val[0]=cq;    val[3]=-sq;   val[6]=0.0;  val[9]=a;
    val[1]=sq*ca; val[4]=cq*ca; val[7]=-sa;  val[10]=-sa*d;
    val[2]=sq*sa; val[5]=cq*sa; val[8]=ca;   val[11]=ca*d;
  }
  trafo(){}

  //!apply trafo: res=R*vec+t
  inline void apply(const double* vec, double* res) const{
    double vec0=vec[0];
    double vec1=vec[1];
    double vec2=vec[2];

    res[0]=val[0]*vec0+val[3]*vec1+val[6]*vec2+val[9];
    res[1]=val[1]*vec0+val[4]*vec1+val[7]*vec2+val[10];
    res[2]=val[2]*vec0+val[5]*vec1+val[8]*vec2+val[11];
  }

  //!apply tafo to n vectors
  inline void apply(const double* vec, double* res, const int n) const{
    for(i=0;i<3*n;i+=3){
      int i1=i+1;
      int i2=i+2;
      double vec0=vec[i];
      double vec1=vec[i1];
      double vec2=vec[i2];
      res[i]=val[0]*vec0+val[3]*vec1+val[6]*vec2+val[9];
      res[i1]=val[1]*vec0+val[4]*vec1+val[7]*vec2+val[10];
      res[i2]=val[2]*vec0+val[5]*vec1+val[8]*vec2+val[11];
    }
  }

  //!apply inverse trafo: res=R^T*(vec-t)
  inline void apply_inv(const double* vec, double* res)const{
    double vec0=vec[0]-val[9];
    double vec1=vec[1]-val[10];
    double vec2=vec[2]-val[11];

    res[0]=val[0]*vec0+val[1]*vec1+val[2]*vec2;
    res[1]=val[3]*vec0+val[4]*vec1+val[5]*vec2;
    res[2]=val[6]*vec0+val[7]*vec1+val[8]*vec2;
  }

  inline void apply_rot(const double* vec, double* res)const{
    double vec0=vec[0];
    double vec1=vec[1];
    double vec2=vec[2];

    res[0]=val[0]*vec0+val[3]*vec1+val[6]*vec2;
    res[1]=val[1]*vec0+val[4]*vec1+val[7]*vec2;
    res[2]=val[2]*vec0+val[5]*vec1+val[8]*vec2;
  }

  inline void apply_rot_inv(const double* vec, double* res)const{
    double vec0=vec[0];
    double vec1=vec[1];
    double vec2=vec[2];

    res[0]=val[0]*vec0+val[1]*vec1+val[2]*vec2;
    res[1]=val[3]*vec0+val[4]*vec1+val[5]*vec2;
    res[2]=val[6]*vec0+val[7]*vec1+val[8]*vec2;
  }


  inline void apply(const trafo &t, trafo &tres) const{
    const double* tval=&(t.val[0]);
    double* res=&(tres.val[0]);

    for(int i=0;i<9;i+=3){
      int i0,i4,i8;
      double vec0=tval[i0=i];
      double vec1=tval[i4=i+1];
      double vec2=tval[i8=i+2];
      res[i0]=val[0]*vec0+val[3]*vec1+val[6]*vec2;
      res[i4]=val[1]*vec0+val[4]*vec1+val[7]*vec2;
      res[i8]=val[2]*vec0+val[5]*vec1+val[8]*vec2;
    }
    double vec0=tval[9];
    double vec1=tval[10];
    double vec2=tval[11];
    res[9]= val[0]*vec0+val[3]*vec1+val[6]*vec2+val[9];
    res[10]=val[1]*vec0+val[4]*vec1+val[7]*vec2+val[10];
    res[11]=val[2]*vec0+val[5]*vec1+val[8]*vec2+val[11];
  }

  inline void apply(double q, double alpha, double a, double d){
    //TODO: optimize!!
    trafo tres;
    trafo t(q,alpha,a,d);
    apply(t,tres);
    set(tres);
  }

  inline const double* get_translation() const{
    return val+9;
  }

  inline trafo& operator=(const trafo& t){
    //unroll...
    for(int i=0;i<12;++i)val[i]=t.val[i];
  }

  inline void set(const trafo& t){
    //unroll...
    for(int i=0;i<12;++i)val[i]=t.val[i];
  }

private:
  double val[12];
};

#endif // TRAFO_H
