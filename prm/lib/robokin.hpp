#ifndef ROBOKIN_HPP
#define ROBOKIN_HPP


#include <math.h>
#include <vector>
#include "ssv.hpp"
#include "util.hpp"

using namespace std;
using namespace ssv;

typedef vector<double> vecd;
typedef vector<int> veci;

class robogeo;
class robokin;

class trafomat{
public:
  //!Trafo nach Denavit Hartenberg
  trafomat(double q, double alpha, double a, double d){
    double cq=cos(q);
    double sq=sin(q);
    double ca=cos(alpha);
    double sa=sin(alpha);

    val[0]=cq;    val[1]=-sq;   val[2]=0.0;  val[3]=a;
    val[4]=sq*ca; val[5]=cq*ca; val[6]=-sa;  val[7]=-sa*d;
    val[8]=sq*sa; val[9]=cq*sa; val[10]=ca;  val[11]=ca*d;
  }
  trafomat(){}

  inline void apply(double* vec, double* res){
    double vec0=vec[0];
    double vec1=vec[1];
    double vec2=vec[2];

    res[0]=val[0]*vec0+val[1]*vec1+val[2]*vec2+val[3];
    res[1]=val[4]*vec0+val[5]*vec1+val[6]*vec2+val[7];
    res[2]=val[8]*vec0+val[9]*vec1+val[10]*vec2+val[11];
  }

  inline void apply(const trafomat &trafo, trafomat &restrafo){
    double* tval=&(trafo.val);
    double* res=&(restrafo->val);

    for(int i=0;i<3;++i){
      int i0,i4,i8;
      double vec0=tval[i0=0+i];
      double vec1=tval[i4=4+i];
      double vec2=tval[i8=8+i];
      res[i0]=val[0]*vec0+val[1]*vec1+val[2]*vec2;
      res[i4]=val[4]*vec0+val[5]*vec1+val[6]*vec2;
      res[i8]=val[8]*vec0+val[9]*vec1+val[10]*vec2;
    }
    double vec0=tval[3];
    double vec1=tval[7];
    double vec2=tval[11];
    res[3]=val[0]*vec0+val[1]*vec1+val[2]*vec2+val[3];
    res[7]=val[4]*vec0+val[5]*vec1+val[6]*vec2+val[7];
    res[11]=val[8]*vec0+val[9]*vec1+val[10]*vec2+val[11];
    return
  }

private:
  double val[12];
};



//!storage class for denavit-hartenberg parameters
class robogeo{
  friend class robokin;
public:
  robogeo(const vecd &alpha_, const vecd &a_, const vecd &d_, const vector<ssvele> &ssvs_, const veci &disp_, const veci &num_,
          vector<ssvele> obsts_, vector<ssvele*> pairs_):
    dim(alpha_.size()),nssv(ssvs_.size()),alpha(alpha_), a(a_),d(d_),ssvs(ssvs_),disp(disp_),num(num_), obsts(obsts_), pairs(pairs_)
  {
    check(alpha_.size()==dim);
    check(a_.size()==dim);
    check(d_.size()==dim);
    check(disp_.size()==dim+1);
    check(num_.size()==dim+1);
    int numsum=0;
    for(int i=0;i<dim+1;++i)numsum+=num[i];
    check(numsum==nssv);
    check(pairs.size()%2==0);

  }
  ~robogeo();

protected:
  int dim;
  vecd alpha;
  vecd a;
  vecd d;

  //!robot parts ssv elements, ordered by dof:
  //i-th dof: ssvs[disp[i]],...,ssvs[disp[i]+num[i]-1]
  //num[0]+...+num[dim]=nssz
  int nssv;
  vector<ssvele> ssvs;
  veci disp;
  veci num;
  //!obstacle ssvs in global coordinates
  vector<ssvele> obsts;
  //!ssv pairs (pairs[0],pairs[1]),(pairs[2],pairs[3]),... which should be tested
  vector<ssvele*> pairs;
};

class robokin{

public:
  //!calculate
  static void calc_kinematics(double* q, const robogeo &geo){
  }

  /*
  robokin(robogeo* geo_){geo=geo_;}
  ~robokin();

private:
  robogeo* geo;
  */
};


#endif
