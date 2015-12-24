#ifndef SSV_ELEMENT_HPP
#define SSV_ELEMENT_HPP

#include <math.h>

namespace ssv{

  class ssvele;
  class ssv1;
  class ssv2;
  class ssv3;

  class ssvele
  {
  public:
    ssvele(){}
  };

  class ssv1: ssvele{
  public:
    ssv1();

    inline double dist(ssv1& ele);
    inline double dist(ssv2& ele);
    inline double dist(ssv3& ele);
  private:
    double p[3];
    double r2;
  };

  class ssv2: ssvele{
  public:
    ssv2();

    inline double dist(ssv1& ele);
    inline double dist(ssv2& ele);
    inline double dist(ssv3& ele);

  private:
    double p[3];
    double d[3];
    double l;
    double r2;
  };

  class ssv3: ssvele{
    ssv3();

    inline double dist(ssv1& ele);
    inline double dist(ssv2& ele);
    inline double dist(ssv3& ele);

  private:
    double p1[3];
    double d1[3];
    double d2[3];
    double l1, l2;
    double r2;
  };



  inline double ssv1::dist(ssv1& ele){
    return p[0]*ele.p[0]+p[1]*ele.p[1]+p[2]*ele.p[2];
  }

  inline double ssv1::dist(ssv2& ele){

  }

  inline double ssv1::dist(ssv3& ele){

  }

  inline double ssv2::dist(ssv1& ele){
    return ele.dist(*this);
  }

  inline double ssv2::dist(ssv2& ele){

  }

  inline double ssv2::dist(ssv3& ele){

  }

  inline double ssv3::dist(ssv1& ele){
    return ele.dist(*this);
  }

  inline double ssv3::dist(ssv2& ele){
    return ele.dist(*this);
  }

  inline double ssv3::dist(ssv3& ele){

  }















}


#endif // SSV_ELEMENT_H
