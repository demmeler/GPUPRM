#ifndef SSV_ELEMENT_HPP
#define SSV_ELEMENT_HPP

#include <math.h>
#include <trafo.hpp>

namespace ssv{

  class ssvele;
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
    double r2;
  };

  inline double dist(ssvele ele1, ssvele ele2);

}


#endif // SSV_ELEMENT_H
