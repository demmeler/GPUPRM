#ifndef POLYTOPE4_H
#define POLYTOPE4_H

#include "geo4.h"


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

}

#endif // POLYTOPE4_H
