#ifndef POLYTOPE_H
#define POLYTOPE_H

namespace collision{

  struct Polytope{
    double* vertices;
    int n;
    //!edges saved in crs format
    int* dsp;
    int* cnt;
    int* dest;
    int m;
  };

}


#endif // POLYTOPE_H
