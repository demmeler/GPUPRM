#ifndef POLYTOPE_H
#define POLYTOPE_H


  struct polytope{
    struct vec4{
      double x,y,z,w;
    };
    vec4* vertices; //array of structs
    int n;
    //!edges saved in crs format
    int* dsp;
    int* cnt;
    int* dest;
    int m;
  };



#endif // POLYTOPE_H
