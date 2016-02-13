#ifndef VEC4_H
#define VEC4_H

#include <cuda.h>

namespace vec4{

  struct trafo{
    struct{ float4 r1,r2,r3,r4; } row;
    float4 data[4];
  };


}






#endif // VEC4_H
