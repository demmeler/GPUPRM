#include <iostream>

//#include "lib/collision.hpp"
//#include "lib/util.hpp"
#include "lib/geo4.h"

using namespace std;
using namespace geo4;
//using namespace collision;

const float pi=3.14159265358;

int main()
{ 
  float4 u=make_float4(sqrt(2),0.0,0.0);
  trafo4 T(0.0,0.0,pi/4,1.0);

 f4print(u);
 t4print(T);

  T.apply(u);

  f4print(u);


  float4 v;

  f4print(v);

  T.apply(u,v);

  f4print(u);
  f4print(v);

  cout<<"test.cpp"<<endl;
  return 0;
}

