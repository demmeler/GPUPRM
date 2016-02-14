#include <iostream>

#include "lib/geo4.h"

#include "lib/robot.h"
#include "lib/collision4.h"

using namespace std;
using namespace geo4;
using namespace collision4;

const float pi=3.14159265358;

int main()
{ 

#if 0
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
#endif

#if 0
  Robot<3> robot;
  robot.a[0]=0.0;
  robot.alpha[0]=0.0;
  robot.q[0]=0.0;
  robot.d[0]=0.0;
  robot.types[0]=rotational;

  robot.a[1]=1.0;
  robot.alpha[1]=0.0;
  robot.q[1]=0.0;
  robot.d[1]=0.0;
  robot.types[1]=rotational;

  robot.a[2]=0;
  robot.alpha[2]=0;
  robot.q[2]=0;
  robot.d[2]=0;
  robot.types[2]=prismatic;

  float q[3]={0.0,pi/4.0,0.9};

  Kinematics<3> kin;
  kin.robot=&robot;

  kin.calculate(&q[0]);

  t4print(kin.trafos[0]);
  t4print(kin.trafos[1]);
  t4print(kin.trafos[2]);

  float4 p=make_float4(1.0,0.0,0.0);

  f4print(p);

  kin.trafos[2].apply(p);

  f4print(p);
#endif



  cout<<"test.cpp"<<endl;
  return 0;
}

