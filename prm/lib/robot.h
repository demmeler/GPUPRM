#ifndef ROBOT_H
#define ROBOT_H

#undef qualifier
#ifdef CUDA_IMPLEMENTATION
  #include <cuda.h>
  #define qualifier __device__
  #define cudaonly(x) x
  #define devonly(x)
#else
  #define qualifier inline
  #define cudaonly(x)
  #define devonly(x) x
#endif


#include "geo4.h"
using namespace geo4;

enum jointtype{rotational, prismatic};

template<int ndof>
struct Robot{
  float a[ndof];
  float alpha[ndof];
  float q[ndof];
  float d[ndof];

  jointtype types[ndof];

  //Robot(int ndof_, float* a_, float* alpha_, float* q_, float* d_):
  //  ndof(ndof_),a(a_),alpha(alpha_),q(q_),d(d_){}
};

template<int ndof>
class Kinematics{
public:
  const Robot<ndof>* robot;
  trafo4 trafos[ndof];

  qualifier void calculate(float* q);
};



///   **************************
///   *       Kinematics       *
///   *    implementations     *
///   **************************

template<int ndof>
qualifier void Kinematics<ndof>::calculate(float* q){
  float qset=robot->q[0], dset=robot->d[0];
  if(robot->types[0]==rotational) qset=q[0];
  else dset=q[0];
  trafos[0].set(robot->a[0],robot->alpha[0],qset,dset);
  for(int i=1;i<ndof;++i){
      float qset=robot->q[i], dset=robot->d[i];
      if(robot->types[i]==rotational) qset=q[i];
      else dset=q[i];
      trafos[i].set(robot->a[i],robot->alpha[i],qset,dset);
      trafos[i].lapply(trafos[i-1]);
  }
}


#endif // ROBOT_H
