#ifndef ROBOT_H
#define ROBOT_H

#include "cuda_head.h"


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
  trafo4 trafos[ndof+1];

  qualifier Kinematics(const Robot<ndof>* robot_): robot(robot_){}

  qualifier void calculate(float* q, int offset);
};



///   **************************
///   *       Kinematics       *
///   *    implementations     *
///   **************************
template<int ndof>
qualifier void Kinematics<ndof>::calculate(float* q, int offset){
  trafos[0].set(0.0, 0.0, 0.0, 0.0);
  for(int i=0;i<ndof;++i){
      float qset=robot->q[i], dset=robot->d[i];
      if(robot->types[i]==rotational) qset=q[i];
      else dset=q[i*offset];
      trafos[i+1].set(robot->a[i],robot->alpha[i],qset,dset);
      trafos[i+1].lapply(trafos[i]);
  }
}


#endif // ROBOT_H
