#ifndef KINEMATICS_H
#define KINEMATICS_H


#include "cuda_head.h"

#include "geo4.h"
//using namespace geo4;

#include "robot.h"

template<int ndof>
class Kinematics{
public:
  const Robot<ndof>* robot;
  geo4::trafo4 trafos[ndof+1];

  qualifier Kinematics(const Robot<ndof>* robot_): robot(robot_){}

  qualifier void calculate(const float* q, int offset);
};



///   **************************
///   *       Kinematics       *
///   *    implementations     *
///   **************************
template<int ndof>
qualifier void Kinematics<ndof>::calculate(const float* q, int offset){
  trafos[0].set(0.0, 0.0, 0.0, 0.0);
  for(int i=0;i<ndof;++i){
      float qset=robot->q[i], dset=robot->d[i];
      if(robot->types[i]==rotational) qset=q[i];
      else dset=q[i*offset];
      trafos[i+1].set(robot->a[i],robot->alpha[i],qset,dset);
      trafos[i+1].lapply(trafos[i]);
  }
}

#endif // KINEMATICS_H
