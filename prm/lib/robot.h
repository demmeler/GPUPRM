#ifndef ROBOT_H
#define ROBOT_H


//enum jointtype{rotational=0, prismatic=1};
const int rotational=0;
const int prismatic=1;


template<int ndof>
struct Robot{
  float a[ndof];
  float alpha[ndof];
  float q[ndof];
  float d[ndof];

  //jointtype types[ndof];
  int types[ndof];

  //Robot(int ndof_, float* a_, float* alpha_, float* q_, float* d_):
  //  ndof(ndof_),a(a_),alpha(alpha_),q(q_),d(d_){}
};


#endif // ROBOT_H
