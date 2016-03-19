#ifndef ROBOTCONFIGSPACE_H
#define ROBOTCONFIGSPACE_H

#include "configspace.hpp"
#include "robot.h"
#include "polytope.h"

template<int ndof> class Kinematics;

namespace collision4{
  class polytope4;
  class polytope4data;
}

template<int ndof>
class RobotConfigspace : public Configspace<ndof>
{
public:
  //! robot_: Robot Object with Denavit-Hartenberg data
  //! polydata_: Object containing compressed Polytope data of robot arm and surrounding
  //! mins[ndof]: maximal q values
  //! maxs[ndof]: minimal   -"-
  //! dq:         edge resolution (stepsize)
  RobotConfigspace(const Robot<ndof>* robot_,
                   const polytope *polys_,
                   const int* sys_,
                   const int N_,
                   const float* mins_, const float* maxs_, const float dq_,
                   const int nbuf_, const int numthreadsmax_);
  RobotConfigspace(const Robot<ndof>* robot_,
                   const polytope *polys_,
                   const int* sys_,
                   const int N_,
                   const int *from_, const int *to_,
                   const int M_,
                   const float* mins_, const float* maxs_, const float dq_,
                   const int nbuf_, const int numthreadsmax_);
private:
  void construct(const Robot<ndof>* robot_,
                 const polytope *polys_,
                 const int* sys_,
                 const int N_,
                 const float* mins_, const float* maxs_, const float dq_,
                 const int nbuf_, const int numthreadsmax_);
public:

  //!initialization function copy polydata to gpu etc..
  int init();
  int clear();

  //!
  //! indicator functions on CPU
  //!

  //! indicator function of obstacles
  //! q: length d*N, array of structures: q[N*k+i]= k-th component of i-th q-vector
  //! res: length N
  int indicator(const float* q, int* res, const int N, const int offset);
  //!case N=1
  int indicator(const float* q);

  //!
  //! indicator functions on GPU
  //!

  //! checks if indicator function = 1 somewhere on the line between (qs[i],qs[i+N]) and (qe[i],qe[i+N]) for all i=1,...,N-1
  //! res is return value
  //! number of pairs
  int indicator2(const float* qs, const float* qe, int *res, const int N, const int offset);
  //! same paircheck as above, but with compressed storage:
  //! checks pairs: (qs[i],...) ->  (qe(posqe[i]),...) , ...., (qe[posqe[i]+numqe[i]-1],...) for i=0,...,M-1
  int indicator2(const float* qs, const int M, const float* qe, int *res, const int *posqe, const int *numqe, const int offset);


  //! boundary check on CPU

  //! structure like indicator function
  //! returns if lies in boundaries
  int check_boundaries(const float* q, int* res, const int N, const int offset);
  //! for N=1
  int check_boundaries(const float* q);

  //! get functions

  //! minimal value of qi, i=0...d-1
  float min(int i){return mins[i];}
  //! maximal value of qi, i=0...d-1
  float max(int i){return maxs[i];}
  //! dq used for indicator 2
  float deltaq(){return dq;}

  //! dimension d
  int dim(){return ndof;}

private:

  Kinematics<ndof>* kin;

  //! robot data
  const Robot<ndof>* robot; //host
  Robot<ndof>* robotdev;    //GPU

  //!polytope data
  struct{
    const polytope* polys;
    const int* sys;
    int N;
    const int *from, *to;
    int M;
  }polylist;

  collision4::polytope4data* polydata;            //data on host
  collision4::polytope4data* polydatadev_hostref; //help obj on host
  collision4::polytope4data* polydatadev;         //pointer to GPU obj

  //!flag if GPU initialized
  bool devloaded;

  //!device storage for indicator querries
  float* qdevbufferfrom; //GPU length nbufqfrom
  int nbufqfrom;
  float* qdevbufferto;   //GPU length nbufqto
  int nbufqto;
  int* resdevbuffer;    //GPU length nbufres
  int nbufres;

  int* testnumdev; //GPU length nbuftest
  int* testposdev; //GPU length nbuftest
  int nbuftest;

  int* resdevbufferext; //GPU length numthreadsmax
  int* resbufferext;    //host length numthreadsmax -> for implementation without cuda
  int numthreadsmax;


  float mins[ndof];
  float maxs[ndof];
  float dq;

};

#endif // ROBOTCONFIGSPACE_H
