#ifndef ARRAYCONFIGSPACE_HPP
#define ARRAYCONFIGSPACE_HPP

#include "configspace.hpp"

class ArrayConfigspace : public Configspace
{
public:
  ArrayConfigspace(const int* array, float minx, float maxx, float miny, float maxy);
  ~ArrayConfigspace();

  //!initialization function
  int init();
  //! indicator function of obstacles
  //! q: length d*N, array of structures: q[N*k+i]= k-th component of i-th q-vector
  //! res: length N
  int indicator(const float* q, int* res, int N);
  //! structure like indicator function
  //! returns if lies in boundaries
  virtual int check_boundaries(const float* q, int* res, int N);
  //! for N=1
  virtual int check_boundaries(const float* q, int* res)=0;

  //! minimal value of qi, i=0...d-1
  float min(int i);
  //! maximal value of qi, i=0...d-1
  float max(int i);
  //! dimension d
  int dim(){return 2;}
private:
  float mins[2];
  float maxs[2];
  float factor[2];
  int* array;
  int b;
  int h;
  int n;
};

#endif // ARRAYCONFIGSPACE_HPP
