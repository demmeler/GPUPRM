#ifndef ARRAYCONFIGSPACE_HPP
#define ARRAYCONFIGSPACE_HPP

#include <prm/configspace.hpp>

class ArrayConfigspace : public Configspace<2>
{
public:
  ArrayConfigspace(const int* array_, int b_, int h_, float minx, float maxx, float miny, float maxy);
  ~ArrayConfigspace();

  //!initialization function
  int init(const int ressource_rank=0, const int ressource_size=1);
  //! indicator function of obstacles
  //! q: length d*N, array of structures: q[N*k+i]= k-th component of i-th q-vector
  //! res: length N
  int indicator(const float* q, int* res, const int N, const int offset);
  //!case N=1
  int indicator(const float* q);
  //! checks if indicator function = 1 somewhere on the line between (qs[i],qs[i+N]) and (qe[i],qe[i+N])
  //! res is return value
  int indicator2(const float* qs, const float* qe, int *res, const int N, const int offset);
  int indicator2_async(const float* qs, const float* qe, int *res, const int N, const int offset, int &request);
  int indicator2_async_wait(int request);

  //! same paircheck as above, but with compressed storage:
  //! checks pairs: (qs[i],...) ->  (qe(posqe[i]),...) , ...., (qe[posqe[i]+numqe[i]-1],...) for i=0,...,M-1
  int indicator2(const float* qs, const int M, const float* qe, int *res, const int *posqe, const int *numqe, const int offset);
  //! structure like indicator function
  //! returns if lies in boundaries
  int check_boundaries(const float* q, int* res, int N, const int offset);
  //! for N=1
  int check_boundaries(const float* q);

  //! minimal value of qi, i=0...d-1
  float min(int i);
  //! maximal value of qi, i=0...d-1
  float max(int i);
  //! dq used for indcator2
  float deltaq(){return dq;}
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
  float dq;
};

#endif // ARRAYCONFIGSPACE_HPP
