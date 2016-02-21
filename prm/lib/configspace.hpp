#ifndef CONFIGSPACE_HPP
#define CONFIGSPACE_HPP

class Configspace
{
public:

  //!initialization function
  virtual int init()=0;
  //! indicator function of obstacles
  //! q: length d*N, array of structures: q[N*k+i]= k-th component of i-th q-vector
  //! res: length N
  virtual int indicator(const float* q, int* res, int N)=0;
  //! checks if indicator function = 1 somewhere on the line between (qs[i],qs[i+N]) and (qe[i],qe[i+N])
  //! res is return value
  //! number of pairs
  virtual int indicator2(const float* qs, const float* qe, float *res, const int N)=0;
  //! structure like indicator function
  //! returns if lies in boundaries
  virtual int check_boundaries(const float* q, int* res, int N)=0;
  //! for N=1
  virtual int check_boundaries(const float* q)=0;

  //! minimal value of qi, i=0...d-1
  virtual float min(int i)=0;
  //! maximal value of qi, i=0...d-1
  virtual float max(int i)=0;
  //! dq used for indicator 2
  virtual float deltaq()=0;

  //! dimension d
  virtual int dim()=0;


};

#endif // CONFIGSPACE_HPP
