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
  //! structure like indicator function
  //! returns if lies in boundaries
  virtual int check_boundaries(const float* q, int* res, int N)=0;
  //! for N=1
  virtual int check_boundaries(const float* q, int* res)=0;

  //! minimal value of qi, i=0...d-1
  virtual float min(int i)=0;
  //! maximal value of qi, i=0...d-1
  virtual float max(int i)=0;
  //! dimension d
  virtual int dim()=0;

};

#endif // CONFIGSPACE_HPP
