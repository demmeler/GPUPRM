#ifndef CONFIGSPACE_HPP
#define CONFIGSPACE_HPP

template<int d>
class Configspace
{
public:
  virtual Configspace();
  virtual ~Configspace();

  virtual int indicator(double*q, int N);

};

#endif // CONFIGSPACE_HPP
