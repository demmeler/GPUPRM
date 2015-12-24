#ifndef CONFIGSPACE_HPP
#define CONFIGSPACE_HPP

class Configspace
{
public:
  virtual Configspace(int d_){d=d_;}
  virtual ~Configspace();

  virtual int indicator(double*q, int N);

private:
  int d;
};

#endif // CONFIGSPACE_HPP
