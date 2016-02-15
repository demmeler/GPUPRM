#include <iostream>

#include "lib/util.hpp"

#include "lib/configspace.hpp"
#include "lib/arrayconfigspace.hpp"

using namespace std;


int main()
{
  int h=480, b=640, n=h*b;
  int* array=new int[n];
  read_file("array.bin",array, n);

  ArrayConfigspace space(array, b, h, 0.0, 1.0, 0.0, 1.0);

  printvar(space.dim());
  printvar(space.min(0));
  printvar(space.max(0));
  printvar(space.min(1));
  printvar(space.max(1));

  float qs[2]={0.5,0.4};
  float qe[2]={0.99,0.4};

  printvar(space.indicator2(&qs[0],&qe[0],0.01));

  //printarr(array,n);
}
