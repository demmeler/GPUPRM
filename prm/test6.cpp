
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
  printvar(space.deltaq());

  float qs[4]={0.8,0.1,
               0.4,0.4};
  float qe[6]={0.99,0.9,0.2,
               0.4,0.4,0.4};
  int posqe[2]={0,1};
  int numqe[2]={1,2};
  int res[2]={7,7};
  space.indicator2(&qs[0],2,&qe[0],&res[0],&posqe[0],&numqe[0],3);
  //space.indicator2(&qs[0],&qe[0],&res[0],2);
  printarr(res,3);

  //printarr(array,3000);

  //printarr(array,n);
}
