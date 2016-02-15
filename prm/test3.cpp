#include <iostream>

#include "lib/util.hpp"

#include "lib/configspace.hpp"
#include "lib/arrayconfigspace.hpp"

using namespace std;


int main()
{
  int h=3, b=4, n=h*b;
  int* array=new int[n];
  read_file("array.bin",array, n);


  printarr(array,n);
}
