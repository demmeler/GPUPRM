


#include <iostream>

#include "lib/util.hpp"

#include "lib/vertexlist.h"

#include <chrono>
#include <cstdlib>
#include <ctime>


using namespace std;

int pow(int x, int n){
    if(n>1)return x*pow(x,n-1);
    else if(n==1) return x;
    else if(n==0) return 1;
    else return 0;
}

#define tick(t1) std::chrono::high_resolution_clock::time_point t1=std::chrono::high_resolution_clock::now();
#define tock(t1,t2) std::chrono::high_resolution_clock::time_point t2=std::chrono::high_resolution_clock::now(); \
                  {int time=std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();              \
                  std::cout << #t1 << " to "<< #t2 << ": \t" << time << "ms" << endl; }

int main()
{

  tick(t1);

  int x=0;
  for(int i=0;i<100000000;++i){
      x+=i*i;
  }
  printvar(x);

  tock(t1,t2);


  return 0;
}

