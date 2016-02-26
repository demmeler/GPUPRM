


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



int main()
{

  int seed=0;

  for(int l=1;l<100;++l){
  srand(0);

  msg("*********");

  tick(t1);

  float D=0.1;
  float H=l*0.0005;

  printvar(D);
  printvar(H);

  vertexlist<2> vlist(H,D);

  int nbuf=10000;
  float* qlist=new float[2*nbuf]();
  float qtemp[2];

  for(int i=0;i<10000;++i){
    qtemp[0]=-5.0+10.0*(float)rand()/RAND_MAX;
    qtemp[1]=-5.0+10.0*(float)rand()/RAND_MAX;
#if 0
    int x=vlist.get_near_vertices(qtemp,qlist,nbuf,nbuf);
    if(i%10000==0){
      printvar(x);
      tock(t1);
    }
#endif
    vlist.insert(qtemp);
  }


  tock(t1);
  tick(t2);

  int sum=0;

  for(int i=0;i<100000;++i){
    qtemp[0]=-5.0+10.0*(float)rand()/RAND_MAX;
    qtemp[1]=-5.0+10.0*(float)rand()/RAND_MAX;
    sum+=vlist.get_near_vertices(qtemp,qlist,nbuf,nbuf);
    //printarr(qlist,2*nbuf);
#if 0
    if(i%5000==0){
      printvar(i);
      printvar(sum);
      tock(t2);
    }
#endif
  }

  printvar(sum);
  tock(t2);

  }

  return 0;
}

