


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

  tick(t1);

  float D=0.5;

  vertexlist<2> vlist(D);


  int nbuf=10000;
  float* qlist=new float[2*nbuf]();

  float qtemp[2];

  for(int i=0;i<1000000;++i){
    qtemp[0]=-5.0+10.0*(float)rand()/RAND_MAX;
    qtemp[1]=-5.0+10.0*(float)rand()/RAND_MAX;
    if(i%1000==0)
    //printvar(vlist.get_near_vertices(qtemp,qlist,nbuf,nbuf));
    vlist.insert(qtemp);
  }


  tock(t1,t2);
  tick(t3);

  int sum=0;

  for(int i=0;i<100000;++i){
      qtemp[0]=-5.0+10.0*(float)rand()/RAND_MAX;
      qtemp[1]=-5.0+10.0*(float)rand()/RAND_MAX;
      sum+=vlist.get_near_vertices(qtemp,qlist,nbuf,nbuf);
      //printarr(qlist,2*nbuf);
  }

  printvar(sum);

  tock(t3,t4);


  return 0;
}

