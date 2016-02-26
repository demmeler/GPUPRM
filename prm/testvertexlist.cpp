


#include <iostream>

#include "lib/util.hpp"

#include "lib/vertexlist2.h"

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

  for(int l=1;l<200;++l){
  srand(0);

  msg("*********");

  tick(t1);

  float D=0.01;
  float H=l*0.01;//0.001;

  printvar(D);
  printvar(H);

  const int ndof=8;

  vertexlist<ndof> vlist(H,D);

  int nbuf=100000;
  float* qlist=new float[ndof*nbuf]();
  float qtemp[ndof];

  for(int i=0;i<100000;++i){
    for(int i=0;i<ndof;++i){
      qtemp[i]=-5.0+10.0*(float)rand()/RAND_MAX;
    }
#if 0
    int x=vlist.get_near_vertices(qtemp,qlist,nbuf,nbuf);
    if(i%10000==0){
      printvar(x);
      tock(t1);
    }
#endif
    vlist.insert(qtemp);
  }


  //tock(t1);
  tick(t2);

  int sum=0;

  for(int i=0;i<100000;++i){
    for(int i=0;i<ndof;++i){
      qtemp[i]=-5.0+10.0*(float)rand()/RAND_MAX;
    }

    sum+=vlist.get_near_vertices(qtemp,qlist,nbuf,nbuf);
    //printarr(qlist,2*nbuf);
#if 0
    if(i%50==0){
      printvar(i);
      printvar(sum);
      tock(t2);
    }
#endif
  }

  //printvar(sum);
  tock(t2);

  }

  return 0;
}

