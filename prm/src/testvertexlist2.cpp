



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

  srand(seed);

  msg("*********");

  tick(t1);

  float D=0.5;
  float H=1.0;//0.001;

  printvar(D);
  printvar(H);

  const int ndof=2;

  vertexlist<ndof> vlist(H,D);

  int nbuf=100;
  float* qlist=new float[ndof*nbuf]();
  float qtemp[ndof];

  for(int i=0;i<10000;++i){
    for(int i=0;i<ndof;++i){
      qtemp[i]=-5.0+10.0*(float)rand()/RAND_MAX;
    }
    vlist.insert(qtemp);
  }


  //tock(t1);
  tick(t2);

  int sum=0;

  for(int i=0;i<ndof;++i){
    qtemp[i]=-5.0+10.0*(float)rand()/RAND_MAX;
  }

  printvar(vlist.get_near_vertices(qtemp,qlist,nbuf,nbuf));

  printvar(sum);
  printarr(qtemp,2);
  printarr(qlist,nbuf*ndof);

  tock(t2);



  return 0;
}

