
#include <iostream>

#include "lib/util.hpp"

#include "lib/configspace.hpp"
#include "lib/arrayconfigspace.hpp"

#include "lib/vertexlist3.h"

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


  vertexlist<2> prm(0.1,0.1,&space);

  float qstart[2]={0.1,0.4};
  float qend[2]={0.9,0.4};

  printvar(space.indicator(&qstart[0]));
  printvar(space.indicator(&qend[0]));

  prm.init(&qstart[0],&qend[0]);

  prm.print();

  int ndof=2;
  int num=2;
  float *qnew=new float[ndof*num];
  int nbuf=4;
  float *qlist=new float[ndof*nbuf];
  int *resbuf=new int[nbuf];
  int offset=nbuf;

  prm.processing_step(qnew,num,qlist,resbuf,nbuf,offset);



  //printarr(array,3000);

  //printarr(array,n);
}
