
#include <iostream>

#define SILENT

#include "lib/util.hpp"

#include "lib/configspace.hpp"
#include "lib/robotconfigspace.h"
#include "lib/geo4.h"
using namespace geo4;

#define NO_IO

#include "lib/vertexlist3.h"

using namespace std;
using namespace collision4;


//void build_example1(Robot<1>* &robot, polytope4data* &polydata)
void build_example1(Robot<2>* &robot, polytope4* &polys, int* &sys, int &N){
  robot=new Robot<2>;

  robot->a[0]=1.1;
  robot->alpha[0]=0.0;
  robot->q[0]=0.0;
  robot->d[0]=0.0;
  robot->types[0]=prismatic;

  robot->a[1]=0.0;
  robot->alpha[1]=0.0;
  robot->q[1]=0.0;
  robot->d[1]=0.0;
  robot->types[1]=rotational;


  polys=new polytope4[2];
  sys=new int[2];
  N=2;

  trafo4 t0(0.0, 0.0, 0.0, 0.0);
  generate_simplex(polys[0], 1.0, 1.0, 1.0);
  transform(polys[0],t0);
  sys[0]=0;

  trafo4 t1(0.0, 0.0, 0.0, 0.0);
  generate_simplex(polys[1], 1.0, 1.0, 1.0);
  transform(polys[1],t1);
  sys[1]=2;

  //trafo4 t2(0.0, -pi/2.0, 0.0, 0.0);
  //generate_quader(P[2], 1.0, 1.0, 1.0);
  //transform(P[2],t2);

}


int main()
{

  srand(clock());

  tick(t0);

  float mins[2]={-1.1,-0.1};
  float maxs[2]={1.5,4.1};
  float dq=0.01;
  int confignbuf=500;
  int numthreadsmax=1024*1024;

  Robot<2>* robot;
  polytope4* polys;
  int* sys;
  int N;
  build_example1(robot, polys, sys, N);

  RobotConfigspace<2> space(robot,
                            polys, sys, N,
                            mins, maxs, dq,
                            confignbuf, numthreadsmax);

  space.init();

  printvar(space.dim());
  printvar(space.min(0));
  printvar(space.max(0));
  printvar(space.min(1));
  printvar(space.max(1));
  printvar(space.deltaq());

  Kinematics<2> kin(robot);

  float qtest[2]={3.1415926/2.0,0.0};
  kin.calculate(&qtest[0],1);


  p4print(polys[0],kin.trafos[0]);
  p4print(polys[1],kin.trafos[1]);



  float qs[6]={0.0, 0.0, 0.0,
               0.0, 0.0, 0.0};
  float qe[6]={0.0, 1.0, 1.0,
               4.0, 0.0, 0.5};
  int res[3]={7,7,7};

  space.indicator2(&qs[0],&qe[0],&res[0],1,3);
  printarr(res,3);



  vertexlist<2> prm(0.5,0.5,&space);

  float qstart[2]={0.0,0.0};
  float qend[2]={0.0,4.0};

  printvar(space.indicator(&qstart[0]));
  printvar(space.indicator(&qend[0]));

  prm.init(&qstart[0],&qend[0]);

  //prm.print();

  int ndof=2;
  int num=4;
  float *qnew=new float[ndof*num];
  int nbuf=500;
  float *qstartlist=new float[ndof*nbuf];
  float *qendlist=new float[ndof*nbuf];
  int *resbuf=new int[nbuf];
  for(int i=0;i<nbuf;++i){resbuf[i]=27;}
  int offset=nbuf;

  tock(t0);
  tick(t1);

  for(int i=0;i<1000;++i){
    int flag=prm.processing_step(qnew,num,qstartlist,qendlist,resbuf,nbuf,offset);
    if(flag==1){
      msg("connection found");
      printvar(i);
      break;
    }
  }

  tock(t1);

  //prm.print();

  tick(t2);

  system("rm -rf prmoutput");
  system("mkdir prmoutput");
  prm.store_graphs("prmoutput");

  tock(t2);

  //printarr(array,3000);

  //printarr(array,n);
}
