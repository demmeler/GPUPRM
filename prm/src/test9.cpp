#include <iostream>
#include <time.h>
#include "lib/config.h"

#include "lib/util.hpp"
#include "lib/tictoc.h"

#include "lib/configspace.hpp"
#include "lib/robotconfigspace.h"

#include "lib/polytope.h"

#include "lib/vertexlist3.h"

using namespace std;


template<int ndof>
int load_config(std::string path, Robot<ndof>* &robot, polytope* &polys, int* &sys, int &N, int* &from, int* &to, int& M, bool printmsg=false){

    //! DH params


  int ndof_loaded=0;
  read_file(path+"/ndof.bin",&ndof_loaded,1);
  if(ndof_loaded!=ndof){
    msg("error: wrong ndof");
    return -1;
  }

  robot=new Robot<ndof>;

  read_file(path+"/dh/a.bin",&(robot->a[0]),ndof);
  read_file(path+"/dh/alpha.bin",&(robot->alpha[0]),ndof);
  read_file(path+"/dh/q.bin",&(robot->q[0]),ndof);
  read_file(path+"/dh/d.bin",&(robot->d[0]),ndof);
  read_file(path+"/dh/types.bin",(int*)&(robot->types[0]),ndof);

  if(printmsg){
      printvar(ndof);
      printarr(robot->a,ndof);
      printarr(robot->alpha,ndof);
      printarr(robot->q,ndof);
      printarr(robot->d,ndof);
      printarr(robot->types,ndof);
  }

    //! Polytopes

  read_file(path+"/polys/N.bin",&N,1);
  if(N<=0){
    msg("error: N<=0");
    return -1;
  }

  sys=new int[N];
  read_file(path+"/polys/sys.bin",sys,N);

  read_file(path+"/pairs/M.bin",&M, 1);   //--> check machen ob file existiert
  printvar(M);
  check(M>0);
  from=new int[M];
  to=new int[M];
  read_file(path+"/pairs/from.bin", from, M);
  read_file(path+"/pairs/to.bin", to, M);
  printarr(from, M);
  printarr(to,M);

  if(printmsg) printarr(sys,N);

  polys=new polytope[N];
  for(int i=0;i<N;++i){
    std::stringstream polypathstream;
    polypathstream<<path<<"/polys/poly"<<i;
    std::string polypath=polypathstream.str();
    int size[2];
    read_file(polypath+"/size.bin",&(size[0]),2);
    int n,m;
    n=polys[i].n=size[0];
    m=polys[i].m=size[1];
    polys[i].dsp=new int[n];
    polys[i].cnt=new int[n];
    polys[i].dest=new int[m];
    polys[i].vertices=new polytope::vec4[4*n];
    float* vertices=new float[3*n];
    read_file(polypath+"/vertices.bin",vertices,3*n);
    for(int j=0;j<n;++j){
        polys[i].vertices[j].x=vertices[3*j];
        polys[i].vertices[j].y=vertices[3*j+1];
        polys[i].vertices[j].z=vertices[3*j+2];
        polys[i].vertices[j].w=1.0;
    }
    read_file(polypath+"/dsp.bin",polys[i].dsp,n);
    read_file(polypath+"/cnt.bin",polys[i].cnt,n);
    read_file(polypath+"/dest.bin",polys[i].dest,m);


    if(printmsg){
      msg("-----");
      printvar(i);
      //geo4::trafo4 t(0.0,0.0,0.0,0.0);
      //p4print(polys[i],t);

      printvar(n);
      printvar(m);
      printarr(polys[i].dsp,n);
      printarr(polys[i].cnt,n);
      printarr(polys[i].dest,m);
    }

  }
}



int main()
{

  srand(time(NULL));
  //srand(clock());

  tick(t0);


  const int ndof=4;
  const int pi=3.1415926535897;

  float mins[ndof];
  float maxs[ndof];
  for(int i=0;i<ndof;++i){
    mins[i]=-pi; maxs[i]=1.5*pi;
  }
  float dq=0.01;
  int confignbuf=500;
  int numthreadsmax=1024*1024;

  Robot<ndof>* robot;
  polytope* polys;
  int* sys;
  int N;

  int *from, *to;
  int M;

  //build_example1(robot, polys, sys, N);
  load_config<ndof>("config1",robot,polys,sys,N,from, to, M, true);

  RobotConfigspace<ndof> space(robot,
                            polys, sys, N,
                            from, to, M,
                            mins, maxs, dq,
                            confignbuf, numthreadsmax);

  space.init();

  printvar(space.dim());
  for(int dof=0;dof<ndof;++dof){
    printvar(dof);
    printvar(space.min(dof));
    printvar(space.max(dof));
  }
  printvar(space.deltaq());

  vertexlist<ndof> prm(0.5,0.5,&space);

  float qstart[ndof];
  float qend[ndof];
  for(int dof=0;dof<ndof;++dof){
    qstart[dof]=0.0;
    qend[dof]=pi;
  }

  printvar(space.indicator(&qstart[0]));
  printvar(space.indicator(&qend[0]));

  prm.init(&qstart[0],&qend[0]);

  //prm.print();

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

  for(int i=0;i<10000;++i){
    int flag=prm.processing_step(qnew,num,qstartlist,qendlist,resbuf,nbuf,offset);
    if(flag==1){
      msg("connection found");
      printvar(i);
      break;
    }else if(i%50==0){
      printvar(i);
    }
  }

  tock(t1);

  //prm.print();

  tick(t2);

  int ret1=system("rm -rf prmoutput");
  int ret2=system("mkdir prmoutput");
  prm.store_graphs("prmoutput");

  tock(t2);

  //printarr(array,3000);

  //printarr(array,n);
  printvar(clock());
  printvar(rand());
}
