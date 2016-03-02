

#include <iostream>

#define SILENT
#define NO_IO

#include "lib/util.hpp"

#include "lib/configspace.hpp"
#include "lib/robotconfigspace.h"
#include "lib/geo4.h"
using namespace geo4;


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
  generate_quader(polys[0], 1.0, 1.0, 1.0);
  transform(polys[0],t0);
  sys[0]=0;

  trafo4 t1(0.0, 0.0, 0.0, 0.0);
  generate_quader(polys[1], 1.0, 1.0, 1.0);
  transform(polys[1],t1);
  sys[1]=2;

  //trafo4 t2(0.0, -pi/2.0, 0.0, 0.0);
  //generate_quader(P[2], 1.0, 1.0, 1.0);
  //transform(P[2],t2);

}

template<int ndof>
int load_config(std::string path, Robot<ndof>* &robot, polytope4* &polys, int* &sys, int &N, int* &from, int* &to, int& M, bool printmsg=false){

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

  read_file(path+"/polys/pairs/M.bin",&M, 1);   //--> check machen ob file existiert
  printvar(M);
  check(M>0);
  from=new int[M];
  to=new int[M];
  read_file(path+"/polys/pairs/from.bin", from, M);
  read_file(path+"/polys/pairs/to.bin", to, M);

  if(printmsg) printarr(sys,N);

  polys=new polytope4[N];
  for(int i=0;i<N;++i){
    std::string polypath=path+"/polys/poly"+std::to_string(i);
    int size[2];
    read_file(polypath+"/size.bin",&(size[0]),2);
    int n,m;
    n=polys[i].n=size[0];
    m=polys[i].m=size[1];
    polys[i].dsp=new int[n];
    polys[i].cnt=new int[n];
    polys[i].dest=new int[m];
    polys[i].vertices=new float4[4*n];
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
      trafo4 t(0.0,0.0,0.0,0.0);
      p4print(polys[i],t);

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

  srand(clock());

  tick(t0);

  float mins[2]={-5.0,-5.0};
  float maxs[2]={5.0,5.0};
  float dq=0.01;
  int confignbuf=500;
  int numthreadsmax=1024*1024;

  Robot<2>* robot;
  polytope4* polys;
  int* sys;
  int N;

  int *from, *to;
  int M;

  //build_example1(robot, polys, sys, N);
  load_config<2>("config1",robot,polys,sys,N,from, to, M, true);

  RobotConfigspace<2> space(robot,
                            polys, sys, N,
                            from, to, M,
                            mins, maxs, dq,
                            confignbuf, numthreadsmax);

  space.init();

  printvar(space.dim());
  printvar(space.min(0));
  printvar(space.max(0));
  printvar(space.min(1));
  printvar(space.max(1));
  printvar(space.deltaq());

#if 0

  Kinematics<2> kin(robot);

  float qtest[2]={3.1415926/2.0,0.0};
  kin.calculate(&qtest[0],1);


  //p4print(polys[0],kin.trafos[0]);
  //p4print(polys[1],kin.trafos[1]);



  float qs[6]={0.0, 0.0, 0.0,
               0.0, 0.0, 0.0};
  float qe[6]={0.0, 1.0, 1.0,
               4.0, 0.0, 0.5};
  int res[3]={7,7,7};

  space.indicator2(&qs[0],&qe[0],&res[0],3,3);
  printarr(res,3);

#endif

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
    }else if(i%50==0){
      printvar(i);
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
