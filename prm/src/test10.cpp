#include <iostream>
#include <time.h>
#include <mpi.h>

#define MPI_CODE

#include "lib/config.h"
#include "lib/util.hpp"
#include "lib/tictoc.h"

#include "lib/robotconfigspace.h"
#include "lib/polytope.h"
#include "lib/vertexlist3.h"


using namespace std;


template<int ndof>
int load_config(std::string path, Robot<ndof>* &robot, polytope* &polys,
                int* &sys, int &N, int* &from, int* &to, int& M, bool printmsg=false){

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
  check(M>0);
  from=new int[M];
  to=new int[M];
  read_file(path+"/pairs/from.bin", from, M);
  read_file(path+"/pairs/to.bin", to, M);


  if(printmsg){
    printvar(M);
    printarr(from, M);
    printarr(to,M);
    printarr(sys,N);
  }

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



int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int rank=0, size=1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //srand(time(NULL));
  //srand(clock());
  //srand(rank+time(NULL));
  int seed=(argc>=6 ? atoi(argv[5]) : 0 );
  srand(seed+rank*10);
  int firstrand=rand();

  tick(tinit);


  const int ndof=4;
  const int pi=3.1415926535897;

  float mins[ndof];
  float maxs[ndof];
  for(int i=0;i<ndof;++i){
    mins[i]=-pi; maxs[i]=1.5*pi;
  }
  float dq=(argc>=4 ? atof(argv[3]) : 0.01);//0.01;
  int confignbuf=(argc>=3 ? 2*atoi(argv[2]) : 4096);
  int numthreadsmax=(argc>=7 ? atoi(argv[6]) : 1024*1024);

  Robot<ndof>* robot;
  polytope* polys;
  int* sys;
  int N;

  int *from, *to;
  int M;

  //build_example1(robot, polys, sys, N);
  load_config<ndof>("config1",robot,polys,sys,N,from, to, M, false);

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


  float D=(argc>=5 ? atof(argv[4]) : 1.0);
  vertexlist<ndof> prm(D,D,&space);

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


  tock(tinit);
  tick(trun);

  int num=(argc>=2 ? atoi(argv[1]) : 32 );//128;
  int nbuf=(argc>=3 ? atoi(argv[2]) : 2048);//2048;
  int maxsteps=100000;
  //prm.process_mpi(num,nbuf,maxsteps);
  prm.process_mpi2(num,nbuf,maxsteps, seed);

  tock(trun);

  //prm.print();

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
      tick(twrite);

      prm.store_graphs("prmoutput");

      tock(twrite);

      printvar(rand());
      printvar(firstrand);

      printvar(num);
      printvar(nbuf);
      printvar(dq);
      printvar(D);
      printvar(seed);
      printvar(numthreadsmax);
  }

  MPI_Finalize();
}
