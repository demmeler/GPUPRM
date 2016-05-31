#include <iostream>
#include <time.h>
#include <mpi.h>

#define MPI_CODE

#include "lib/config.h"
#include "lib/util.hpp"
#include "lib/tictoc.h"

#include "lib/robotspace/robotconfigspace.h"
#include "lib/prm/prmsolver.h"


using namespace std;


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int rank=0, size=1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  const int ndof=4;


  tick(tinit);

  //! Parameters

  int num=(argc>=10 ? atoi(argv[9]) : 32 );
  int nbuf=(argc>=3 ? atoi(argv[2]) : 2048);
  int confignbuf=(argc>=3 ? 2*atoi(argv[2]) : 4096);
  float dq=(argc>=4 ? atof(argv[3]) : 0.01);
  float D=(argc>=5 ? atof(argv[4]) : 1.0);
  int seed=(argc>=6 ? atoi(argv[5]) : 0 );
  int blocksize=(argc>=7 ? atoi(argv[6]) : 256);
  int prmversion=(argc>=8 ? atoi(argv[7]) : 5 );
  int store=(argc>=9 ? atoi(argv[8]) : 1 );
  int maxstorage=1024*1024;
  int maxsteps=10000;//(argc>=2 ? atoi(argv[1]) : 100000 );

  bool new_kernel=(argc>=2 ? atoi(argv[1]) : 1 );

  //srand(time(NULL));
  //srand(clock());
  //srand(rank+time(NULL));
  srand(seed+rank*10);
  int firstrand=rand();
  //printvar(firstrand);


  //! Initialization

  RobotConfigspace<ndof> space("../config1",dq, confignbuf);
  space.init_(rank,size,new_kernel);

  PRMSolver<ndof> prm(&space, D, D, maxstorage, blocksize);


  float qstart[ndof];
  float qend[ndof];
  for(int dof=0;dof<ndof;++dof){
    qstart[dof]=0.0;
    qend[dof]=3;
  }

  if(rank==0){
    assert(space.indicator(&qstart[0])==0);
    assert(space.indicator(&qend[0])==0);
  }

  prm.init(&qstart[0],&qend[0]);

  tock(tinit);


  //! run

  tick(trun);

  int version=-1;
  if(prmversion==1){
    prm.process_mpi(num,nbuf,maxsteps);
    version=1;
  }else if(prmversion==2){
    prm.process_mpi2(num,nbuf,maxsteps, seed);
    version=-1;
  }else if(prmversion==3){
    prm.process_mpi3(num,nbuf,maxsteps, seed);
    version=-1;
  }else if(prmversion==4){
    prm.process_mpi4(num,nbuf,maxsteps, seed);
    version=2;
  }else if(prmversion==5){
    prm.process_mpi5(num,nbuf,maxsteps, seed,1); //workerversion 1
    version=3;
  }else if(prmversion==6){
    prm.process_mpi5(num,nbuf,maxsteps, seed,2); //workerversion 2
    version=4;
  }else {
    msg("Error: prmversion not valid");
    printvar(prmversion);
    return 0;
  }

  int time;
  tockval(trun,time);


  //! write output

  //prm.print();

  MPI_Barrier(MPI_COMM_WORLD);
  int numthreads=space.get_numthreads_all();
  int numthreadsall;
  MPI_Reduce(&numthreads, &numthreadsall,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank==0){

      if(store){
        tick(twrite);
        prm.store_results("../prmoutput");
        tock(twrite);
      }

      printvar(numthreadsall);

      printvar(rand());
      printvar(firstrand);

      printvar(num);
      printvar(nbuf);
      printvar(dq);
      printvar(D);
      printvar(seed);
      printvar(blocksize);
      printvar(prmversion);
      printvar(store);
      printvar(new_kernel);

      ofstream file;
      file.open("stats.txt",ios::app);
      file<<version<<"\t"<<size<<"\t"<<seed<<"\t"<<time<<"\t"<<numthreadsall<<"\t"<<((long)time*1000000)/numthreadsall<<"\n";
      file.close();

  }

  MPI_Finalize();
}
