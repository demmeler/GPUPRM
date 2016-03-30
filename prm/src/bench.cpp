
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


  //! Parameters

  int num=(argc>=2 ? atoi(argv[1]) : 32 );
  int nbuf=(argc>=3 ? atoi(argv[2]) : 2048);
  int confignbuf=(argc>=3 ? 2*atoi(argv[2]) : 4096);
  float dq=(argc>=4 ? atof(argv[3]) : 0.01);
  float D=(argc>=5 ? atof(argv[4]) : 1.0);
  int seed=(argc>=6 ? atoi(argv[5]) : 0 );
  int blocksize=(argc>=7 ? atoi(argv[6]) : 256);
  int prmversion=(argc>=8 ? atoi(argv[7]) : 5 );
  int maxstorage=1024*1024;
  int maxsteps=100000;



  srand(seed+rank*10);
  int firstrand=rand();


  //! Initialization

  RobotConfigspace<ndof> space("config1",dq, confignbuf);
  space.init(rank,size);

  PRMSolver<ndof> prm(&space, D, D, maxstorage, blocksize);


  float qstart[ndof];
  float qend[ndof];
  for(int dof=0;dof<ndof;++dof){
    qstart[dof]=0.0;
    qend[dof]=3;
  }

  prm.init(&qstart[0],&qend[0]);


  //! run

  tick(trun);


  int version=-1;

  if(prmversion==1){
    prm.process_mpi(num,nbuf,maxsteps);
    version=1;
  }if(prmversion==2){
    prm.process_mpi2(num,nbuf,maxsteps, seed);
    version=-1;
  }else if(prmversion==3){
    prm.process_mpi3(num,nbuf,maxsteps, seed);
    version=-1;
  }else if(prmversion==4){
    prm.process_mpi4(num,nbuf,maxsteps, seed);
    version=2;
  }else if(prmversion==5){
    prm.process_mpi5(num,nbuf,maxsteps, seed);
    version=3;
  }else {
    msg("Error: prmprocess not valid");
    return 0;
  }

  tock(trun);

  //! write output

  MPI_Barrier(MPI_COMM_WORLD);
  int numthreads=space.get_numthreads_all();
  int numthreadsall;
  MPI_Reduce(&numthreads, &numthreadsall,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank==0){
    ofstream file;
    file.open("results.txt",ios::app);
    file<<size<<"\t"<<prmversion<<"\n";
    file.close();
  }

  MPI_Finalize();
}
