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



int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int rank=0, size=1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  const int ndof=4;


  tick(tinit);

  //! Parameters

  int num=(argc>=2 ? atoi(argv[1]) : 32 );
  int nbuf=(argc>=3 ? atoi(argv[2]) : 2048);
  int confignbuf=(argc>=3 ? 2*atoi(argv[2]) : 4096);
  float dq=(argc>=4 ? atof(argv[3]) : 0.01);
  float D=(argc>=5 ? atof(argv[4]) : 1.0);
  int seed=(argc>=6 ? atoi(argv[5]) : 0 );
  int maxsteps=100000;

  //srand(time(NULL));
  //srand(clock());
  //srand(rank+time(NULL));
  srand(seed+rank*10);
  int firstrand=rand();


  //! Initialization

  RobotConfigspace<ndof> space("config1",dq, confignbuf);
  space.init(rank,size);

  vertexlist<ndof> prm(D,D,&space);


  float qstart[ndof];
  float qend[ndof];
  for(int dof=0;dof<ndof;++dof){
    qstart[dof]=0.0;
    qend[dof]=3;
  }

  if(rank==0){
    printvar(space.indicator(&qstart[0]));
    printvar(space.indicator(&qend[0]));
  }

  prm.init(&qstart[0],&qend[0]);

  //prm.print();

  tock(tinit);


  //! run

  tick(trun);


  //prm.process_mpi(num,nbuf,maxsteps);
  prm.process_mpi3(num,nbuf,maxsteps, seed);

  tock(trun);


  //! write output

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
  }

  MPI_Finalize();
}
