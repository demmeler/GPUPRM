#include <iostream>
#include <time.h>
#include <mpi.h>

//for parsing parameters
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>


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
  int num=32;
  int nbuf=2048;
  float dq=0.01;
  float distD=1.0;
  float distH=1.0;
  int seed=0;
  int blocksize=256;
  int prmversion=5;
  int workerversion=1;
  int store=1;
  int maxstorage=1024*1024;
  int maxsteps=10000;
  bool new_kernel=0;

  //! Parse inputs
  static struct option long_options[] = {
      {"num",      required_argument,        0,  'n' },
      {"nbuf", required_argument,            0,  'b' },
      {"dq",    required_argument,           0,  'q' },
      {"seed",   required_argument,          0,  's' },
      {"prmversion", required_argument,      0,  'V' },
      {"workerversion", required_argument,   0 , 'v' },
      {"store",      required_argument,      0,  'w' },
      {"maxsteps",   required_argument,      0,  'N' },
      {"new_kernel", required_argument,      0,  'k' },
      {"blocksize", required_argument,       0,  'B' },
      {"distD",     required_argument,       0,  'D' },
      {"distH",     required_argument,       0,  'H' },
      {0,           0,                       0,   0  }
  };

  char opt= 0;
  int long_index =0;
  while ((opt = getopt_long(argc, argv,"",
                            long_options, &long_index )) != -1) {
      switch (opt) {
      case 'n' : num = atoi(optarg);
          break;
      case 'b' : nbuf = atoi(optarg);
          break;
      case 'q' : dq = atof(optarg);
          break;
      case 's' : seed = atoi(optarg);
          break;
      case 'V' : prmversion = atoi(optarg);
          break;
      case 'v' : workerversion = atoi(optarg);
          break;
      case 'w' : store = atoi(optarg);
          break;
      case 'N' : maxsteps = atoi(optarg);
          break;
      case 'k': new_kernel = atoi(optarg);
          break;
      case 'B': blocksize = atoi(optarg);
          break;
      case 'D': distD = atof(optarg);
          break;
      case 'H': distH = atof(optarg);
          break;

      default:
          cout<<"bad argument: "<< opt <<endl;
          exit(EXIT_FAILURE);
      }
  }



  int confignbuf=2*nbuf;


  //srand(time(NULL));
  //srand(clock());
  //srand(rank+time(NULL));
  srand(seed+rank*10);
  int firstrand=rand();
  //printvar(firstrand);


  //! Initialization

  RobotConfigspace<ndof> space("../config1",dq, confignbuf);
  space.init_(rank,size,new_kernel);

  PRMSolver<ndof> prm(&space, distH, distD, maxstorage, blocksize);


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
    prm.process_mpi5(num,nbuf,maxsteps, seed, workerversion);
    version=3;
  }else {
    msg("Error: prmversion not valid");
    printvar(prmversion);
    return 0;
  }

  int time;
  tockval(trun,time);

  printvar(prm.rand_nodes_all);
  printvar(prm.rand_nodes_accepted);
  printvar(prm.rand_nodes_dismissed);
  printvar(prm.rand_nodes_dismissed_prob);
  printvar(prm.rand_nodes_dismissed_indicator);
  printvar(prm.rand_block_numbers);
  printvar(prm.num_nodes);
  printvar(prm.num_blocks);
  printvar(prm.num_primary_blocks);
  printvar(prm.num_secondary_blocks);
  printvar(prm.num_nodes/blocksize);

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
      printvar(distD);
      printvar(distH);
      printvar(seed);
      printvar(blocksize);
      printvar(prmversion);
      printvar(workerversion);
      printvar(store);
      printvar(new_kernel);

      ofstream file;
      file.open("stats.txt",ios::app);
      file<<version<<"\t"<<size<<"\t"<<seed<<"\t"<<time<<"\t"<<numthreadsall<<"\t"<<((long)time*1000000)/numthreadsall<<"\n";
      file.close();

  }

  MPI_Finalize();
}
