#define MPI_CODE
#include "prmsolver.h"

#include <cmath>
#include <cfloat>
#include <stdlib.h>
#include <queue>

#include <util.hpp>
#include <tictoc.h>


#define NO_IO


    //! *******************
    //! *      class      *
    //! *******************

  template<int ndof>
  PRMSolver<ndof>::PRMSolver(Configspace<ndof> *space_, float H_, float D_, int N, int blocksize)
    :N(N),
    blocksize(blocksize)
  {

    int connectivity=20; //pi*daumen wert fuer reservierung

    graphl.qstorage.resize(ndof*N);
    graphl.surrnum.resize(N);
    graphl.edgelists.resize(N);
    graphl.edgeweights.resize(N);
    graphl.blocks.resize(N/blocksize);
    graphl.from.reserve(N*connectivity);
    graphl.to.reserve(N*connectivity);
    graphl.weights.reserve(N*connectivity);

    graphl.newblockpos=0;
    graphl.blocknum=0;

    graphr.qstorage.resize(ndof*N);
    graphr.surrnum.resize(N);
    graphr.edgelists.resize(N);
    graphr.edgeweights.resize(N);
    graphr.blocks.resize(N/blocksize);
    graphr.from.reserve(N*connectivity);
    graphr.to.reserve(N*connectivity);
    graphr.weights.reserve(N*connectivity);
    graphr.newblockpos=0;
    graphr.blocknum=0;


    block b;
    b.next=0x0;
    b.num=0;
    b.pos=0;
    graphl.blocks.assign(N,b);
    graphr.blocks.assign(N,b);

    connection_found=true;

    H=H_;
    D=D_;
    D2=D*D;
    factor=1.0/H;

    space=space_;
  }

  template<int ndof>
  PRMSolver<ndof>::~PRMSolver(){
    //nothing to do
  }

  template<int ndof>
  int PRMSolver<ndof>::init(const float* qstart, const float* qend){
    i0l=insert(qstart,graphl);
    i0r=insert(qend,graphr);
    return 0;
  }























  //! ***************************
  //! *                         *
  //! *   processing routines   *
  //! *                         *
  //! ***************************



  //! *******************
  //! *                 *
  //! *    process 1    *
  //! *                 *
  //! *******************

  template<int ndof>
  int PRMSolver<ndof>::process_mpi(int num, const int nbuf, const int maxsteps){
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printvar(rank);
    printvar(size);
    printvar(num);
    assert(num%size==0);

    int cnt=num/size;
    int dsp=cnt*rank;


    float *qnew=new float[cnt*ndof];
    float *qnewall=new float[num*ndof];
    float *qstartlist=new float[ndof*nbuf];
    float *qendlist=new float[ndof*nbuf];
    int *resbuf=new int[nbuf];
    int offset=nbuf;

    int *leftconn=new int[cnt];
    int *rightconn=new int[cnt];
    int *leftconnall=new int[num];
    int *rightconnall=new int[num];

    int *surrnumnew=new int[cnt];
    int *surrnumnewall=new int[num];

    int *positionsall=new int[num];

    for(int i=0;i<maxsteps;++i){
      int flag=processing_step(rank,size,
                               qnew, qnewall,num,dsp,cnt,
                               surrnumnew,surrnumnewall,positionsall,
                               leftconn, leftconnall,rightconn, rightconnall,
                               qstartlist,qendlist,resbuf,nbuf,offset);
      if(flag==1){
        msg("connection found");
        printvar(i);
        break;
      }else if(i%50==0){
        printvar(i);
      }
    }

    delete qnew;
    delete qnewall;
    delete qstartlist;
    delete qendlist;
    delete resbuf;
    delete leftconn;
    delete rightconn;
    delete leftconnall;
    delete rightconnall;
    delete surrnumnew;
    delete surrnumnewall;
    delete positionsall;

    return 0;
  }


  //! get list of all vertices nearer than D
  //! qlist: buffer to store neighbour candidates, struct of arrays
  //! offset: i-th component of k-th vec in qlist is qlist[i*offset+k]
  //! nbuf: length(qlist)
  //! D: maximal distance
  //! return: 0 - no errors, 1 - connection found
  template<int ndof>
  int PRMSolver<ndof>::processing_step(const int mpirank, const int mpisize,
                             float* qnew, float* qnewall, const int num, const int dsp, const int cnt,
                             int *surrnumnew, int *surrnumnewall, int *positionsall,
                             int* leftconn, int *leftconnall, int *rightconn, int *rightconnall,
                             float* qstart, float* qend, int* resbuf, const int nbuf, const int offset){
    //float* qnew=qnewall+ndof*dsp;

    //!
    //! create nodes qnew randomly
    //!

    //!from left graph
    for(int j=0;j<cnt/2;++j){
      bool dismiss;
      do{
        //!choose random block
        int k;
        block *b;
        do{
            k=rand()%graphl.blocknum;
            b=&(graphl.blocks[k]);
        }while(b->num==0);
        //!chose random vertex
        int m=(b->pos+rand()%b->num);
        int l;
        int x=1+graphl.surrnum[m];
        int prob=RAND_MAX/(x*x*x);
        if(rand()>prob){
          dismiss=true;
        }else{
          l=ndof*m;
          for(int i=0;i<ndof;++i){
            qnew[ndof*j+i]=graphl.qstorage[l+i]-D+2*D*((float)rand()/RAND_MAX);
          }

          dismiss=space->indicator(&qnew[ndof*j])!=0;
        }
#ifndef NO_IO
        printvar(l);
        printarr(&qnew[ndof*j],ndof);
        printvar(dismiss);
#endif
      }while(dismiss);
    }
    //!from right graph
    for(int j=cnt/2;j<cnt;++j){
      bool dismiss;
      do{
        //!choose random block
        int k;
        block *b;
        do{
            k=rand()%graphr.blocknum;
            b=&(graphr.blocks[k]);
        }while(b->num==0);
        //!chose random vertex
        int m=(b->pos+rand()%b->num);
        int l;
        int x=1+graphr.surrnum[m];
        int prob=RAND_MAX/(x*x*x);
        if(rand()>prob){
          dismiss=true;
        }else{
          l=ndof*m;
          for(int i=0;i<ndof;++i){
            qnew[ndof*j+i]=graphr.qstorage[l+i]-D+2*D*((float)rand()/RAND_MAX);
          }
          dismiss=space->indicator(&qnew[ndof*j])!=0;
        }
#ifndef NO_IO
        printvar(l);
        printarr(&qnew[ndof*j],ndof);
        printvar(dismiss);
#endif
      }while(dismiss);
    }

    // --> start sending and recieving new nodes (async)

    //MPI_Allgather(snd.data(), snd.size(), MPI_DOUBLE, rcv.data(), rcv.size()/size, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Request qnewrequest;
    MPI_Iallgather(qnew, ndof*cnt, MPI_FLOAT, qnewall, ndof*cnt, MPI_FLOAT, MPI_COMM_WORLD,&qnewrequest);



    //!
    //! make list of potential neighbours for all new nodes
    //!


    //!posqlist[i]: position in qlist buffer, at which neighbours of qnew[i] start
    int posqlist[cnt];
    int numqlist[cnt];
    int numqlistleft[cnt];
    int Nqlist;             //sum(numqlist)

    int *poslist=new int[nbuf];
    float *distlist=new float[nbuf];
    int index=0;
    float *qstartp=qstart;
    float *qendp=qend;
    int *poslistp=poslist;
    float *distlistp=distlist;
    int nbufrest=nbuf;

    for(int j=0;j<cnt;++j){
      posqlist[j]=index;

      int writtenleft=get_near_vertices(&qnew[ndof*j],1,qendp,poslistp,distlistp,nbufrest,offset,graphl);

      for(int i=0;i<ndof;++i)
      for(int k=0;k<writtenleft;++k){
        qstartp[k+offset*i]=qnew[ndof*j+i];
      }
      numqlistleft[j]=writtenleft;
      index+=writtenleft;

      if(writtenleft>=nbufrest){
        for(int l=j+1;l<cnt;++l){
          posqlist[l]=nbuf;
          numqlist[l]=0;
          numqlistleft[l]=0;
          msg("buffer full");
        }
        numqlist[j]=writtenleft;
        break;
      }

      qstartp+=writtenleft;
      qendp+=writtenleft;
      poslistp+=writtenleft;
      distlistp+=writtenleft;
      nbufrest-=writtenleft;

      int writtenright=get_near_vertices(&qnew[ndof*j],1,qendp,poslistp,distlistp,nbufrest,offset,graphr);

      for(int i=0;i<ndof;++i)
      for(int k=0;k<writtenright;++k){
        qstartp[k+offset*i]=qnew[ndof*j+i];
      }
      numqlist[j]=writtenleft+writtenright;
      index+=writtenright;

      if(writtenright>=nbufrest){
        for(int l=j+1;l<cnt;++l){
          posqlist[l]=nbuf;
          numqlist[l]=0;
          numqlistleft[l]=0;
          msg("buffer full");
        }
        break;
      }

      qstartp+=writtenright;
      qendp+=writtenright;
      poslistp+=writtenright;
      distlistp+=writtenright;
      nbufrest-=writtenright;
    }

    Nqlist=index;


    //!
    //! edge tests between new nodes neglected
    //!

    //     ---

    //!
    //! calculate which edges exist
    //!
    //...... --> call indicator function on GPU
    //space->indicator2(qnew,num,qlist,resbuf,posqlist,numqlist,offset);
    space->indicator2(qstart,qend,resbuf,Nqlist,offset);


#ifndef NO_IO
    printarr(qnew,ndof*cnt);
    printvar(cnt);
    printvar(num);
    printarr(qstart,ndof*nbuf);
    printarr(qend,ndof*nbuf);
    printvar(Nqlist);
    printarr(posqlist,cnt);
    printarr(numqlist,cnt);
    printarr(numqlistleft,cnt);
    printvar(offset);
    printarr(distlist,nbuf);
    printarr(poslist,nbuf);
    printarr(resbuf,nbuf);
#endif


    //!
    //! insert nodes and edges
    //!

    //calculate connectivities

    //int *leftconn=leftconnall+dsp;
    //int *rightconn=rightconnall+dsp;
    for(int j=0;j<cnt;++j){

      //! left:  min,...,maxleft-1
      //! right: maxleft,...,max-1
      int min=posqlist[j];
      int maxleft=min+numqlistleft[j];
      int max=min+numqlist[j];

      leftconn[j]=0;
      for(int i=min;i<maxleft;++i){
        if(resbuf[i]==0){
          leftconn[j]=1;
          break;
        }
      }
      rightconn[j]=0;
      for(int i=maxleft;i<max;++i){
        if(resbuf[i]==0){
          rightconn[j]=1;
          break;
        }
      }

    }


    // --> finish send, recieve nodes
    MPI_Status qnewstatus;
    MPI_Wait(&qnewrequest,&qnewstatus);

    // --> send rightconn, leftconn here
    MPI_Allgather(leftconn,cnt,MPI_INT,leftconnall,cnt,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(rightconn,cnt,MPI_INT,rightconnall,cnt,MPI_INT,MPI_COMM_WORLD);

    // insert all nodes of previous processes
    for(int j=0;j<dsp;++j){
      if(leftconnall[j])positionsall[j]=insert(&qnewall[ndof*j],1,graphl);
      if(rightconnall[j])positionsall[j]=insert(&qnewall[ndof*j],1,graphr);
    }

    int connected=mpisize;

    for(int j=0;j<cnt;++j){

      int min=posqlist[j];
      int maxleft=min+numqlistleft[j];
      int max=min+numqlist[j];

      int positionl;
      if(leftconn[j]){
        //!
        //! connection to left graph exists -> insert in left graph
        //!
        positionl=insert(&qnew[ndof*j],1,graphl);
        //std::vector<int> *v=&(graphl.edgelists[positionl]);
        //std::vector<float> *w=&(graphl.edgeweights[positionl]);
        surrnumnew[j]=maxleft-min;
        graphl.surrnum[positionl]+=surrnumnew[j];
        for(int i=min;i<maxleft;++i){
          int goalpos=poslist[i];
          ++(graphl.surrnum[goalpos]);
          if(resbuf[i]==0){
            graphl.from.push_back(positionl);
            graphl.to.push_back(goalpos);
            float dist=distlist[i];
            graphl.weights.push_back(dist);
            //v->push_back(goalpos);
            //w->push_back(dist);
            //graphl.edgelists[goalpos].push_back(positionl);
            //graphl.edgeweights[goalpos].push_back(dist);
          }
        }
      }

      int positionr;
      if(rightconn[j]){
         //!
         //! connection to right graph exists -> insert in right graph
         //!
         positionr=insert(&qnew[ndof*j],1,graphr);
         //std::vector<int> *v=&(graphr.edgelists[positionr]);
         //std::vector<float> *w=&(graphr.edgeweights[positionr]);
         surrnumnew[j]=max-maxleft;
         graphr.surrnum[positionr]+=surrnumnew[j];
         for(int i=maxleft;i<max;++i){
           int goalpos=poslist[i];
           ++(graphr.surrnum[goalpos]);
           if(resbuf[i]==0){
             graphr.from.push_back(positionr);
             graphr.to.push_back(goalpos);
             float dist=distlist[i];
             graphr.weights.push_back(dist);
             //v->push_back(goalpos);
             //w->push_back(dist);
             //graphr.edgelists[goalpos].push_back(positionr);
             //graphr.edgeweights[goalpos].push_back(dist);
           }
         }
       }

       if(leftconn[j] && rightconn[j] && connected==mpisize){
         //!
         //!  Connection found! abort
         //!
         for(int l=0;l<ndof;++l)connection.q[l]=qnew[ndof*j+l];

         connection.index_left=positionl;
         connection.index_right=positionr;


         printvar(connection.index_left);
         printvar(connection.index_right);


         connected=mpirank;
       }

    }//for

    MPI_Request surrnumrequest;
    MPI_Iallgather(surrnumnew, cnt, MPI_INT, surrnumnewall, cnt, MPI_INT, MPI_COMM_WORLD,&surrnumrequest);


    //insert all nodes of following processes
    for(int j=dsp+cnt;j<num;++j){
      if(leftconnall[j])positionsall[j]=insert(&qnewall[ndof*j],1,graphl);
      if(rightconnall[j])positionsall[j]=insert(&qnewall[ndof*j],1,graphr);
    }

    MPI_Status statussurr;
    MPI_Wait(&surrnumrequest,&statussurr);

    //update surrnums for other processes
    for(int j=0;j<dsp;++j){
      if(leftconnall[j])(graphl.surrnum[positionsall[j]])+=surrnumnewall[j];
      if(rightconnall[j])(graphr.surrnum[positionsall[j]])+=surrnumnewall[j];
    }
    for(int j=dsp+cnt;j<num;++j){
      if(leftconnall[j])(graphl.surrnum[positionsall[j]])+=surrnumnewall[j];
      if(rightconnall[j])(graphr.surrnum[positionsall[j]])+=surrnumnewall[j];
    }


    //!
    //! connection checking
    //!

    //check if any process found a connection
    int connectedany=mpisize;
    //MPI_Allreduce(snd.data(), rcv.data(), snd.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&connected,&connectedany,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);


    if(connectedany<mpisize){
      tick(tconnect);

      printvar(connectedany);
      msg("gathering edges...");

      int root=0;

      build_edges(graphl,mpirank,mpisize,root);
      build_edges(graphr,mpirank,mpisize,root);

      if(connectedany==mpirank && mpirank!=root){
        msg("Sending");

        MPI_Send(&connection.index_left,1,MPI_INT,root,1000,MPI_COMM_WORLD);
        MPI_Send(&connection.index_right,1,MPI_INT,root,1001,MPI_COMM_WORLD);
        MPI_Send(&connection.q[0],ndof,MPI_FLOAT,root,1002,MPI_COMM_WORLD);
        msg("Sendt");
      }
      else if(mpirank==root && connectedany !=mpirank){
        msg("recieving");
        MPI_Status status[3];
        MPI_Recv(&connection.index_left,1,MPI_INT,connectedany,1000,MPI_COMM_WORLD,&status[0]);
        MPI_Recv(&connection.index_right,1,MPI_INT,connectedany,1001,MPI_COMM_WORLD,&status[1]);
        MPI_Recv(&connection.q[0],ndof,MPI_FLOAT,connectedany,1002,MPI_COMM_WORLD,&status[2]);
        msg("recieved");
      }

      printvar(graphl.newblockpos);
      printvar(graphr.newblockpos);
      printvar(blocksize);
      printvar(graphl.from.size());

      if(mpirank==root){
        int res0=do_dijkstra(graphl,dijkstral,i0l,connection.index_left);
        if(res0==0){msg("ERROR: no path found by dijkstra in graphl");}
        int res1=do_dijkstra(graphr,dijkstrar,connection.index_right,i0r);
        if(res1==0){msg("ERROR: no path found by dijkstra in graphr");}
      }

      tock(tconnect);
      return 1;
    }

    delete poslist;
    delete distlist;
    return 0;

  }

  template<int ndof>
  void PRMSolver<ndof>::build_edges(graph &g, int mpirank, int mpisize, int root){
      int numedges=g.from.size();
      int numedgesall;

      MPI_Reduce(&numedges, &numedgesall, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

      std::vector<int> numedgesvec(mpisize,0);
      MPI_Gather(&numedges, 1, MPI_FLOAT, numedgesvec.data(), 1, MPI_FLOAT, root, MPI_COMM_WORLD);

      if(mpirank==root){
        std::vector<int> fromall(numedgesall);
        std::vector<int> toall(numedgesall);
        std::vector<float> weightsall(numedgesall);
        std::vector<int> dsps(mpisize);
        dsps[0]=0;
        for(int i=0;i<mpisize-1;++i) dsps[i+1]=dsps[i]+numedgesvec[i];
        MPI_Gatherv(g.from.data(),numedges,MPI_INT,fromall.data(),numedgesvec.data(),dsps.data(),MPI_INT,root,MPI_COMM_WORLD);
        MPI_Gatherv(g.to.data(),numedges,MPI_INT,toall.data(),numedgesvec.data(),dsps.data(),MPI_INT,root,MPI_COMM_WORLD);
        MPI_Gatherv(g.weights.data(),numedges,MPI_INT,weightsall.data(),numedgesvec.data(),dsps.data(),MPI_INT,root,MPI_COMM_WORLD);

        for(int i=0;i<numedgesall;++i){
          int f=fromall[i];
          int t=toall[i];
          float w=weightsall[i];
          g.edgelists[f].push_back(t);
          g.edgeweights[f].push_back(w);
          g.edgelists[t].push_back(f);
          g.edgeweights[t].push_back(w);
        }
      }else{
        MPI_Gatherv(g.from.data(),numedges,MPI_INT,0x0,0x0,0x0,MPI_INT,root,MPI_COMM_WORLD);
        MPI_Gatherv(g.to.data(),numedges,MPI_INT,0x0,0x0,0x0,MPI_INT,root,MPI_COMM_WORLD);
        MPI_Gatherv(g.weights.data(),numedges,MPI_INT,0x0,0x0,0x0,MPI_INT,root,MPI_COMM_WORLD);
      }
  }






























  //! *******************
  //! *                 *
  //! *    process 2    *
  //! *                 *
  //! *******************



  //seed only relevant for root
  template<int ndof>
  int PRMSolver<ndof>::process_mpi2(const int num, const int nbuf, const int maxsteps, int seed){
    //int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm );
    MPI_Bcast( &seed, 1, MPI_INT, 0, MPI_COMM_WORLD );
    srand(seed);

    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    assert(num%size==0);

    float *qnew=new float[num*ndof];
    float *qstartlist=new float[ndof*nbuf];
    float *qendlist=new float[ndof*nbuf];
    int *resbuf=new int[nbuf];
    int *resbufloc=new int[nbuf];
    int offset=nbuf;

    int *leftconn=new int[num];
    int *rightconn=new int[num];

    for(int i=0;i<maxsteps;++i){
      int flag=processing_step2(rank,size,
                               qnew, num,
                               leftconn, rightconn,
                               qstartlist,qendlist,resbuf, resbufloc,nbuf,offset);
      if(flag==1){
        msg("connection found");
        printvar(i);
        break;
      }else if(i%50==0){
        printvar(i);
      }
    }

    delete qnew;
    delete qstartlist;
    delete qendlist;
    delete resbuf;
    delete resbufloc;
    delete leftconn;
    delete rightconn;

    return 0;
  }


  template<int ndof>
  int PRMSolver<ndof>::processing_step2(const int mpirank, const int mpisize,
                             float* qnew, const int num,
                             int *leftconn, int *rightconn,
                             float* qstart, float* qend, int* resbuf, int *resbufloc, const int nbuf, const int offset)
  {
      //!
      //! create nodes qnew randomly
      //!

      //!from left graph
      for(int j=0;j<num/2;++j){
        bool dismiss;
        do{
          //!choose random block
          int k;
          block *b;
          do{
              k=rand()%graphl.blocknum;
              b=&(graphl.blocks[k]);
          }while(b->num==0);
          //!chose random vertex
          int m=(b->pos+rand()%b->num);
          int l;
          int x=1+graphl.surrnum[m];
          int prob=RAND_MAX/(x*x*x);
          if(rand()>prob){
            dismiss=true;
          }else{
            l=ndof*m;
            for(int i=0;i<ndof;++i){
              qnew[ndof*j+i]=graphl.qstorage[l+i]-D+2*D*((float)rand()/RAND_MAX);
            }

            dismiss=space->indicator(&qnew[ndof*j])!=0;
          }
  #ifndef NO_IO
          printvar(l);
          printarr(&qnew[ndof*j],ndof);
          printvar(dismiss);
  #endif
        }while(dismiss);
      }
      //!from right graph
      for(int j=num/2;j<num;++j){
        bool dismiss;
        do{
          //!choose random block
          int k;
          block *b;
          do{
              k=rand()%graphr.blocknum;
              b=&(graphr.blocks[k]);
          }while(b->num==0);
          //!chose random vertex
          int m=(b->pos+rand()%b->num);
          int l;
          int x=1+graphr.surrnum[m];
          int prob=RAND_MAX/(x*x*x);
          if(rand()>prob){
            dismiss=true;
          }else{
            l=ndof*m;
            for(int i=0;i<ndof;++i){
              qnew[ndof*j+i]=graphr.qstorage[l+i]-D+2*D*((float)rand()/RAND_MAX);
            }
            dismiss=space->indicator(&qnew[ndof*j])!=0;
          }
  #ifndef NO_IO
          printvar(l);
          printarr(&qnew[ndof*j],ndof);
          printvar(dismiss);
  #endif
        }while(dismiss);
      }





      //!
      //! make list of potential neighbours for all new nodes
      //!


      //!posqlist[i]: position in qlist buffer, at which neighbours of qnew[i] start
      int posqlist[num];
      int numqlist[num];
      int numqlistleft[num];
      int Nqlist;             //sum(numqlist)

      int *poslist=new int[nbuf];
      float *distlist=new float[nbuf];
      int index=0;
      float *qstartp=qstart;
      float *qendp=qend;
      int *poslistp=poslist;
      float *distlistp=distlist;
      int nbufrest=nbuf;

      for(int j=0;j<num;++j){
        posqlist[j]=index;

        int writtenleft=get_near_vertices(&qnew[ndof*j],1,qendp,poslistp,distlistp,nbufrest,offset,graphl);

        for(int i=0;i<ndof;++i)
        for(int k=0;k<writtenleft;++k){
          qstartp[k+offset*i]=qnew[ndof*j+i];
        }
        numqlistleft[j]=writtenleft;
        index+=writtenleft;

        if(writtenleft>=nbufrest){
          for(int l=j+1;l<num;++l){
            posqlist[l]=nbuf;
            numqlist[l]=0;
            numqlistleft[l]=0;
            msg("buffer full");
          }
          numqlist[j]=writtenleft;
          break;
        }

        qstartp+=writtenleft;
        qendp+=writtenleft;
        poslistp+=writtenleft;
        distlistp+=writtenleft;
        nbufrest-=writtenleft;

        int writtenright=get_near_vertices(&qnew[ndof*j],1,qendp,poslistp,distlistp,nbufrest,offset,graphr);

        for(int i=0;i<ndof;++i)
        for(int k=0;k<writtenright;++k){
          qstartp[k+offset*i]=qnew[ndof*j+i];
        }
        numqlist[j]=writtenleft+writtenright;
        index+=writtenright;

        if(writtenright>=nbufrest){
          for(int l=j+1;l<num;++l){
            posqlist[l]=nbuf;
            numqlist[l]=0;
            numqlistleft[l]=0;
            msg("buffer full");
          }
          break;
        }

        qstartp+=writtenright;
        qendp+=writtenright;
        poslistp+=writtenright;
        distlistp+=writtenright;
        nbufrest-=writtenright;
      }

      Nqlist=index;


      //!
      //! edge tests between new nodes neglected
      //!

      //     ---

      //!
      //! calculate which edges exist
      //!
      //...... --> call indicator function on GPU
      //space->indicator2(qnew,num,qlist,resbuf,posqlist,numqlist,offset);
      //space->indicator2(qstart,qend,resbuf,Nqlist,offset);

      int r = Nqlist % mpisize, q = Nqlist / mpisize;
      int *counts=new int[mpisize];
      int *disps=new int[mpisize];
      for(int rank=0;rank<mpisize;++rank){
        int count = q;
        int disp = rank * q;
        count += rank < r ? 1 : 0;
        disp += rank < r ? rank : r;
        counts[rank]=count;
        disps[rank]=disp;
      }
      int disp=disps[mpirank];
      int count=counts[mpirank];
      space->indicator2(qstart+disp,qend+disp,resbufloc,count,offset);

      MPI_Allgatherv(resbufloc,count,MPI_INT,resbuf,counts,disps,MPI_INT,MPI_COMM_WORLD);

      delete counts;
      delete disps;

  #ifndef NO_IO
      printarr(qnew,ndof*num);
      printvar(num);
      printarr(qstart,ndof*nbuf);
      printarr(qend,ndof*nbuf);
      printvar(Nqlist);
      printarr(posqlist,num);
      printarr(numqlist,num);
      printarr(numqlistleft,num);
      printvar(offset);
      printarr(distlist,nbuf);
      printarr(poslist,nbuf);
      printarr(resbuf,nbuf);
  #endif



      //!
      //! insert nodes and edges
      //!


      for(int j=0;j<num;++j){

        //! left:  min,...,maxleft-1
        //! right: maxleft,...,max-1
        int min=posqlist[j];
        int maxleft=min+numqlistleft[j];
        int max=min+numqlist[j];

        leftconn[j]=0;
        for(int i=min;i<maxleft;++i){
          if(resbuf[i]==0){
            leftconn[j]=1;
            break;
          }
        }
        rightconn[j]=0;
        for(int i=maxleft;i<max;++i){
          if(resbuf[i]==0){
            rightconn[j]=1;
            break;
          }
        }

      }




      for(int j=0;j<num;++j){

        int min=posqlist[j];
        int maxleft=min+numqlistleft[j];
        int max=min+numqlist[j];

        int positionl;
        if(leftconn[j]){
          //!
          //! connection to left graph exists -> insert in left graph
          //!
          positionl=insert(&qnew[ndof*j],1,graphl);
          int surrnump=0;
          std::vector<int> *v=&(graphl.edgelists[positionl]);
          std::vector<float> *w=&(graphl.edgeweights[positionl]);
          for(int i=min;i<maxleft;++i){
            if(resbuf[i]==0){
              int goalpos=poslist[i];
              float dist=distlist[i];
              ++surrnump;
              v->push_back(goalpos);
              w->push_back(dist);
              ++(graphl.surrnum[goalpos]);
              graphl.edgelists[goalpos].push_back(positionl);
              graphl.edgeweights[goalpos].push_back(dist);
            }
          }
          graphl.surrnum[positionl]+=surrnump;
        }

        int positionr;
        if(rightconn[j]){
           //!
           //! connection to right graph exists -> insert in right graph
           //!
           positionr=insert(&qnew[ndof*j],1,graphr);
           int surrnump=0;
           std::vector<int> *v=&(graphr.edgelists[positionr]);
           std::vector<float> *w=&(graphr.edgeweights[positionr]);
           for(int i=maxleft;i<max;++i){
             if(resbuf[i]==0){
               int goalpos=poslist[i];
               float dist=distlist[i];
               ++surrnump;
               v->push_back(goalpos);
               w->push_back(dist);
               ++(graphr.surrnum[goalpos]);
               graphr.edgelists[goalpos].push_back(positionr);
               graphr.edgeweights[goalpos].push_back(dist);
             }
           }
           graphr.surrnum[positionr]+=surrnump;
         }

         if(leftconn[j] && rightconn[j]){
           //!
           //!  Connection found! abort
           //!
           for(int l=0;l<ndof;++l)connection.q[l]=qnew[ndof*j+l];

           connection.index_left=positionl;
           connection.index_right=positionr;

           //tick(tdijkstra);
           int res0=do_dijkstra(graphl,dijkstral,i0l,connection.index_left);
           int res1=do_dijkstra(graphr,dijkstrar,connection.index_right,i0r);
           if(res0==0){msg("ERROR: no path found by dijkstra in graphl");}
           if(res1==0){msg("ERROR: no path found by dijkstra in graphr");}
           //tock(tdijkstra);

           delete poslist;
           delete distlist;
           return 1;
         }

      }//for



      delete poslist;
      delete distlist;
      return 0;

    }
































  //! *******************
  //! *                 *
  //! *    process 3    *
  //! *                 *
  //! *******************



  //seed only relevant for root
  template<int ndof>
  int PRMSolver<ndof>::process_mpi3(const int num, const int nbuf, const int maxsteps, int seed){
    //int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm );
    MPI_Bcast( &seed, 1, MPI_INT, 0, MPI_COMM_WORLD );
    srand(seed);

    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    assert(num%size==0);
    int cnt_=num/size;
    int *dsp=new int[size];
    int *cnt=new int[size];
    for(int r=0; r<size;++r){
      dsp[r]=cnt_*r;
      cnt[r]=cnt_;
    }

    float *qnew=new float[ndof*num];
    float *qstartlist=new float[ndof*nbuf];
    float *qendlist=new float[ndof*nbuf];
    int *resbuf=new int[nbuf];
    int *resbufloc=new int[nbuf];
    int offset=nbuf;

    int *leftconn=new int[num];
    int *rightconn=new int[num];
    int *poslist=new int[nbuf];
    float *distlist=new float[nbuf];

    for(int i=0;i<maxsteps;++i){
      int flag=processing_step3(rank,size,
                               qnew, num, dsp, cnt,
                               leftconn, rightconn,
                               poslist, distlist,
                               qstartlist,qendlist,resbuf, resbufloc,nbuf,offset);
      if(flag==1){
        msg("connection found");
        printvar(i);
        break;
      }else if(i%50==0){
        printvar(i);
      }
    }

    delete qnew;
    delete qstartlist;
    delete qendlist;
    delete resbuf;
    delete resbufloc;
    delete leftconn;
    delete rightconn;
    delete poslist;
    delete distlist;
    delete dsp;
    delete cnt;

    return 0;
  }


  template<int ndof>
  int PRMSolver<ndof>::processing_step3(const int mpirank, const int mpisize,
                             float* qnew, const int num, const int *dsp, const int *cnt,
                             int *leftconn, int *rightconn,
                             int *poslist, float *distlist,
                             float* qstart, float* qend, int* resbuf, int *resbufloc, const int nbuf, const int offset)
  {
      //!
      //! create nodes qnew randomly
      //!

      get_random_nodes(graphl,0,num/2,qnew,D,space);
      get_random_nodes(graphr,num/2,num,qnew,D,space);


      //!
      //! make list of potential neighbours for all new nodes
      //!


      //!posqlist[i]: position in qlist buffer, at which neighbours of qnew[i] start
      int posqlist[num];
      int numqlist[num];
      int numqlistleft[num];
      int Nqlist;             //sum(numqlist)

      int index=0;
      float *qstartp=qstart;
      float *qendp=qend;
      int *poslistp=poslist;
      float *distlistp=distlist;
      int nbufrest=nbuf;

      for(int j=0;j<num;++j){
        posqlist[j]=index;

        int writtenleft=get_near_vertices(&qnew[ndof*j],1,qendp,poslistp,distlistp,nbufrest,offset,graphl);

        for(int i=0;i<ndof;++i)
        for(int k=0;k<writtenleft;++k){
          qstartp[k+offset*i]=qnew[ndof*j+i];
        }
        numqlistleft[j]=writtenleft;
        index+=writtenleft;

        if(writtenleft>=nbufrest){
          for(int l=j+1;l<num;++l){
            posqlist[l]=nbuf;
            numqlist[l]=0;
            numqlistleft[l]=0;
            msg("buffer full");
          }
          numqlist[j]=writtenleft;
          break;
        }

        qstartp+=writtenleft;
        qendp+=writtenleft;
        poslistp+=writtenleft;
        distlistp+=writtenleft;
        nbufrest-=writtenleft;

        int writtenright=get_near_vertices(&qnew[ndof*j],1,qendp,poslistp,distlistp,nbufrest,offset,graphr);

        for(int i=0;i<ndof;++i)
        for(int k=0;k<writtenright;++k){
          qstartp[k+offset*i]=qnew[ndof*j+i];
        }
        numqlist[j]=writtenleft+writtenright;
        index+=writtenright;

        if(writtenright>=nbufrest){
          for(int l=j+1;l<num;++l){
            posqlist[l]=nbuf;
            numqlist[l]=0;
            numqlistleft[l]=0;
            msg("buffer full");
          }
          break;
        }

        qstartp+=writtenright;
        qendp+=writtenright;
        poslistp+=writtenright;
        distlistp+=writtenright;
        nbufrest-=writtenright;
      }

      Nqlist=index;


      //!
      //! edge tests between new nodes neglected
      //!

      //     ---

      //!
      //! calculate which edges exist
      //!
      //...... --> call indicator function on GPU
      //space->indicator2(qnew,num,qlist,resbuf,posqlist,numqlist,offset);
      //space->indicator2(qstart,qend,resbuf,Nqlist,offset);

      int *counts=new int[mpisize];
      int *disps=new int[mpisize];
#if 0
      int r = Nqlist % mpisize, q = Nqlist / mpisize;
      for(int rank=0;rank<mpisize;++rank){
        int count = q;
        int disp = rank * q;
        count += rank < r ? 1 : 0;
        disp += rank < r ? rank : r;
        counts[rank]=count;
        disps[rank]=disp;
      }
#else
      for(int rank=0;rank<mpisize-1;++rank){
        int dsp_=dsp[rank];
        disps[rank]=posqlist[dsp_];
        counts[rank]=posqlist[dsp[rank+1]]-posqlist[dsp_];
      }
      disps[mpisize-1]=posqlist[dsp[mpisize-1]];
      counts[mpisize-1]=Nqlist-disps[mpisize-1];
#endif
      int disp=disps[mpirank];
      int count=counts[mpirank];

      printvar(disp);
      printvar(count);
      printarr(disps,mpisize);
      printarr(counts,mpisize);
      printvar(Nqlist);

      space->indicator2(qstart+disp,qend+disp,resbufloc,count,offset);

      MPI_Request resrequest;
      MPI_Iallgatherv(resbufloc,count,MPI_INT,resbuf,counts,disps,MPI_INT,MPI_COMM_WORLD,&resrequest);

      int dsp_=dsp[mpirank];
      int cnt_=cnt[mpirank];
      calc_conn(resbufloc-disp, posqlist, numqlistleft, numqlist, leftconn, rightconn, dsp_, dsp_+cnt_);


      MPI_Status resstatus;
      MPI_Wait(&resrequest,&resstatus);

      delete counts;
      delete disps;


  #ifndef NO_IO
      printarr(qnew,ndof*num);
      printvar(num);
      printarr(qstart,ndof*nbuf);
      printarr(qend,ndof*nbuf);
      printvar(Nqlist);
      printarr(posqlist,num);
      printarr(numqlist,num);
      printarr(numqlistleft,num);
      printvar(offset);
      printarr(distlist,nbuf);
      printarr(poslist,nbuf);
      printarr(resbuf,nbuf);
  #endif



      //!
      //! insert nodes and edges
      //!

      //calc_conn(resbufloc-disp, posqlist, numqlistleft, numqlist, leftconn, rightconn, dsp_, dsp_+cnt_);
      calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, 0, dsp_);
      calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, dsp_+cnt_, num);

      //calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, 0, num);



      for(int j=0;j<num;++j){

        int min=posqlist[j];
        int maxleft=min+numqlistleft[j];
        int max=min+numqlist[j];

        int positionl;
        if(leftconn[j]){
          //!
          //! connection to left graph exists -> insert in left graph
          //!
          positionl=insert(&qnew[ndof*j],1,graphl);
          int surrnump=0;
          std::vector<int> *v=&(graphl.edgelists[positionl]);
          std::vector<float> *w=&(graphl.edgeweights[positionl]);
          for(int i=min;i<maxleft;++i){
            if(resbuf[i]==0){
              int goalpos=poslist[i];
              float dist=distlist[i];
              ++surrnump;
              v->push_back(goalpos);
              w->push_back(dist);
              ++(graphl.surrnum[goalpos]);
              graphl.edgelists[goalpos].push_back(positionl);
              graphl.edgeweights[goalpos].push_back(dist);
            }
          }
          graphl.surrnum[positionl]+=surrnump;
        }

        int positionr;
        if(rightconn[j]){
           //!
           //! connection to right graph exists -> insert in right graph
           //!
           positionr=insert(&qnew[ndof*j],1,graphr);
           int surrnump=0;
           std::vector<int> *v=&(graphr.edgelists[positionr]);
           std::vector<float> *w=&(graphr.edgeweights[positionr]);
           for(int i=maxleft;i<max;++i){
             if(resbuf[i]==0){
               int goalpos=poslist[i];
               float dist=distlist[i];
               ++surrnump;
               v->push_back(goalpos);
               w->push_back(dist);
               ++(graphr.surrnum[goalpos]);
               graphr.edgelists[goalpos].push_back(positionr);
               graphr.edgeweights[goalpos].push_back(dist);
             }
           }
           graphr.surrnum[positionr]+=surrnump;
         }

         if(leftconn[j] && rightconn[j]){
           //!
           //!  Connection found! abort
           //!
           for(int l=0;l<ndof;++l)connection.q[l]=qnew[ndof*j+l];

           connection.index_left=positionl;
           connection.index_right=positionr;

           //tick(tdijkstra);
           int res0=do_dijkstra(graphl,dijkstral,i0l,connection.index_left);
           int res1=do_dijkstra(graphr,dijkstrar,connection.index_right,i0r);
           if(res0==0){msg("ERROR: no path found by dijkstra in graphl");}
           if(res1==0){msg("ERROR: no path found by dijkstra in graphr");}
           //tock(tdijkstra);

           return 1;
         }

      }//for

      return 0;

    }


































  //! *******************
  //! *                 *
  //! *    process 4    *
  //! *                 *
  //! *******************



  //seed only relevant for root
  template<int ndof>
  int PRMSolver<ndof>::process_mpi4(const int num, const int nbuf, const int maxsteps, int seed){
    //int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm );
    MPI_Bcast( &seed, 1, MPI_INT, 0, MPI_COMM_WORLD );
    srand(seed);

    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    assert(num%size==0);
    int cnt_=num/size;
    int *dsp=new int[size];
    int *cnt=new int[size];
    for(int r=0; r<size;++r){
      dsp[r]=cnt_*r;
      cnt[r]=cnt_;
    }

    float *qnew=new float[ndof*num];
    float *qstartlist=new float[ndof*nbuf];
    float *qendlist=new float[ndof*nbuf];
    int *resbuf=new int[nbuf];
    int *resbufloc=new int[nbuf];
    int offset=nbuf;

    int *leftconn=new int[num];
    int *rightconn=new int[num];
    int *poslist=new int[nbuf];
    float *distlist=new float[nbuf];

    for(int i=0;i<maxsteps;++i){
      int flag=processing_step4(rank,size,
                               qnew, num, dsp, cnt,
                               leftconn, rightconn,
                               poslist, distlist,
                               qstartlist,qendlist,resbuf, resbufloc,nbuf,offset);
      if(flag==1){
        msg("connection found");
        printvar(i);
        break;
      }else if(i%50==0){
        printvar(i);
      }
    }

    delete qnew;
    delete qstartlist;
    delete qendlist;
    delete resbuf;
    delete resbufloc;
    delete leftconn;
    delete rightconn;
    delete poslist;
    delete distlist;
    delete dsp;
    delete cnt;

    return 0;
  }



  template<int ndof>
  int PRMSolver<ndof>::processing_step4(const int mpirank, const int mpisize,
                             float* qnew, const int num, const int *dsp, const int *cnt,
                             int *leftconn, int *rightconn,
                             int *poslist, float *distlist,
                             float* qstart, float* qend, int* resbuf, int *resbufloc, const int nbuf, const int offset)
  {

      int dsp_=dsp[mpirank];
      int cnt_=cnt[mpirank];

      //!
      //! create nodes qnew randomly
      //!

      get_random_nodes(graphl,0,num/2,qnew,D,space);
      get_random_nodes(graphr,num/2,num,qnew,D,space);


      //!
      //! make list of potential neighbours for all new nodes
      //!


      //!posqlist[i]: position in qlist buffer, at which neighbours of qnew[i] start
      int posqlist[num];
      int numqlist[num];
      int numqlistleft[num];
      int Nqlist;             //sum(numqlist)

      int index=0;
      float *qstartp=qstart;
      float *qendp=qend;
      int *poslistp=poslist;
      float *distlistp=distlist;
      int nbufrest=nbuf;


      find_neighbours(  dsp_, dsp_+cnt_,
                        qnew, qstartp, qendp,
                        poslistp, distlistp, nbufrest, index,
                        &posqlist[0], &numqlistleft[0], &numqlist[0],
                        nbuf, offset
                      );



      int disp=posqlist[dsp[mpirank]]; //==0
      int count=index-disp;

      int configrequest;
      space->indicator2_async(qstart+disp,qend+disp,resbufloc,count,offset,configrequest);

      int *counts=new int[mpisize];
      int *disps=new int[mpisize];

      MPI_Allgather(&count,1,MPI_INT,counts,1,MPI_INT,MPI_COMM_WORLD);


      disps[mpirank]=0;
      int disps_=count;
      for(int rank=0;rank<mpisize;++rank){
        if(rank!=mpirank){
          disps[rank]=disps_;
          disps_+=counts[rank];
        }
      }


      find_neighbours(  0, dsp_,
                        qnew, qstartp, qendp,
                        poslistp, distlistp, nbufrest, index,
                        &posqlist[0], &numqlistleft[0], &numqlist[0],
                        nbuf, offset
                      );
      find_neighbours(  dsp_+cnt_, num,
                        qnew, qstartp, qendp,
                        poslistp, distlistp, nbufrest, index,
                        &posqlist[0], &numqlistleft[0], &numqlist[0],
                        nbuf, offset
                      );

      Nqlist=index;


      //!
      //! calculate which edges exist
      //!

      space->indicator2_async_wait(configrequest);

      MPI_Request resrequest;
      MPI_Iallgatherv(resbufloc,count,MPI_INT,resbuf,counts,disps,MPI_INT,MPI_COMM_WORLD,&resrequest);

      calc_conn(resbufloc-disp, posqlist, numqlistleft, numqlist, leftconn, rightconn, dsp_, dsp_+cnt_);

      MPI_Status resstatus;
      MPI_Wait(&resrequest,&resstatus);

      delete counts;
      delete disps;



  #ifndef NO_IO
      printarr(dsp,mpisize);
      printarr(cnt,mpisize);
      printvar(disp);
      printvar(count);
      printarr(disps,mpisize);
      printarr(counts,mpisize);
      printvar(Nqlist);
  #endif



      //!
      //! insert nodes and edges
      //!

      calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, 0, dsp_);
      calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, dsp_+cnt_, num);
      //calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, 0, num);



      for(int j=0;j<num;++j){

        int min=posqlist[j];
        int maxleft=min+numqlistleft[j];
        int max=min+numqlist[j];

        int positionl;
        if(leftconn[j]){
          //!
          //! connection to left graph exists -> insert in left graph
          //!
          positionl=insert(&qnew[ndof*j],1,graphl);
          int surrnump=0;
          std::vector<int> *v=&(graphl.edgelists[positionl]);
          std::vector<float> *w=&(graphl.edgeweights[positionl]);
          for(int i=min;i<maxleft;++i){
            if(resbuf[i]==0){
              int goalpos=poslist[i];
              float dist=distlist[i];
              ++surrnump;
              v->push_back(goalpos);
              w->push_back(dist);
              ++(graphl.surrnum[goalpos]);
              graphl.edgelists[goalpos].push_back(positionl);
              graphl.edgeweights[goalpos].push_back(dist);
            }
          }
          graphl.surrnum[positionl]+=surrnump;
        }

        int positionr;
        if(rightconn[j]){
           //!
           //! connection to right graph exists -> insert in right graph
           //!
           positionr=insert(&qnew[ndof*j],1,graphr);
           int surrnump=0;
           std::vector<int> *v=&(graphr.edgelists[positionr]);
           std::vector<float> *w=&(graphr.edgeweights[positionr]);
           for(int i=maxleft;i<max;++i){
             if(resbuf[i]==0){
               int goalpos=poslist[i];
               float dist=distlist[i];
               ++surrnump;
               v->push_back(goalpos);
               w->push_back(dist);
               ++(graphr.surrnum[goalpos]);
               graphr.edgelists[goalpos].push_back(positionr);
               graphr.edgeweights[goalpos].push_back(dist);
             }
           }
           graphr.surrnum[positionr]+=surrnump;
         }

         if(leftconn[j] && rightconn[j]){
           //!
           //!  Connection found! abort
           //!
           for(int l=0;l<ndof;++l)connection.q[l]=qnew[ndof*j+l];

           connection.index_left=positionl;
           connection.index_right=positionr;

           //tick(tdijkstra);
           int res0=do_dijkstra(graphl,dijkstral,i0l,connection.index_left);
           int res1=do_dijkstra(graphr,dijkstrar,connection.index_right,i0r);
           if(res0==0){msg("ERROR: no path found by dijkstra in graphl");}
           if(res1==0){msg("ERROR: no path found by dijkstra in graphr");}
           //tock(tdijkstra);

           return 1;
         }

      }//for

      return 0;

    }























  //! *******************
  //! *                 *
  //! *    process 5    *
  //! *                 *
  //! *******************



  //seed only relevant for root
  template<int ndof>
  int PRMSolver<ndof>::process_mpi5(const int numall, const int nbuf, const int maxsteps, int seed){
    //int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm );
    MPI_Bcast( &seed, 1, MPI_INT, 0, MPI_COMM_WORLD );
    srand(seed);

    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    MPI_Barrier(MPI_COMM_WORLD);

    assert(numall%(size*2)==0);

    int num=numall/2;

    int cnt_=num/size;
    int *dsp=new int[size];
    int *cnt=new int[size];
    for(int r=0; r<size;++r){
      dsp[r]=cnt_*r;
      cnt[r]=cnt_;
    }

    int offset=nbuf;

    float *qnew=new float[ndof*num];
    int *leftconn=new int[num];
    int *rightconn=new int[num];
    float *qstartlist=new float[ndof*nbuf];
    float *qendlist=new float[ndof*nbuf];
    int *resbuf=new int[nbuf];
    int *resbufloc=new int[nbuf];
    int *poslist=new int[nbuf];
    float *distlist=new float[nbuf];


    float *qnew2=new float[ndof*num];
    int *leftconn2=new int[num];
    int *rightconn2=new int[num];
    float *qstartlist2=new float[ndof*nbuf];
    float *qendlist2=new float[ndof*nbuf];
    int *resbuf2=new int[nbuf];
    int *resbufloc2=new int[nbuf];
    int *poslist2=new int[nbuf];
    float *distlist2=new float[nbuf];


    worker processor1(rank,size,
                         qnew, num, dsp, cnt,
                         leftconn, rightconn,
                         poslist, distlist,
                         qstartlist,qendlist,resbuf, resbufloc,nbuf,offset,
                         this);

    worker processor2(rank,size,
                         qnew2, num, dsp, cnt,
                         leftconn2, rightconn2,
                         poslist2, distlist2,
                         qstartlist2,qendlist2,resbuf2, resbufloc2,nbuf,offset,
                         this);



    processor1.processing_step_part1();
    processor2.processing_step_part1();

    //tick(tloop);
    for(int i=0;i<maxsteps;++i){
      //tick(evaluating1);
      processor1.processing_step_part2();
      int flag1=processor1.processing_step_part3();
      //tock(evaluating1);
      if(flag1==1){
        msg("connection found 1");
        printvar(i);
        break;
      }
      //tick(setting1);
      processor1.processing_step_part1();
      //tock(setting1);

      //tick(evaluating2);
      processor2.processing_step_part2();
      int flag2=processor2.processing_step_part3();
      //tock(evaluating2);
      if(flag2==1){
        msg("connection found 2");
        printvar(i);
        break;
      }
      //tick(setting2);
      processor2.processing_step_part1();
      //tock(setting2);

      if(i%10==0){
        printvar(i);
      }
    }
    //tock(tloop);

    delete qnew;
    delete qstartlist;
    delete qendlist;
    delete resbuf;
    delete resbufloc;
    delete leftconn;
    delete rightconn;
    delete poslist;
    delete distlist;
    delete dsp;
    delete cnt;

    delete qnew2;
    delete qstartlist2;
    delete qendlist2;
    delete resbuf2;
    delete resbufloc2;
    delete leftconn2;
    delete rightconn2;
    delete poslist2;
    delete distlist2;


    //msg("process_mpi5 finished");
    return 0;
  }



  template<int ndof>
  PRMSolver<ndof>::worker::worker( const int mpirank_, const int mpisize_,
                                   float* qnew_, const int num_, const int *dsp_, const int *cnt_,
                                   int *leftconn_, int *rightconn_,
                                   int *poslist_, float *distlist_,
                                   float* qstart_, float* qend_, int* resbuf_, int *resbufloc_, const int nbuf_, const int offset_,
                                   PRMSolver *instance_):
        mpirank(mpirank_), mpisize(mpisize_), num(num_), nbuf(nbuf_), offset(offset_),
        in(instance_), D(in->D), space(in->space), graphl(in->graphl), graphr(in->graphr)

    {
        qnew=qnew_; dsp=dsp_; cnt=cnt_;
        leftconn=leftconn_; rightconn=rightconn_;
        poslist=poslist_;
        distlist=distlist_;
        qstart=qstart_; qend=qend_; resbuf=resbuf_; resbufloc=resbufloc_;


        this->dsp_=dsp[mpirank];
        this->cnt_=cnt[mpirank];
        posqlist=new int[num];
        numqlist=new int[num];
        numqlistleft=new int[num];
        counts=new int[mpisize];
        disps=new int[mpisize];
    }

    template<int ndof>
    PRMSolver<ndof>::worker::~worker(){
        delete posqlist;
        delete numqlist;
        delete numqlistleft;
        delete counts;
        delete disps;
    }




    template<int ndof>
    int PRMSolver<ndof>::worker::processing_step_part1()
    {


      //!
      //! create nodes qnew randomly
      //!



      get_random_nodes(graphl,0,num/2,qnew,D,space);
      get_random_nodes(graphr,num/2,num,qnew,D,space);


      //!
      //! make list of potential neighbours for all new nodes
      //!


      //!posqlist[i]: position in qlist buffer, at which neighbours of qnew[i] start

      index=0;
      qstartp=qstart;
      qendp=qend;
      poslistp=poslist;
      distlistp=distlist;
      nbufrest=nbuf;


      in->find_neighbours(  dsp_, dsp_+cnt_,
                        qnew, qstartp, qendp,
                        poslistp, distlistp, nbufrest, index,
                        &posqlist[0], &numqlistleft[0], &numqlist[0],
                        nbuf, offset
                      );



      disp=posqlist[dsp[mpirank]]; //==0
      count=index-disp;

      space->indicator2_async(qstart+disp,qend+disp,resbufloc,count,offset,configrequest);


      MPI_Allgather(&count,1,MPI_INT,counts,1,MPI_INT,MPI_COMM_WORLD);


      disps[mpirank]=0;
      int disps_=count;
      for(int rank=0;rank<mpisize;++rank){
        if(rank!=mpirank){
          disps[rank]=disps_;
          disps_+=counts[rank];
        }
      }


      in->find_neighbours(  0, dsp_,
                        qnew, qstartp, qendp,
                        poslistp, distlistp, nbufrest, index,
                        &posqlist[0], &numqlistleft[0], &numqlist[0],
                        nbuf, offset
                      );
      in->find_neighbours(  dsp_+cnt_, num,
                        qnew, qstartp, qendp,
                        poslistp, distlistp, nbufrest, index,
                        &posqlist[0], &numqlistleft[0], &numqlist[0],
                        nbuf, offset
                      );

      Nqlist=index;


      //!
      //! calculate which edges exist
      //!


      return 0;
  }



    template<int ndof>
    int PRMSolver<ndof>::worker::processing_step_part2(){

#ifndef NO_IO
    printarr(dsp,mpisize);
    printarr(cnt,mpisize);
    printvar(disp);
    printvar(count);
    printarr(disps,mpisize);
    printarr(counts,mpisize);
    printvar(Nqlist);
#endif

      //tick(waiting);
      space->indicator2_async_wait(configrequest);
      //tock(waiting);


      MPI_Request resrequest;
      MPI_Iallgatherv(resbufloc,count,MPI_INT,resbuf,counts,disps,MPI_INT,MPI_COMM_WORLD,&resrequest);

      calc_conn(resbufloc-disp, posqlist, numqlistleft, numqlist, leftconn, rightconn, dsp_, dsp_+cnt_);

      MPI_Status resstatus;
      MPI_Wait(&resrequest,&resstatus);



      //!
      //! insert nodes and edges
      //!

      calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, 0, dsp_);
      calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, dsp_+cnt_, num);
      //calc_conn(resbuf, posqlist, numqlistleft, numqlist, leftconn, rightconn, 0, num);

      return 0;
    }



    template<int ndof>
    int PRMSolver<ndof>::worker::processing_step_part3(){

      for(int j=0;j<num;++j){

        int min=posqlist[j];
        int maxleft=min+numqlistleft[j];
        int max=min+numqlist[j];

        int positionl;
        if(leftconn[j]){
          //!
          //! connection to left graph exists -> insert in left graph
          //!
          positionl=in->insert(&qnew[ndof*j],1,graphl);
          int surrnump=0;
          std::vector<int> *v=&(graphl.edgelists[positionl]);
          std::vector<float> *w=&(graphl.edgeweights[positionl]);
          for(int i=min;i<maxleft;++i){
            if(resbuf[i]==0){
              int goalpos=poslist[i];
              float dist=distlist[i];
              ++surrnump;
              v->push_back(goalpos);
              w->push_back(dist);
              ++(graphl.surrnum[goalpos]);
              graphl.edgelists[goalpos].push_back(positionl);
              graphl.edgeweights[goalpos].push_back(dist);
            }
          }
          graphl.surrnum[positionl]+=surrnump;
        }

        int positionr;
        if(rightconn[j]){
           //!
           //! connection to right graph exists -> insert in right graph
           //!
           positionr=in->insert(&qnew[ndof*j],1,graphr);
           int surrnump=0;
           std::vector<int> *v=&(graphr.edgelists[positionr]);
           std::vector<float> *w=&(graphr.edgeweights[positionr]);
           for(int i=maxleft;i<max;++i){
             if(resbuf[i]==0){
               int goalpos=poslist[i];
               float dist=distlist[i];
               ++surrnump;
               v->push_back(goalpos);
               w->push_back(dist);
               ++(graphr.surrnum[goalpos]);
               graphr.edgelists[goalpos].push_back(positionr);
               graphr.edgeweights[goalpos].push_back(dist);
             }
           }
           graphr.surrnum[positionr]+=surrnump;
         }

         if(leftconn[j] && rightconn[j]){
           //!
           //!  Connection found! abort
           //!

          connection_data &connection=in->connection;
          dijkstra_result &dijkstral=in->dijkstral;
          dijkstra_result &dijkstrar=in->dijkstrar;
          int &i0l=in->i0l, &i0r=in->i0r;

           for(int l=0;l<ndof;++l)connection.q[l]=qnew[ndof*j+l];

           connection.index_left=positionl;
           connection.index_right=positionr;

           //tick(tdijkstra);
           int res0=in->do_dijkstra(graphl,dijkstral,i0l,connection.index_left);
           int res1=in->do_dijkstra(graphr,dijkstrar,connection.index_right,i0r);
           if(res0==0){msg("ERROR: no path found by dijkstra in graphl");}
           if(res1==0){msg("ERROR: no path found by dijkstra in graphr");}
           //tock(tdijkstra);

           return 1;
         }

      }//for

      return 0;

    }






  //! **********************
  //! *                    *
  //! *    help methods    *
  //! *                    *
  //! **********************

  //!insert node q, data of surrnum and edges are not modified -> position returned for this
  template<int ndof>
  int PRMSolver<ndof>::insert(const float* q, const int offset, graph& g){
    int key=calc_key(q[0]);
    piterator it = g.map.find(key);
    block *b;
    if(it==g.map.end()){
      b=&(g.blocks[g.blocknum++]);
      g.map[key]=b;
      b->pos=g.newblockpos;
      g.newblockpos+=blocksize;
      if(g.newblockpos>N) return -1;
      b->num=0;
    }else{
      b=it->second;
      while(b->num>=blocksize){
        b=b->next;
      }
    }
    int position=b->pos+b->num++;
    int qposition=ndof*position;
    for(int i=0;i<ndof;++i){
      g.qstorage[qposition+i]=q[offset*i];
    }
    if(b->num>=blocksize){
      block *bnew=&(g.blocks[g.blocknum++]);
      b->next=bnew;
      bnew->pos=g.newblockpos;
      g.newblockpos+=blocksize;
      if(g.newblockpos>N) return -1;
      bnew->num=0;
    }
    return position;
  }


  //! get list of all vertices nearer than D
  //! qref: vertex to process
  //! qlist: buffer to store neighbour candidates, struct of arrays
  //! offset: i-th component of k-th vec in qlist is qlist[i*offset+k]
  //! nbuf: length(qlist)
  //! D: maximal distance
  template<int ndof>
  int PRMSolver<ndof>::get_near_vertices(const float* qref, float* qlist, const int& nbuf, const int& offset, const graph& g) const{
    const int keylower=calc_key(qref[0]-D);
    const_piterator begin=g.map.lower_bound(keylower);
    const int keyupper=calc_key(qref[0]+D);
    const_piterator end=g.map.upper_bound(keyupper);
    int index[ndof];
    for(int i=0;i<ndof;++i){
        index[i]=i*offset;
    }
    for(;!(begin==end);++begin){
      const block *b=begin->second;
      bool more=true;
      while(more){
        more=b->num>=blocksize;
        printvar(b->num);
        const int pos=ndof*b->pos;
        const int max=pos+ndof*b->num;
        for(int k=pos;k<max;k+=ndof){
          //!calc norm
          float normsq=0;
          for(int i=0;i<ndof;++i){
            float diff=g.qstorage[k+i]-qref[i];
            normsq+=diff*diff;
          }
          //!compare norm
          if(normsq<D2){
            //! store neighbour in buffer
            for(int i=0;i<ndof;++i){
              qlist[index[i]++]=g.qstorage[k+i];
            }
            //! terminate if buffer full
            if(index[0]==nbuf)return nbuf;
          }
        }//for
        b=b->next;
      }
    }//for
    //!number of written q's
    return index[0];
  }

  //!
  //! \brief get_near_vertices  get nodes nearer than D to base node
  //! \param qref               base node
  //! \param offsetref          offset for base node
  //! \param qlist              near nodes data
  //! \param posqlist           near nodes indices
  //! \param distlist           distances
  //! \param nbuf               maximal number of near nodes (arrays must have length nbuf resp. ndof*nbuf)
  //! \param offset             offset of buffers
  //! \param g                  the graph
  //! \return
  //!
  template<int ndof>
  int PRMSolver<ndof>::get_near_vertices(const float* qref, const int& offsetref, float* qlist, int* posqlist, float* distlist, const int& nbuf, const int& offset, const graph& g) const{
    const int keylower=calc_key(qref[0]-D);
    const_piterator begin=g.map.lower_bound(keylower);
    const int keyupper=calc_key(qref[0]+D);
    const_piterator end=g.map.upper_bound(keyupper);
    int index[ndof];
    for(int i=0;i<ndof;++i){
        index[i]=i*offset;
    }
    for(;!(begin==end);++begin){
      const block *b=begin->second;
      bool more=true;
      while(more){
        more=b->num>=blocksize;
        const int pos=b->pos;
        const int max=pos+b->num;
        for(int l=pos;l<max;++l){
          int k=ndof*l;
          //!calc norm
          float normsq=0;
          for(int i=0;i<ndof;++i){
            float diff=g.qstorage[k+i]-qref[offsetref*i];
            normsq+=diff*diff;
          }
          //!compare norm
          if(normsq<D2){
            //! store neighbour in buffer
            posqlist[index[0]]=l;
            distlist[index[0]]=sqrt(normsq);
            for(int i=0;i<ndof;++i){
              qlist[index[i]++]=g.qstorage[k+i];
            }
            //! terminate if buffer full
            if(index[0]==nbuf)return nbuf;
          }
        }//for
        b=b->next;
      }
    }//for
    //!number of written q's
    return index[0];
  }

  //!
  //! \brief get_random_nodes   get random feasible nodes from graph
  //! \param g                  the graph
  //! \param start              start index in qnew
  //! \param end                end (excl.) index
  //! \param qnew               node storage (array of structs)
  //! \return
  //!
  template<int ndof>
  int PRMSolver<ndof>::get_random_nodes(const graph &g, const int start, const int end, float *qnew, float D, Configspace<ndof> *space){
      for(int j=start;j<end;++j){
        bool dismiss;
        do{
          //!choose random block
          int k;
          const block *b;
          do{
              k=rand()%g.blocknum;
              b=&(g.blocks[k]);
          }while(b->num==0);
          //!chose random vertex
          int m=(b->pos+rand()%b->num);
          int l;
          int x=1+g.surrnum[m];
          int prob=RAND_MAX/(x*x*x);
          if(rand()>prob){
            dismiss=true;
          }else{
            l=ndof*m;
            for(int i=0;i<ndof;++i){
              qnew[ndof*j+i]=g.qstorage[l+i]-D+2*D*((float)rand()/RAND_MAX);
            }

            dismiss=space->indicator(&qnew[ndof*j])!=0;
          }
  #ifndef NO_IO
          printvar(l);
          printarr(&qnew[ndof*j],ndof);
          printvar(dismiss);
  #endif
        }while(dismiss);
      }
  }



  template<int ndof>
  void PRMSolver<ndof>::find_neighbours(int from, int to,
                              float *qnew, float *&qstartp, float *&qendp,
                              int *&poslistp, float *&distlistp, int &nbufrest, int &index,
                              int *posqlist, int *numqlistleft, int *numqlist,
                              int nbuf, int offset
                              ) const{

      for(int j=from;j<to;++j){
        posqlist[j]=index;

        int writtenleft=get_near_vertices(&qnew[ndof*j],1,qendp,poslistp,distlistp,nbufrest,offset,graphl);

        for(int i=0;i<ndof;++i)
        for(int k=0;k<writtenleft;++k){
          qstartp[k+offset*i]=qnew[ndof*j+i];
        }
        numqlistleft[j]=writtenleft;
        index+=writtenleft;

        assert(writtenleft<=nbufrest);

        qstartp+=writtenleft;
        qendp+=writtenleft;
        poslistp+=writtenleft;
        distlistp+=writtenleft;
        nbufrest-=writtenleft;

        int writtenright=get_near_vertices(&qnew[ndof*j],1,qendp,poslistp,distlistp,nbufrest,offset,graphr);

        for(int i=0;i<ndof;++i)
        for(int k=0;k<writtenright;++k){
          qstartp[k+offset*i]=qnew[ndof*j+i];
        }
        numqlist[j]=writtenleft+writtenright;
        index+=writtenright;

        assert(writtenright<=nbufrest);


        qstartp+=writtenright;
        qendp+=writtenright;
        poslistp+=writtenright;
        distlistp+=writtenright;
        nbufrest-=writtenright;
      }

  }


  //!
  //! \brief calc_conn      calculate, which nodes have edges and therefore can be inserted in graphl/graphr
  //! \param resbuf         0=edge, other value=no edge
  //! \param posqlist       displacement array
  //! \param numqlistleft   count array 1
  //! \param numqlist       count array 2
  //! \param leftconn       result graphl
  //! \param rightconn      result graphr
  //! \param from           start index to treat
  //! \param to             end (excl.) index
  //!
  template<int ndof>
  void PRMSolver<ndof>::calc_conn(const int *resbuf, const int *posqlist, const int *numqlistleft, const int *numqlist, int *leftconn, int *rightconn, const int from, const int to){
        for(int j=from;j<to;++j){

          //! left:  min,...,maxleft-1
          //! right: maxleft,...,max-1
          int min=posqlist[j];
          int maxleft=min+numqlistleft[j];
          int max=min+numqlist[j];

          leftconn[j]=0;
          for(int i=min;i<maxleft;++i){
            if(resbuf[i]==0){
              leftconn[j]=1;
              break;
            }
          }
          rightconn[j]=0;
          for(int i=maxleft;i<max;++i){
            if(resbuf[i]==0){
              rightconn[j]=1;
              break;
            }
          }
        }
    }








    //!
    //! \brief do_dijkstra  Dijkstra algorithm implementation
    //! \param g            the graph
    //! \param d            the result
    //! \param from         start index
    //! \param to           goal index
    //! \return
    //!
  template<int ndof>
  int PRMSolver<ndof>::do_dijkstra(graph& g, dijkstra_result& d, int from, int to){
    d.dist.assign(g.newblockpos,FLT_MAX);
    d.dist[from]=0.0;
    d.parent.resize(g.newblockpos,-1);
    d.parent[from]=from;
    std::set<std::pair<float,int>> queue;

    queue.insert(std::pair<float,int>(0.0,from));

    while(!queue.empty()){
      auto it=queue.begin();
      float dist=it->first;
      int pos=it->second;
      queue.erase(it);

      if(pos==to){
        //!
        //! fertig
        //!
        d.path.assign(1,to);
        int pos=to;
        while(pos!=from){
          pos=d.parent[pos];
          d.path.push_back(pos);
        }
        return 1;
      }

      for(int i=0;i<g.edgelists[pos].size();++i){
        int dest=g.edgelists[pos][i];
        float weight=g.edgeweights[pos][i];
        if(d.dist[dest]>dist+weight){
          queue.erase(std::pair<float,int>(d.dist[dest],dest));
          d.dist[dest]=dist+weight;
          queue.insert(std::pair<float,int>(d.dist[dest],dest));
          d.parent[dest]=pos;
        }
      }//for

    }//while
    return 0; //no path
  }









  //! ************************
  //! *                      *
  //! *    output methods    *
  //! *                      *
  //! ************************


  template<int ndof>
  void PRMSolver<ndof>::print(){
    msg("-------PRMSolver-------");
    msg("graphl:");
    for(piterator it=graphl.map.begin();it!=graphl.map.end();++it){
      block *b=it->second;
      bool more=false;
      do{
        msg("--");
        printvar(b->pos);
        printvar(b->num);
        for(int i=0;i<b->num;++i){

          float* qi=&(graphl.qstorage.data()[ndof*(b->pos+i)]);
          printarr(qi,ndof);
          int* goalsi=graphl.edgelists[b->pos+i].data();
          printarr(goalsi,graphl.edgelists[b->pos+i].size());
          float* weightsi=graphl.edgeweights[b->pos+i].data();
          printarr(weightsi,graphl.edgeweights[b->pos+i].size());
        }
        more=(b->num>=blocksize);
        b=b->next;
      }while(more);
    }
    msg("----------");
    msg("graphr:");
    for(piterator it=graphr.map.begin();it!=graphr.map.end();++it){
      block *b=it->second;
      bool more=false;
      do{
        msg("--");
        printvar(b->pos);
        printvar(b->num);
        for(int i=0;i<b->num;++i){

          float* qi=&(graphr.qstorage.data()[ndof*(b->pos+i)]);
          printarr(qi,ndof);
          int* goalsi=graphr.edgelists[b->pos+i].data();
          printarr(goalsi,graphr.edgelists[b->pos+i].size());
          float* weightsi=graphr.edgeweights[b->pos+i].data();
          printarr(weightsi,graphr.edgeweights[b->pos+i].size());
        }
        more=(b->num>=blocksize);
        b=b->next;
      }while(more);
    }

    msg("------------------------");
  }

  template<int ndof>
  void PRMSolver<ndof>::store_results(std::string path){
    int rank, root=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank!=root)return;

    msg("Storing results....");
    printvar(connection.index_left);
    printvar(connection.index_right);

    int ret0=system(("rm -rf "+path).c_str());
    if(ret0!=0)printvar(ret0);

    std::string pathl=path+"/graphl";
    //system("rm -rf "+pathl);
    int ret1=system(("mkdir -p "+pathl).c_str());
    if(ret1!=0)printvar(ret1);
    store_graph(pathl,graphl,dijkstral,i0l,connection.index_left);

    std::string pathr=path+"/graphr";
    //system("rm -rf "+pathr);
    int ret2=system(("mkdir "+pathr).c_str());
    if(ret2!=0)printvar(ret2);
    store_graph(pathr,graphr,dijkstrar,connection.index_right,i0r);

    write_file(path+"/connection.bin",(int*)&connection_found,1);

  }

  template<int ndof>
  void PRMSolver<ndof>::store_graph(std::string path, graph& g, dijkstra_result& d, int start, int end) const {
    write_file(path + "/qstorage.bin",g.qstorage.data(),ndof*g.newblockpos);
    write_file(path + "/surrnum.bin",g.surrnum.data(), g.newblockpos);

    std::vector<int> edgesfrom;
    std::vector<int> edgesto;
    std::vector<int> edgesweight;

    std::vector<int> blockpos;
    std::vector<int> blocknum;

    //!compressed arrays
    std::vector<int> edgesfromc;
    std::vector<int> edgestoc;
    std::vector<float> qstoragec;

    std::vector<int> posmapping(N,0);

    int index=0;
    for(int l=0;l<g.blocknum;++l){
      block *b=&(g.blocks[l]);
      blockpos.push_back(b->pos);
      blocknum.push_back(b->num);
      for(int l=b->pos;l<b->pos+b->num;++l){
        for(int i=0;i<ndof;++i){
          qstoragec.push_back(g.qstorage[ndof*l+i]);
        }
        posmapping[l]=index;
        for(int i=0;i<g.edgelists[l].size();++i){
          edgesfrom.push_back(l);
          edgesto.push_back(g.edgelists[l][i]);
          edgesweight.push_back(g.edgeweights[l][i]);
        }
        ++index;
      }
    }
    int numc=index;

    for(int l=0;l<g.blocknum;++l){
      block *b=&(g.blocks[l]);
      for(int l=b->pos;l<b->pos+b->num;++l){
        for(int i=0;i<g.edgelists[l].size();++i){
          edgesfromc.push_back(posmapping[l]);
          edgestoc.push_back(posmapping[g.edgelists[l][i]]);
        }
      }
    }

    write_file(path + "/edgesfrom.bin",edgesfrom.data(),edgesfrom.size());
    write_file(path + "/edgesto.bin",edgesto.data(),edgesto.size());
    write_file(path + "/edgesfromc.bin",edgesfromc.data(),edgesfromc.size());
    write_file(path + "/edgestoc.bin",edgestoc.data(),edgestoc.size());
    write_file(path + "/edgesweight.bin",edgesweight.data(),edgesweight.size());

    write_file(path + "/qstoragec.bin",qstoragec.data(),qstoragec.size());


    write_file(path + "/blockpos.bin",blockpos.data(),blockpos.size());
    write_file(path + "/blocknum.bin",blocknum.data(),blocknum.size());

    write_file(path + "/newblockpos.bin",&(g.newblockpos),1);
    write_file(path + "/blocknum.bin",&(g.blocknum),1);
    write_file(path + "/numc.bin",&numc,1);


    write_file(path + "/startc.bin",&(posmapping[start]),1);
    write_file(path + "/endc.bin",&(posmapping[end]),1);

    //!plot dijkstra
    if(connection_found){
      std::vector<int> pathc;
      for(int i=0;i<d.path.size();++i){
        pathc.push_back(posmapping[d.path[i]]);
      }

      write_file(path +"/path.bin",d.path.data(), d.path.size());
      write_file(path +"/pathc.bin",pathc.data(), pathc.size());
    }

  }




  template class PRMSolver<4>;
  template class PRMSolver<2>;




