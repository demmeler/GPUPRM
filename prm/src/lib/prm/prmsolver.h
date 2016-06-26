#ifndef PRMSOLVER_H
#define PRMSOLVER_H

#include <mpi.h>
#include <map>
#include <vector>
#include <set>
#include <string>

#define MPI_CODE

#include <config.h>

#include "configspace.hpp"





template<int ndof>
class PRMSolver{

    //! *******************
    //! *    subtypes     *
    //! *******************

  struct block;

  struct block{
    int pos; //! position of first vertex
    int num; //! current number of vertices stored
    int main_num;
    block* next; //!if num==blocksize -> pointer to next block
  };

  typedef std::map<int,block*> pmap;
  typedef typename pmap::iterator piterator;
  typedef typename pmap::const_iterator const_piterator;


  struct graph{
    pmap map; //map for sorting on high level
    std::vector<block> blocks;
    std::vector<int> main_keys;

    std::vector<float> qstorage;                 //length ndof*N
    std::vector<int> surrnum;                    //length N
    std::vector<std::vector<int>> edgelists;     //length N
    std::vector<std::vector<float>> edgeweights; //length N

    std::vector<int> from;
    std::vector<int> to;
    std::vector<int> weights;

    int newblockpos;  //position of next block in size-N-arrays
    int blocknum;     //number of used blocks

    std::string label;
  };

  struct dijkstra_result{
    std::vector<int> parent; //length newblockpos
    std::vector<float> dist; //length newblockpos
    std::vector<int> path; //computed path
  };

    //! *******************
    //! *      class      *
    //! *******************


public:

  PRMSolver(Configspace<ndof> *space_, float H_, float D_, int N=1024*1024, int blocksize=256);
  ~PRMSolver();
  int init(const float* qstart, const float* qend);

  //! ***************************
  //! *                         *
  //! *   processing routine    *
  //! *        versions         *
  //! *                         *
  //! ***************************



  //! *******************
  //! *                 *
  //! *    process 1    * (=version 1 in paper)
  //! *                 *
  //! *******************

  int process_mpi(int num, const int nbuf, const int maxsteps);

  //! get list of all vertices nearer than D
  //! qlist: buffer to store neighbour candidates, struct of arrays
  //! offset: i-th component of k-th vec in qlist is qlist[i*offset+k]
  //! nbuf: length(qlist)
  //! D: maximal distance
  //! return: 0 - no errors, 1 - connection found
  inline int processing_step(const int mpirank, const int mpisize,
                             float* qnew, float* qnewall, const int num, const int dsp, const int cnt,
                             int *surrnumnew, int *surrnumnewall, int *positionsall,
                             int* leftconn, int *leftconnall, int *rightconn, int *rightconnall,
                             float* qstart, float* qend, int* resbuf, const int nbuf, const int offset);

  void build_edges(graph &g, int mpirank, int mpisize, int root);


  //! *******************
  //! *                 *
  //! *    process 2    *
  //! *                 *
  //! *******************



  //seed only relevant for root
  int process_mpi2(const int num, const int nbuf, const int maxsteps, int seed);

  inline int processing_step2(const int mpirank, const int mpisize,
                             float* qnew, const int num,
                             int *leftconn, int *rightconn,
                             float* qstart, float* qend, int* resbuf, int *resbufloc, const int nbuf, const int offset);


  //! *******************
  //! *                 *
  //! *    process 3    *
  //! *                 *
  //! *******************



  //seed only relevant for root
  int process_mpi3(const int num, const int nbuf, const int maxsteps, int seed);

  inline int processing_step3(const int mpirank, const int mpisize,
                             float* qnew, const int num, const int *dsp, const int *cnt,
                             int *leftconn, int *rightconn,
                             int *poslist, float *distlist,
                             float* qstart, float* qend, int* resbuf, int *resbufloc, const int nbuf, const int offset);


  //! *******************
  //! *                 *
  //! *    process 4    * (=version 2 in paper)
  //! *                 *
  //! *******************

  //seed only relevant for root
  int process_mpi4(const int num, const int nbuf, const int maxsteps, int seed);



  inline int processing_step4(const int mpirank, const int mpisize,
                             float* qnew, const int num, const int *dsp, const int *cnt,
                             int *leftconn, int *rightconn,
                             int *poslist, float *distlist,
                             float* qstart, float* qend, int* resbuf, int *resbufloc, const int nbuf, const int offset);

  //! *******************
  //! *                 *
  //! *    process 5    * (=version 3 in paper)
  //! *                 *
  //! *******************

  //seed only relevant for root
  int process_mpi5(const int numall, const int nbuf, const int maxsteps, int seed, int workerversion=1);


  friend class worker;

  class worker{
  protected:
      const int mpirank;
      const int mpisize;
      float* qnew;
      const int num;
      const int *dsp;
      const int *cnt;
      int *leftconn;
      int *rightconn;
      int *poslist;
      float *distlist;
      float* qstart;
      float* qend;
      int* resbuf;
      int *resbufloc;
      const int nbuf;
      const int offset;

      int dsp_, cnt_;
      int *posqlist;
      int *numqlist;
      int *numqlistleft;
      int Nqlist;             //sum(numqlist)

      int index;
      float *qstartp;
      float *qendp;
      int *poslistp;
      float *distlistp;
      int nbufrest;

      int configrequest;
      int disp,count;
      int *counts, *disps;

      PRMSolver<ndof> *in;

      float &D;
      Configspace<ndof>* &space;
      graph &graphl;
      graph &graphr;

      int (worker::*step1ptr)();
      int (worker::*step2ptr)();
      int (worker::*step3ptr)();

  public:
    worker( const int mpirank_, const int mpisize_,
                                   float* qnew_, const int num_, const int *dsp_, const int *cnt_,
                                   int *leftconn_, int *rightconn_,
                                   int *poslist_, float *distlist_,
                                   float* qstart_, float* qend_, int* resbuf_, int *resbufloc_, const int nbuf_, const int offset_,
                                   PRMSolver *instance_, int workerversion=1);
    ~worker();

    int processing_step_part1();
    int processing_step_part2();
    int processing_step_part3();

  protected:
    //! workerversion 1
    int processing_step_part1a();
    int processing_step_part2a();
    int processing_step_part3a();

    //! workerversion 2
    int processing_step_part1b();
    int processing_step_part2b();

    //! workerversion 3
    int processing_step_part1c();
  };



  //! **********************
  //! *                    *
  //! *    help methods    *
  //! *                    *
  //! **********************

  //!
  //! \brief calc_key   mapping q -> key
  //! \param component  q[0]
  //! \return
  //!
  inline int calc_key(const float& component) const{
    return (int)(component*factor);
  }


  //!insert node q, data of surrnum and edges are not modified -> position returned for this
  inline int insert(const float* q, graph& g){
    return insert(q,1,g);
  }

  //!insert node q, data of surrnum and edges are not modified -> position returned for this
  int insert(const float* q, const int offset, graph& g);

  //! get list of all vertices nearer than D
  //! qref: vertex to process
  //! qlist: buffer to store neighbour candidates, struct of arrays
  //! offset: i-th component of k-th vec in qlist is qlist[i*offset+k]
  //! nbuf: length(qlist)
  //! D: maximal distance
  inline int get_near_vertices(const float* qref, float* qlist, const int& nbuf, const int& offset, const graph& g) const;

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
  inline int get_near_vertices(const float* qref, const int& offsetref, float* qlist, int* posqlist, float* distlist, const int& nbuf, const int& offset, const graph& g) const;

  //!
  //! \brief get_random_nodes   get random feasible nodes from graph
  //! \param g                  the graph
  //! \param start              start index in qnew
  //! \param end                end (excl.) index
  //! \param qnew               node storage (array of structs)
  //! \return
  //!
  static inline int get_random_nodes(const graph &g, const int start, const int end, float *qnew, float D, Configspace<ndof> *space);
  static inline int get_random_nodes2(const graph &g, const int start, const int end, float *qnew, float D, Configspace<ndof> *space);

  //!
  //! \brief find_neighbours
  //! \param from
  //! \param to
  //! \param qnew
  //! \param qstartp
  //! \param qendp
  //! \param poslistp
  //! \param distlistp
  //! \param nbufrest
  //! \param index
  //! \param posqlist
  //! \param numqlistleft
  //! \param numqlist
  //! \param nbuf
  //! \param offset
  //!
  inline void find_neighbours(int from, int to,
                              float *qnew, float *&qstartp, float *&qendp,
                              int *&poslistp, float *&distlistp, int &nbufrest, int &index,
                              int *posqlist, int *numqlistleft, int *numqlist,
                              int nbuf, int offset
                              ) const;


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
    static inline void calc_conn(const int *resbuf, const int *posqlist, const int *numqlistleft, const int *numqlist, int *leftconn, int *rightconn, const int from, const int to);



    //!
    //! \brief do_dijkstra  Dijkstra algorithm implementation
    //! \param g            the graph
    //! \param d            the result
    //! \param from         start index
    //! \param to           goal index
    //! \return
    //!
  int do_dijkstra(graph& g, dijkstra_result& d, int from, int to);




  //! ************************
  //! *                      *
  //! *    output methods    *
  //! *                      *
  //! ************************

  void print();
  void store_results(std::string path);
  void store_graph(std::string path, graph& g, dijkstra_result& d, int start, int end) const ;


  //for some statistics:
  static int rand_nodes_all;
  static int rand_nodes_dismissed;
  static int rand_nodes_dismissed_indicator;
  static int rand_nodes_dismissed_prob;
  static int rand_nodes_accepted;
  static int rand_block_numbers;

  static int num_nodes;
  static int num_blocks;

private:
  float H;
  float D;
  float D2;
  float factor;

  Configspace<ndof> *space;


  const int N;          //whole capacity: how many nodes can be stored
  const int blocksize;  //size of blocks

  graph graphl, graphr;

  dijkstra_result dijkstral, dijkstrar;

  int i0l, i0r; //indices of start and goal

  struct connection_data{
    float q[ndof];   //connecting node
    int index_left;  //index of connected node in graphl
    int index_right; //index of connected node in graphr
  };
  connection_data connection;


  bool connection_found;

};


#endif // PRMSOLVER3_H
