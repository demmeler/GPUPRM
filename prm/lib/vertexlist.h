#ifndef VERTEXLIST_H
#define VERTEXLIST_H

#include <map>
#include <vector>
#include "util.hpp"

template<int ndof>
class vertexlist{

    //! *******************
    //! *    subtypes     *
    //! *******************

#if 0
  //!lexikographic comparison
  struct comparer{
    bool operator() (const vertex& l, const vertex& r) const
    {
      for(int i=0;i<ndof;++i){
        if(l.val[i]<r.val[i])return true;
      }
      return false;
    }
  };
#endif

  struct block;

  struct edgelist{
    block *blockdest;
    int *vertexdest;
    // wie kanten speichern?
  };

  struct block{
    //!array of structs
    std::vector<float> q;
    std::vector<edgelist> edgelists;
  };

  typedef std::map<int,block> pmap;
  typedef typename pmap::iterator piterator;



    //! *******************
    //! *      class      *
    //! *******************


public:

  vertexlist(float H_, float D_){
    H=H_;
    D=D_;
    D2=D*D;
    factor=1.0/H;
  }

  ~vertexlist(){
    //nicht ganz sicher wie ...
  }


  void insert(const float (&q)[ndof]){ //kanten?
    int key=calc_key(q[0]);
    piterator it = map.find(key);
    block *b;
    if(it==map.end()){
      b=&(map[key]);
      //b->q.clear();
    }else{
      b=&(it->second);
    }
    int count=b->q.size();
    b->q.resize(count+ndof);
    for(int i=0;i<ndof;++i){
      //b->q.push_back(q[i]);
      b->q[count+i]=q[i];
    }
  }


  //! get list of all vertices nearer than D
  //! qref: vertex to process
  //! qlist: buffer to store neighbour candidates, struct of arrays
  //! offset: i-th component of k-th vec in qlist is qlist[i*offset+k]
  //! nbuf: length(qlist)
  //! D: maximal distance
  inline int get_near_vertices(const float (&qref)[ndof], float* qlist, const int& nbuf, const int& offset){
    const int keylower=calc_key(qref[0]-D);
    piterator begin=map.lower_bound(keylower);
    const int keyupper=calc_key(qref[0]+D);
    piterator end=map.upper_bound(keyupper);
    int index[ndof];
    for(int i=0;i<ndof;++i){
        index[i]=i*offset;
    }
    //printvar(begin->first);
    //printvar(end->first);
    for(;!(begin==end);++begin){
      //printvar(begin->first);
      block *b=&(begin->second);
      //printvar(b->q.size()/2);
      int end=b->q.size();
      for(int k=0;k<end;k+=ndof){
        //!calc norm
        float normsq=0;
        for(int i=0;i<ndof;++i){
          float diff=b->q[k+i]-qref[i];
          normsq+=diff*diff;
        }
        //!compare norm
        if(normsq<D2){
          //printvar(k);
          //! store neighbour in buffer
          for(int i=0;i<ndof;++i){
            qlist[index[i]++]=b->q[k+i];
          }
          //! terminate if buffer full
          if(index[0]==nbuf)return nbuf;
        }
      }//for
    }//for
    //!number of written q's
    return index[0];
  }


  inline int calc_key(const float& component){
    return (int)(component*factor);
  }

  inline void set_nodes_randomly(float* qnew, int num){

  }

private:
  float H;
  float D;
  float D2;
  float factor;
  pmap map;
};

#endif // VERTEXLIST_H
