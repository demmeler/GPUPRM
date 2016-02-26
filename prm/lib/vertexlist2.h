#ifndef VERTEXLIST2_H
#define VERTEXLIST2_H

#include <map>
#include <vector>
#include "util.hpp"

template<int ndof>
class vertexlist{

    //! *******************
    //! *    subtypes     *
    //! *******************

    struct key{
        int index[ndof];
    };

#if 1
  //!lexikographic comparison
  struct comparer{
    inline bool operator() (const key& l, const key& r) const
    {
      for(int i=0;i<ndof;++i){
        int diff=l.index[i]-r.index[i];
        if(diff<0)return true;
        else if(diff>0)return false;
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
  //public:
    //block():q(100),edgelists(100){count=0;}
    //!array of structs
    std::vector<float> q;
    std::vector<edgelist> edgelists;
    //int count;
  };

  typedef std::map<key,block,comparer> pmap;
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


  //! get list of all vertices nearer than D
  //! qref: vertex to process
  //! qlist: buffer to store neighbour candidates, struct of arrays
  //! offset: i-th component of k-th vec in qlist is qlist[i*offset+k]
  //! nbuf: length(qlist)
  //! D: maximal distance
  inline int get_near_vertices(const float (&qref)[ndof], float* qlist, const int nbuf, const int offset){
    const key keylower=calc_key_dist(qref,-D);
    piterator begin=map.lower_bound(keylower);
    const key keyupper=calc_key_dist(qref,D);
    piterator end=map.upper_bound(keyupper);
    int index[ndof];
    for(int i=0;i<ndof;++i){
        index[i]=i*offset;
    }

    //printarr(begin->first.index,ndof);
    //printarr(end->first.index,ndof);

    for(;!(begin==end);++begin){
      //printarr(begin->first.index,ndof);
      block *b=&(begin->second);
      //printvar(b->q.size()/2);
      for(int k=0;k<b->q.size();k+=ndof){
        //!calc norm
        float normsq=0;
        float qakt[ndof];
        for(int i=0;i<ndof;++i){
          qakt[i]=b->q[k+i];
          float diff=qakt[i]-qref[i];
          normsq+=diff*diff;
        }
        //!compare norm
        if(normsq<D2){
          //printvar(k);
          //! store neighbour in buffer
          for(int i=0;i<ndof;++i){
            qlist[index[i]++]=qakt[i];
          }
          //! terminate if buffer full
          if(index[0]==nbuf)return nbuf;
        }
      }//for
    }//for
    //!number of written q's
    return index[0];
  }

  void insert(const float (&q)[ndof]){ //kanten?
    //printvar(q[0]);
    key k=calc_key(q);
    //printarr(k.index,ndof);
    piterator it = map.find(k);
    block *b;
    bool end=false;
    if(it==map.end()){
      end=true;
      b=&(map[k]);
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

  inline key calc_key(const float (&q)[ndof]) const{
    key k;
    for(int i=0;i<ndof;++i) k.index[i]=(int)(factor*q[i]);
    return k;
  }
  inline key calc_key_dist(const float (&q)[ndof], const float& D) const{
    key k;
    for(int i=0;i<ndof;++i) k.index[i]=(int)(factor*(q[i]+D));
    return k;
  }

private:
  float H;
  float D;
  float D2;
  float factor;
  pmap map;
};

#endif // VERTEXLIST2_H
