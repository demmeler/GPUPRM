#ifndef VERTEXLIST_H
#define VERTEXLIST_H

#include <map>
#include <vector>

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

  vertexlist(float D_){D=D_; D2=D*D; factor=1.0/D;}

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
    int key=(int)(qref[0]*factor);
    const float keylower=key-1;
    piterator begin=map.lower_bound(keylower);
    const float keyupper=key+1;
    piterator end=map.upper_bound(keyupper);
    int index[ndof];
    for(int i=0;i<ndof;++i){
        index[i]=i*offset;
    }
    for(;!(begin==end);++begin){
      block *b=&(begin->second);
      for(int k=0;k<b->q.size()*ndof;k+=ndof){
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
    int key=calc_key(q[0]);
    piterator it = map.find(key);
    block *b;
    if(it=map.end()){
      b=&(map[key]);
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


private:
  inline int calc_key(const float& component){return (int)(component*factor);}

  float D;
  float D2;
  float factor;
  pmap map;
};

#endif // VERTEXLIST_H
