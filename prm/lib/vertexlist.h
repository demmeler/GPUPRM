#ifndef VERTEXMAP_H
#define VERTEXMAP_H

#include <map>
#include <vector>

template<int ndof, int M>
class vertexmap{

  struct vertex{
    float val[ndof];
  };

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

  struct edgelist{

    // kanten?
  };

  typedef std::vector<block*> pvec;

  struct plist{
      pvec v;
      int act;
  };

  typedef std::map<float,edgelist,comparer> pmap;
  typedef pmap::iterator piterator;

public:

  vertexmap(double* D_){
    D=D_;
  }

  ~vertexmap(){
    //nicht ganz sicher wie ...
  }


  void get(float* q, piterator& begin, piterator& end){
    vertex q;
    float q0=q[0];
    q.val[0]=q0-D;
    for(int i=1;i<ndof;++i) q.val[i]=q[i];
    begin=map.lower_bound(q);
    q.val[0]=q0+D;
    end=map.upper_bound(q);
  }

  void insert(float* q){ //kanten?
    vertex *v=new vertex;
    for(int i=0;i<ndof;++i) v->val[i]=q[i];
    map[*v]=
  }

  pmap map;

private:
  float D;

};

#endif // VERTEXMAP_H
