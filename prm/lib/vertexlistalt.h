#ifndef VERTEXMAP_H
#define VERTEXMAP_H

#include <map>
#include <vector>

template<int ndof, int M>
class vertexmap{

  struct key{
    int val[ndof];
  };

  //!lexikographic comparison
  struct comparer{
    bool operator() (const key& l, const key& r) const
    {
      for(int i=0;i<ndof;++i){
        if(l.val[i]<r.val[i])return true;
      }
      return false;
    }
  };

  struct block{
    //!q: array of structs
    float q[ndof*M];
    int num;

    // kanten?
  };

  typedef std::vector<block*> pvec;

  struct plist{
      pvec v;
      int act;
  };

  typedef std::map<key,plist,comparer> pmap;
  typedef pmap::iterator piterator;

public:

  vertexmap(float l_, float* mins){
    l=l_;
    linv=1.0/l;
    for(int i=0;i<ndof;++i)mins[i]=mins_[i];
  }

  ~vertexmap(){
    //nicht ganz sicher wies geht ...
  }


  void get(float* q, piterator& begin, piterator& end){
    key k;
    for(int i=0;i<ndof;++i) k.val[i]=linv*(q[i]-mins[i])-1;
    begin=map.lower_bound(k);
    for(int i=0;i<ndof;++i) k.val[i]+=2;
    end=map.upper_bound(k);
  }

  void insert(float* q){ //kanten?
    key k;
    for(int i=0;i<ndof;++i) k.val[i]=linv*(q[i]-mins[i]);
    piterator it=map.lower_bound(k);
    if(it->first!=k){
      it=map.insert(it,std::pair<key,plist>(k,plist())).first;
    }
    plist* list=&(it->second);
    pvec* v=list->v;
    int act=list->act;
    block* b=v[act];
    int num=b->num;
    float *qb=&(b->q[ndof*num]);
    for(int i=0;i<ndof;++i){
      qb[i]=q[i];
    }
    ++num;
    if(num>=M) ++(list->act);
    b->num=num;
    //kanten?

  }

  pmap map;

private:

  float l;
  float linv;
  float mins[ndof];

};

#endif // VERTEXMAP_H
