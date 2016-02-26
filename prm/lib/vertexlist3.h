#ifndef VERTEXLIST3_H
#define VERTEXLIST3_H

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

  struct block{
    int pos; //! position of first vertex
    int num; //! current number of vertices stored
    block* next; //!if num==blocksize -> pointer to next block
  };

  typedef std::map<int,block> pmap;
  typedef typename pmap::iterator piterator;

  struct graph{
    pmap map; //map for sorting on high level

    std::vector<float> qstorage;              //length ndof*N
    std::vector<int> surrnum;                 //length N
    std::vector<std::vector<int>> edgelists;  //length N

    const int N=1024*1024;   //whole capacity: how many nodes can be stored
    const int blocksize=256; //size of blocks
    int newblockpos;  //number of used blocks
  };



    //! *******************
    //! *      class      *
    //! *******************


public:

  vertexlist(float H_, float D_)
  {
    qstorage.resize(ndof*N);
    surrnum.resize(N);
    edgelists.resize(N);
    H=H_;
    D=D_;
    D2=D*D;
    factor=1.0/H;
    newblockpos=0;
  }

  ~vertexlist(){
    //nicht ganz sicher wie ...
  }


  //!insert node q, data of surrnum and edges are not modified -> position returned for this
  int insert(const float* q, graph& g){
    int key=calc_key(q[0]);
    piterator it = map.find(key);
    block *b;
    if(it==map.end()){
      b=&(map[key]);
      if(newblockpos>=N) return -1;
      b->pos=newblockpos;
      newblockpos+=blocksize;
      b->num=0;
    }else{
      b=&(it->second);
      while(b->num>=blocksize){
        b=b->next;
      }
    }
    int position=b->pos+b->num++;
    int qposition=ndof*position;
    for(int i=0;i<ndof;++i){
      qstorage[qposition+i]=q[i];
    }
    //surrnum[position]=0;
    if(b->num>=blocksize){
        b->next=new block;
        block *bnew=b->next;
        bnew->pos=newblockpos;
        newblockpos+=blocksize;
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
      //printvar(b->pos);
      bool more=true;
      while(more){
        more=b->num>=blocksize;
        //printvar(b->num);
        const int pos=ndof*b->pos;
        const int max=pos+ndof*b->num;
        for(int k=pos;k<max;k+=ndof){
          //!calc norm
          float normsq=0;
          for(int i=0;i<ndof;++i){
            float diff=qstorage[k+i]-qref[i];
            //printvar(diff);
            normsq+=diff*diff;
          }
          //!compare norm
          if(normsq<D2){
            //printvar(k);
            //! store neighbour in buffer
            for(int i=0;i<ndof;++i){
              qlist[index[i]++]=qstorage[k+i];
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

#if 1
  //! get list of all vertices nearer than D
  //! qlist: buffer to store neighbour candidates, struct of arrays
  //! offset: i-th component of k-th vec in qlist is qlist[i*offset+k]
  //! nbuf: length(qlist)
  //! D: maximal distance
  inline int processing_step(float* qnew, const int num, float* qlist, int* resbuf, const int nbuf, const int offset){
    //!
    //! create random node qref
    //!

    set_nodes_randomly(qnew, num);

    //!
    //! check if qnew is free
    //!

    // ....... --> call indicator function on cpu, wait

    //!
    //! make list of potential neighbours for all new nodes
    //!

    int index[ndof];
    for(int i=0;i<ndof;++i){
        index[i]=i*offset;
    }

    //!posqlist[i]: position in qlist buffer, at which neighbours of qnew[i] start
    int *poslist=new int[nbuf];
    int posqlist[num];
    int numqlist[num];
    int Nqlist;             //sum(numqlist)
    bool buffer_full=false; //abort flag if buffer is full (should not happen normally)

    for(int j=0;j<num && !buffer_full;++j){
      int p=ndof*j;
      posqlist[j]=index[0];

      const int keylower=calc_key(qnew[p]-D);
      piterator begin=map.lower_bound(keylower);
      const int keyupper=calc_key(qnew[p]+D);
      piterator end=map.upper_bound(keyupper);

      //printvar(begin->first);
      //printvar(end->first);
      for(;!(begin==end)&&!buffer_full;++begin){
        //printvar(begin->first);
        block *b=&(begin->second);
        //printvar(b->pos);
        bool more=true;
        while(more){
          const int bnum=b->num;
          if(index[0]+bnum>nbuf){
            buffer_full=true;
            for(int l=j+1;l<num;++l){
              posqlist[l]=nbuf;
              numqlist[l]=0;
            }
            break;
          }
          more=bnum>=blocksize;
          //printvar(b->num);
          const int pos=b->pos;
          const int max=pos+bnum;
          for(int l=pos;l<max;l+=ndof){
            int k=ndof*l;
            //!calc norm
            float normsq=0;
            for(int i=0;i<ndof;++i){
              float diff=qstorage[k+i]-qnew[p+i];
              //printvar(diff);
              normsq+=diff*diff;
            }
            //!compare norm
            if(normsq<D2){
              //printvar(k);
              //! store neighbour in buffer
              poslist[index[0]]=l;
              for(int i=0;i<ndof;++i){
                qlist[index[i]++]=qstorage[k+i];
              }
            }
          }//for
          b=b->next;
        }
      }//for

      numqlist[j]=index[0]-posqlist[j];

    }

    Nqlist=index[0];


    //!
    //! edge tests between new nodes neglected
    //!

    //     ---

    //!
    //! calculate which edges exist
    //!


    //...... --> call indicator function on GPU
    // => resbuf= ......


    //!
    //! insert nodes and edges
    //!

    for(int j=0;j<num;++j){
      int p=ndof*j;

      int position=insert(&q[p]);
      int surrnump=0;
      std::vector<int> *v=&(edgelists[position]);
      int max=posqlist[i]+numqlist[j];
      for(int i=posqlist[i];i<max;++i){
        if(resbuf[i]==0){
          int goalpos=poslist[i];
          ++surrnump;
          v->push_back(goalpos);
          ++(surrnum[goalpos]);
          edgelists[goalpos].push_back(position);
        }
      }
      surrnum[position]+=surrnump;

      //! insert edges in dependance of resbuf

    }
  }
#endif

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

  graph graphl, graphr;

};


#endif // VERTEXLIST3_H
