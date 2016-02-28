#ifndef VERTEXLIST3_H
#define VERTEXLIST3_H

#include <map>
#include <vector>
#include <cmath>
#include <cfloat>
#include <stdlib.h>

#include <queue>
#include <set>

#include "configspace.hpp"
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

  typedef std::map<int,block*> pmap;
  typedef typename pmap::iterator piterator;

  struct graph{
    pmap map; //map for sorting on high level
    std::vector<block> blocks;

    std::vector<float> qstorage;                 //length ndof*N
    std::vector<int> surrnum;                    //length N
    std::vector<std::vector<int>> edgelists;     //length N
    std::vector<std::vector<float>> edgeweights; //length N

    int newblockpos;  //position of next block in size-N-arrays
    int blocknum;     //number of used blocks
  };

    //! *******************
    //! *      class      *
    //! *******************


public:

  vertexlist(float H_, float D_, Configspace<ndof> *space_=0x0)
  {
    graphl.qstorage.resize(ndof*N);
    graphl.surrnum.resize(N);
    graphl.edgelists.resize(N);
    graphl.edgeweights.resize(N);
    graphl.blocks.resize(N/blocksize);

    graphl.newblockpos=0;
    graphl.blocknum=0;

    graphr.qstorage.resize(ndof*N);
    graphr.surrnum.resize(N);
    graphr.edgelists.resize(N);
    graphr.edgeweights.resize(N);
    graphr.blocks.resize(N/blocksize);

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

  ~vertexlist(){
    //nicht ganz sicher wie ...
  }

  int init(const float* qstart, const float* qend){
    i0l=insert(qstart,graphl);
    i0r=insert(qend,graphr);
  }

  //!insert node q, data of surrnum and edges are not modified -> position returned for this
  int insert(const float* q){return insert(q,graphl);}
  int insert(const float* q, graph& g){
    return insert(q,1,g);
  }

  //!insert node q, data of surrnum and edges are not modified -> position returned for this
  int insert(const float* q, const int offset, graph& g){
    int key=calc_key(q[0]);
    piterator it = g.map.find(key);
    block *b;
    int fall=0;
    int size,size2,pos,num;
    if(it==g.map.end()){
      size2=g.blocknum;
      size=g.newblockpos;

      b=&(g.blocks[g.blocknum++]);
      g.map[key]=b;
      b->pos=g.newblockpos;
      g.newblockpos+=blocksize;
      if(g.newblockpos>N) return -1;
      b->num=0;
    }else{
      size=g.newblockpos;
      size2=g.blocknum;
      fall=1;

      b=it->second;
      num=b->num;
      pos=b->pos;
      while(b->num>=blocksize){
        b=b->next;
      }
    }
    //printvar(b->pos);
    //printvar(b->num);
    int position=b->pos+b->num++;
    int qposition=ndof*position;
    for(int i=0;i<ndof;++i){
      g.qstorage[qposition+i]=q[offset*i];
    }
    //surrnum[position]=0;
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
  inline int get_near_vertices(const float* qref, float* qlist, const int& nbuf, const int& offset){
    return get_near_vertices(qref,qlist,nbuf,offset,graphl);
  }

  inline int get_near_vertices(const float* qref, float* qlist, const int& nbuf, const int& offset, graph& g){
    const int keylower=calc_key(qref[0]-D);
    piterator begin=g.map.lower_bound(keylower);
    const int keyupper=calc_key(qref[0]+D);
    piterator end=g.map.upper_bound(keyupper);
    int index[ndof];
    for(int i=0;i<ndof;++i){
        index[i]=i*offset;
    }
    //printvar(begin->first);
    //printvar(end->first);
    for(;!(begin==end);++begin){
      //printvar(begin->first);
      block *b=begin->second;
      //printvar(b->pos);
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
            //printvar(diff);
            normsq+=diff*diff;
          }
          //!compare norm
          if(normsq<D2){
            //printvar(k);
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

  inline int get_near_vertices(const float* qref, const int& offsetref, float* qlist, int* posqlist, float* distlist, const int& nbuf, const int& offset, graph& g){
    const int keylower=calc_key(qref[0]-D);
    piterator begin=g.map.lower_bound(keylower);
    const int keyupper=calc_key(qref[0]+D);
    piterator end=g.map.upper_bound(keyupper);
    int index[ndof];
    for(int i=0;i<ndof;++i){
        index[i]=i*offset;
    }
    //printvar(begin->first);
    //printvar(end->first);
    for(;!(begin==end);++begin){
      //printvar(begin->first);
      block *b=begin->second;
      //printvar(b->pos);
      bool more=true;
      while(more){
        more=b->num>=blocksize;
        //printvar(b->num);
        const int pos=b->pos;
        const int max=pos+b->num;
        for(int l=pos;l<max;++l){
          int k=ndof*l;
          //!calc norm
          float normsq=0;
          for(int i=0;i<ndof;++i){
            float diff=g.qstorage[k+i]-qref[offsetref*i];
            //printvar(diff);
            normsq+=diff*diff;
          }
          //!compare norm
          if(normsq<D2){
            //printvar(k);
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

  //! get list of all vertices nearer than D
  //! qlist: buffer to store neighbour candidates, struct of arrays
  //! offset: i-th component of k-th vec in qlist is qlist[i*offset+k]
  //! nbuf: length(qlist)
  //! D: maximal distance
  //! return: 0 - no errors, 1 - connection found
  inline int processing_step(float* qnew, const int num, float* qlist, int* resbuf, const int nbuf, const int offset){
    //!
    //! create nodes qnew randomly
    //!

    //!from left graph
    for(int j=0;j<num/2;++j){
      float qnewtemp[ndof];
      do{
        //!choose random block
        int k;
        block *b;
        do{
            k=rand()%graphl.blocknum;
            b=&(graphl.blocks[k]);
        }while(b->num==0);
        //!chose random vertex
        int l=ndof*(b->pos+rand()%b->num);
        for(int i=0;i<ndof;++i){
          float x=graphl.qstorage[l+i];
          qnew[j+num*i]=graphl.qstorage[l+i]-D+2*D*((float)rand()/RAND_MAX);
          float y=qnew[j+num*i];
          qnewtemp[i]=qnew[j+num*i];
        }
#ifndef NO_IO
        printvar(l);
        printarr(qnewtemp,ndof);
#endif
      }while(space->indicator(&qnewtemp[0])!=0);
    }
    //!from right graph
    for(int j=num/2;j<num;++j){
      float qnewtemp2[ndof];
      do{
        //!choose random block
        int k;
        block *b;
        do{
            k=rand()%graphr.blocknum;
            b=&(graphr.blocks[k]);
        }while(b->num==0);
        //!chose random vertex
        int l=ndof*(b->pos+rand()%b->num);
        for(int i=0;i<ndof;++i){
          qnew[j+num*i]=graphr.qstorage[l+i]-D+2*D*((float)rand()/RAND_MAX);
          qnewtemp2[i]=qnew[j+num*i];
        }
#ifndef NO_IO
        printvar(l);
        printarr(qnewtemp2,ndof);
#endif
      }while(space->indicator(&qnewtemp2[0])!=0);
    }



    //!
    //! check if qnew is free
    //!

    // ....... --> call indicator function on cpu, wait
    //no, already done in set_nodes_randomly

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
    float *qlistp=qlist;
    int *poslistp=poslist;
    float *distlistp=distlist;
    int nbufrest=nbuf;

    for(int j=0;j<num;++j){
      posqlist[j]=index;

      int writtenleft=get_near_vertices(&qnew[j],num,qlistp,poslistp,distlistp,nbufrest,offset,graphl);

      numqlistleft[j]=writtenleft;

      if(writtenleft>=nbufrest){
        for(int l=j+1;l<num;++l){
          posqlist[l]=nbuf;
          numqlist[l]=0;
          numqlistleft[l]=0;
        }
        numqlist[j]=writtenleft;
        break;
      }


      qlistp+=writtenleft;
      poslistp+=writtenleft;
      distlistp+=writtenleft;
      nbufrest-=writtenleft;
      index+=writtenleft;

      int writtenright=get_near_vertices(&qnew[j],num,qlistp,poslistp,distlistp,nbufrest,offset,graphr);

      if(writtenright>=nbufrest){
        for(int l=j+1;l<num;++l){
          posqlist[l]=nbuf;
          numqlist[l]=0;
          numqlistleft[l]=0;
        }
        numqlist[j]=writtenleft+writtenright;
        break;
      }

      qlistp+=writtenright;
      poslistp+=writtenright;
      distlistp+=writtenright;
      nbufrest-=writtenright;
      index+=writtenright;

      numqlist[j]=writtenleft+writtenright;
    }

    Nqlist=index;


    //!
    //! edge tests between new nodes neglected
    //!

    //     ---

    //!
    //! calculate which edges exist
    //!
#ifndef NO_IO
    printarr(qnew,ndof*num);
    printvar(num);
    printarr(qlist,ndof*nbuf);
    printarr(posqlist,num);
    printarr(numqlist,num);
    printarr(numqlistleft,num);
    printvar(offset);
    printarr(distlist,nbuf);
    printarr(poslist,nbuf);
#endif
    //...... --> call indicator function on GPU
    space->indicator2(qnew,num,qlist,resbuf,posqlist,numqlist,offset);

#ifndef NO_IO
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

      bool leftconn=false;
      for(int i=min;i<maxleft;++i){
        if(resbuf[i]==0){
          leftconn=true;
          break;
        }
      }
      if(leftconn){
        //!
        //! connection to left graph exists -> insert in left graph
        //!
        int position=insert(&qnew[j],num,graphl);
        int surrnump=0;
        std::vector<int> *v=&(graphl.edgelists[position]);
        std::vector<float> *w=&(graphl.edgeweights[position]);
        for(int i=min;i<maxleft;++i){
          if(resbuf[i]==0){
            int goalpos=poslist[i];
            float dist=distlist[i];
            ++surrnump;
            v->push_back(goalpos);
            w->push_back(dist);
            ++(graphl.surrnum[goalpos]);
            graphl.edgelists[goalpos].push_back(position);
            graphl.edgeweights[goalpos].push_back(dist);
          }
        }
        graphl.surrnum[position]+=surrnump;
        for(int i=maxleft;i<max;++i){
          if(resbuf[i]==0){
            //!
            //!  Connection found! abort
            //!
            for(int l=0;l<ndof;++l)connection.q[l]=qnew[j+num*l];
            float minnorm=w->at(0);
            int minpos=0;
            for(int j=0;j<w->size();++j){
              if(minnorm>w->at(j)){minnorm=w->at(j);minpos=j;}
            }
            connection.index_left=v->at(minpos);
            connection.index_right=poslist[i];

            int res0=do_dijkstra(graphl,dijkstral,i0l,connection.index_left);
            int res1=do_dijkstra(graphr,dijkstrar,i0r,connection.index_right);
            if(res0==0){msg("ERROR: no path found by dijkstra in graphl");}
            if(res1==0){msg("ERROR: no path found by dijkstra in graphr");}

            return 1;
          }
        }
      }else{
        bool rightconn=false;
        for(int i=maxleft;i<max;++i){
          if(resbuf[i]==0){
            rightconn=true;
            break;
          }
        }
        if(rightconn){
          //!
          //! connection to right graph exists -> insert in right graph
          //!
          int position=insert(&qnew[j],num,graphr);
          int surrnump=0;
          std::vector<int> *v=&(graphr.edgelists[position]);
          std::vector<float> *w=&(graphr.edgeweights[position]);
          for(int i=maxleft;i<max;++i){
            if(resbuf[i]==0){
              int goalpos=poslist[i];
              float dist=distlist[i];
              ++surrnump;
              v->push_back(goalpos);
              w->push_back(dist);
              ++(graphr.surrnum[goalpos]);
              graphr.edgelists[goalpos].push_back(position);
              graphr.edgeweights[goalpos].push_back(dist);
            }
          }
          graphr.surrnum[position]+=surrnump;
          for(int i=min;i<maxleft;++i){
            if(resbuf[i]==0){
              //!
              //!  Connection found! abort
              //!
              for(int l=0;l<ndof;++l)connection.q[l]=qnew[j+num*l];
              float minnorm=w->at(0);
              int minpos=0;
              for(int j=0;j<w->size();++j){
                if(minnorm>w->at(j)){minnorm=w->at(j);minpos=j;}
              }
              connection.index_right=v->at(minpos);
              connection.index_left=poslist[i];

              int res0=do_dijkstra(graphl,dijkstral,i0l,connection.index_left);
              int res1=do_dijkstra(graphr,dijkstrar,i0r,connection.index_right);
              if(res0==0){msg("ERROR: no path found by dijkstra in graphl");}
              if(res1==0){msg("ERROR: no path found by dijkstra in graphr");}

              return 1;
            }
          }//for
        }
      }

    }//for


    return 0;

  }




  //!for dijkstra at the end
  struct dijkstra_result{
#if 0
    struct node {float dist; int pos;};
    struct comp {
      inline bool operator() (const node& n1, const node& n2) const {
        return n1.dist>n2.dist;
      }
    };
    std::priority_queue<node,std::vector<node>,comp> queue;
#endif

    std::vector<int> parent; //length newblockpos
    std::vector<float> dist; //length newblockpos
    std::vector<int> path;
  };



  int do_dijkstra(graph& g, dijkstra_result& d, int from, int to){
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
      }

    }//while
    return 0; //no path
  }





  inline int calc_key(const float& component){
    return (int)(component*factor);
  }

  void print(){
    msg("-------vertexlist-------");
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

  void store_graphs(std::string path){
    std::string pathl=path+"/graphl";
    //system("rm -rf "+pathl);
    system(("mkdir "+pathl).c_str());
    store_graph(pathl,graphl,i0l,connection.index_left);

    std::string pathr=path+"/graphr";
    //system("rm -rf "+pathr);
    system(("mkdir "+pathr).c_str());
    store_graph(pathr,graphr,i0r,connection.index_right);

    write_file(path+"/connection.bin",(int*)&connection_found,1);
  }

  void store_graph(std::string path, graph& g, int start, int end) const {
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

  }


private:
  float H;
  float D;
  float D2;
  float factor;

  Configspace<ndof> *space;


  const int N=1024*1024;   //whole capacity: how many nodes can be stored
  const int blocksize=256; //size of blocks

  graph graphl, graphr;

  dijkstra_result dijkstral, dijkstrar;

  int i0l, i0r; //indices of start and goal

  struct {
    float q[ndof];   //connecting node
    int index_left;  //index of connected node in graphl
    int index_right; //index of connected node in graphr
  }connection;


  bool connection_found;

};


#endif // VERTEXLIST3_H
