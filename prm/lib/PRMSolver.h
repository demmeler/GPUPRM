#ifndef PRMSOLVER_H
#define PRMSOLVER_H

#include "vertexlist.h"
#include "configspace.hpp"

template<int ndof>
class PRMSolver{
public:
  PRMSolver(Configspace<ndof> *cspace_){
      PRMSolver(cspace_,
                (cspace->max(0)-cspace->min(0))/500.0,
                (cspace->max(0)-cspace->min(0))/50.0);
  }

  PRMSolver(Configspace<ndof> *cspace_, float H_, float D_){
    H=H_; D=D_;
    cspace=cspace_;
    vlistfrom=new vertexlist<ndof>(H,D);
    vlistto=new vertexlist<ndof>(H,D);
  }

  ~PRMSolver(){
    delete vlistfrom, vlistto;
  }

  int reset(){
    delete vlistfrom,vlistto;
    vlistfrom=new vertexlist<ndof>(H,D);
    vlistto=new vertexlist<ndof>(H,D);
  }

  int single_querry(const float (&qstart)[ndof], const float (&qend)[ndof], const int itmax){
    vlistfrom->insert(qstart);
    vlistto->insert(qend);
    for(i=0;i<itmax;++i){
      //..........
    }
  }




private:
  vertexlist<ndof> *vlistfrom;
  vertexlist<ndof> *vlistto;

  Configspace<ndof> *cspace;

  float H,D;
};


#endif // PRMSOLVER_H
