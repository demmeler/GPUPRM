#include <iostream>

#include "lib/collision.hpp"
#include "lib/util.hpp"

using namespace std;

using namespace collision;

int main()
{
  trafo tp(0.50*3.141592653589793,0.0,0.0,0.0);
  trafo tq(0.26*3.1415,0.0,0.5,0.65);
  trafo tr(0.00*3.1415,0.0,0.5,0.0);

  Polytope P;

  double vertp[]={0.0,0.0,0.0,
                  1.0,0.0,0.0,
                  0.0,1.0,0.0,
                  0.0,0.0,1.0};
  int dspp[]={0,3,6,9};
  int cntp[]={3,3,3,3};
  int destp[]={1,2,3,
               0,2,3,
               0,1,3,
               0,1,2};

  P.n=4;
  P.m=12;
  P.vertices=&(vertp[0]);
  P.dsp=&(dspp[0]);
  P.cnt=&(cntp[0]);
  P.dest=&(destp[0]);



  Polytope Q;

  double vertq[]={0.0,0.0,0.0,
                  1.0,0.0,0.0,
                  0.0,1.0,0.0,
                  0.0,0.0,1.0};
  int dspq[]={0,3,6,9};
  int cntq[]={3,3,3,3};
  int destq[]={1,2,3,
               0,2,3,
               0,1,3,
               0,1,2};

  Q.n=4;
  Q.m=12;
  Q.vertices=&(vertq[0]);
  Q.dsp=&(dspq[0]);
  Q.cnt=&(cntq[0]);
  Q.dest=&(destq[0]);

  Polytope R;

  double vertr[]={0.0,0.0,0.0,
                  1.0,0.0,0.0,
                  0.0,1.0,0.0,
                  1.0,1.0,0.0,
                  0.0,0.0,1.0,
                  1.0,0.0,1.0,
                  0.0,1.0,1.0,
                  1.0,1.0,1.0};
  int dspr[]={0,3,6,9,12,15,18,21};
  int cntr[]={3,3,3,3,3,3,3,3};
  int destr[]={1,3,4,
               0,2,5,
               1,3,6,
               0,2,7,
               5,7,0,
               4,6,1,
               5,7,2,
               4,6,3};

  R.n=8;
  R.m=24;
  R.vertices=&(vertr[0]);
  R.dsp=&(dspr[0]);
  R.cnt=&(cntr[0]);
  R.dest=&(destr[0]);



  cout << "inited" <<endl;

  bool coll=seperating_vector_algorithm(P,R,tp,tr);

  cout << "Collision?? -> "<< coll << endl;
  return 0;
}

