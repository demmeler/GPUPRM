#include <iostream>

#include "lib/geo4.h"

#include "lib/kinematics.h"
#include "lib/collision4.h"

using namespace std;
using namespace geo4;
using namespace collision4;

const float pi=3.14159265358;

int main()
{

  trafo4 tp(0.0, 0.0, 0.0, 0.0);
  trafo4 tq(1.0, 0.0, 0.0, -0.01);
  trafo4 tr(-1.0, 0.0, 0.0005, 0.0);


  polytope4 P;
  generate_simplex(P, 1.0, 1.0, 1.0);

  polytope4 Q;
  generate_simplex(Q, 1.0, 1.0, 1.0);

  polytope4 R;
  generate_quader(R, 1.0, 1.0, 1.0);


#if 0
  p4print(P,tp);
  p4print(Q,tq);
  p4print(R,tr);
#endif

  printvar(seperating_vector_algorithm(P,Q,tp,tq));
  printvar(seperating_vector_algorithm(P,R,tp,tr));
  printvar(seperating_vector_algorithm(Q,R,tq,tr));


  msg("finished");
  return 0;
}

