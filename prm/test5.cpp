

#include <iostream>

#include "lib/util.hpp"
#include <map>

using namespace std;



class point{
public:
    point(int x, int y){val[0]=x;val[1]=y;}
    int val[2];
};

#define printpoint(P) printarr(P.val,2);


struct comparer{
    bool operator() (const point& l, const point& r) const
    {
        if(l.val[0]==r.val[0])return l.val[1]<r.val[1];
        else return l.val[0]<r.val[0];
    }
};

typedef map<point,int,comparer> pmap;
typedef pmap::iterator pit;

int main()
{


  pmap m;

  point p1(1,1);
  point p2(0,2);
  point p3(2,0);

  m[p1]=1;
  m[p2]=2;
  m[p3]=3;

  point pt(1,0);

  pit it=m.upper_bound(pt);


  printpoint(it->first);
  printvar(it->second);

  printvar(m.size());

  return 0;
}

