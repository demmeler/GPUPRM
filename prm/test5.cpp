

#include <iostream>

#include "lib/util.hpp"
#include <map>
#include <vector>

#include <math.h>

#include <chrono>
#include <cstdlib>
#include <ctime>


using namespace std;

int pow(int x, int n){
    if(n>1)return x*pow(x,n-1);
    else if(n==1) return x;
    else if(n==0) return 1;
}

class point{
public:
    point(const point& p){
        val[0]=p.val[0];
        val[1]=p.val[1];
    }
    point(){}
    point(int x, int y){val[0]=x;val[1]=y;}
    int val[2];
};

#define printpoint(P) printarr(P.val,2);
#define many(x) x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x \
    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x


struct comparer{
    bool operator() (const point& l, const point& r) const
    {
        if(l.val[0]==r.val[0])return l.val[1]<r.val[1];
        else return l.val[0]<r.val[0];
    }
};

typedef map<point,point,comparer> pmap;
typedef pmap::iterator pit;

int main()
{

#if 0
  pmap m;

  point p1(1,1);
  point p2(0,2);
  point p3(2,0);

  m[p1]=p1;
  m[p2]=p2;
  m[p3]=p3;

  p2.val[1]=7;

  point pt(0,1);

  pit it=m.upper_bound(pt);


  printpoint(it->first);
  printpoint(it->second);

  printvar(m.size());
#endif

  clock_t t0=clock();

  int n=pow(2,19);

  vector<float> v(n);

  clock_t t1=clock();

  v[n-1]=7;

  printvar(n);
  printvar(v[n-1]);

  srand(clock());

  for(int i=0;i<n;++i){
      v[i]=rand();
  }


  clock_t t2=clock();
  std::chrono::high_resolution_clock::time_point t2c=std::chrono::high_resolution_clock::now();

  float x=0;
  for(int i=0;i<n;++i){
      float vi=rand();//v[i];
      x+=vi*vi many(*vi*vi);
      //x+=fabs(vi) many(+fabs(vi));
  }

  clock_t t3=clock();
  std::chrono::high_resolution_clock::time_point t3c=std::chrono::high_resolution_clock::now();
  int delta23=std::chrono::duration_cast<std::chrono::milliseconds>(t3c-t2c).count();


  printvar(x);

  printvar((float)(t1-t0)/CLOCKS_PER_SEC);
  printvar((float)(t2-t1)/CLOCKS_PER_SEC);
  printvar((float)(t3-t2)/CLOCKS_PER_SEC);
  printvar(delta23);



  return 0;
}

