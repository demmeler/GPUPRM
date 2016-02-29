#include "arrayconfigspace.hpp"
#include <math.h>

ArrayConfigspace::ArrayConfigspace(const int* array_, int b_, int h_, float minx, float maxx, float miny, float maxy){
  b=b_; h=h_; n=b*h;
  mins[0]=minx;
  mins[1]=miny;
  maxs[0]=maxx;
  maxs[1]=maxy;
  factor[0]=b/(maxs[0]-mins[0]);
  factor[1]=h/(maxs[1]-mins[1]);
  dq=1.0/factor[0];
  if(1.0/factor[1]<dq)dq=1.0/factor[1];
  array=new int[n];
  for(int i=0;i<n;++i){
      array[i]=array_[i];
  }
}
ArrayConfigspace::~ArrayConfigspace(){
  delete array;
}

//!initialization function
int ArrayConfigspace::init(){
  return 0;
}
//! indicator function of obstacles
//! q: length d*N, array of structures: q[N*k+i]= k-th component of i-th q-vector
//! res: length N
int ArrayConfigspace::indicator(const float* q, int* res, const int N, const int offset){
  int iy=offset;
  for(int ix=0;ix<N;++ix,++iy){
    int x=(int)((q[ix]-mins[0])*factor[0]);
    int y=(int)((q[iy]-mins[1])*factor[1]);
    int resix=-1;
    if(x>=0 && x<b && y>=0 && y<h){
      int index=b*y+x;
      resix=array[index];
    }
    res[ix]=resix;
  }
  return 0;
}
//!case N=1
int ArrayConfigspace::indicator(const float* q){
  int res=0;
  indicator(q,&res,1,1);
  return res;
}

//! checks if indicator function = 1 somewhere on the line between qs and qe
//! res is return value
int ArrayConfigspace::indicator2(const float* qs, const float* qe, int *res, const int N, const int offset){
  for(int k=0;k<N;++k){
    float qs_[]={qs[k],qs[k+offset]};
    float qe_[]={qe[k],qe[k+offset]};
    if(0 != check_boundaries(&qs_[0]) || 0 != check_boundaries(&qe_[0])){
      res[k]=2;
    }else{
      float qsx=qs_[0], qsy=qs_[1], qex=qe_[0], qey=qe_[1];
      float dqx=qex-qsx, dqy=qey-qsy;
      float dist=sqrt(dqx*dqx+dqy*dqy);
      int n=(int)(dist/dq);
      if(n==0)n=1;
      float nqx=dqx/n, nqy=dqy/n;
      res[k]=0;
      for(int i=0;i<n;++i){
        float q[2]={qsx+nqx*i,qsy+nqy*i};
        int x=(int)((q[0]-mins[0])*factor[0]);
        int y=(int)((q[1]-mins[1])*factor[1]);
        int index=b*y+x;
        if(0 != array[index]){
          res[k]=1;
          break;
        }
      }//for
    }
  }//for
  return 0;
}


//! same paircheck as above, but with compressed storage:
//! checks pairs: (qs[i],...) ->  (qe(posqe[i]),...) , ...., (qe[posqe[i]+numqe[i]-1],...) for i=0,...,M-1
int ArrayConfigspace::indicator2(const float* qs, const int M, const float* qe, int *res, const int *posqe, const int *numqe, const int offset){
  for(int k=0;k<M;++k)
  for(int l=posqe[k];l<posqe[k]+numqe[k];++l){
    float qs_[]={qs[k],qs[k+M]};
    float qe_[]={qe[l],qe[l+offset]};
    if(0 != check_boundaries(&qs_[0]) || 0 != check_boundaries(&qe_[0])){
      res[l]=2;
    }else{
      float qsx=qs_[0], qsy=qs_[1], qex=qe_[0], qey=qe_[1];
      float dqx=qex-qsx, dqy=qey-qsy;
      float dist=sqrt(dqx*dqx+dqy*dqy);
      int n=(int)(dist/dq);
      if(n==0)n=1;
      float nqx=dqx/n, nqy=dqy/n;
      res[l]=0;
      for(int i=0;i<n;++i){
        float q[2]={qsx+nqx*i,qsy+nqy*i};
        int x=(int)((q[0]-mins[0])*factor[0]);
        int y=(int)((q[1]-mins[1])*factor[1]);
        int index=b*y+x;
        if(0 != array[index]){
          res[l]=1;
          break;
        }
      }//for
    }
  }//for
  return 0;
}

//! structure like indicator function
//! returns if lies in boundaries
int ArrayConfigspace::check_boundaries(const float* q, int* res, int N, int offset){
  int iy=offset;
  for(int ix=0;ix<N;++ix,++iy){
    float xf=q[ix];
    float yf=q[iy];
    int resix=1;
    if(xf>=mins[0] && xf<maxs[0] && yf>=mins[1] && yf<maxs[1]) resix=0;
    res[ix]=resix;
  }
  return 0;
}
//! for N=1, res is return value
int ArrayConfigspace::check_boundaries(const float* q){
  float xf=q[0];
  float yf=q[1];
  int res;
  if(xf>=mins[0] && xf<maxs[0] && yf>=mins[1] && yf<maxs[1]) res=0;
  else res=1;
  return res;
}


//! minimal value of qi, i=0...d-1
float ArrayConfigspace::min(int i){
  return mins[i];
}
//! maximal value of qi, i=0...d-1
float ArrayConfigspace::max(int i){
  return maxs[i];
}
