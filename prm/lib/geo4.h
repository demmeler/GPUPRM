#ifndef GEO4_H
#define GEO4_H

#include "cuda_head.h"

#include <fstream>
#include <iomanip>


#ifndef CUDA_IMPLEMENTATION
  struct float4{
    float x,y,z,w;
  };
  struct float2{
    float x,y;
  };
#endif

namespace geo4{

  class trafo4{
  public:
    float4 col[4];

    //!from DH-parameters
    qualifier trafo4(float a, float alpha, float q, float d);
    qualifier void set(float a, float alpha, float q, float d);
    qualifier trafo4(){}
    qualifier ~trafo4(){}

    //!u=T*u
    qualifier float4& apply(float4& u) const;
    //!Tu=this*u
    qualifier float4& apply(const float4& u, float4& Tu) const;
    //!u=R*u
    qualifier float4& apply_rot(float4& u) const;
    //!Ru=R*u
    qualifier float4& apply_rot(const float4& u, float4& Ru) const;
    //!u=R'*u
    qualifier float4& apply_rot_inv(float4& u) const;
    //!Rtu=R'*u
    qualifier float4& apply_rot_inv(const float4& u, float4& Rtu) const;

    //!this=this*T
    //qualifier trafo4& apply(const trafo4& T);
    //!Tres=this*T
    qualifier trafo4& apply(const trafo4& T, trafo4& Tres) const;
    //!this=T*this
    qualifier trafo4& lapply(trafo4& T);
    //!Tres=T*this (not needed)
    qualifier trafo4& lapply(const trafo4& T, trafo4& Tres) const;

    qualifier const float4& translation() const {return col[3];}

    void print(std::ostream& out, const std::string& name="");

  };


  ///   **************************
  ///   *        float2          *
  ///   *    implementations     *
  ///   **************************

  qualifier float cross2(float2& v, float2& w){
    return v.x*w.y-v.y*w.x;
  }

  qualifier void add(const float2& u, const float2& v, float2& res){
    res.x=u.x+v.x;
    res.y=u.y+v.y;
  }

  qualifier void sub(const float2& u, const float2& v, float2& res){
    res.x=u.x-v.x;
    res.y=u.y-v.y;
  }

  qualifier void normalize(float2& u){
    float factor=1.0/sqrt_(u.x*u.x+u.y*u.y);
    u.x*=factor; u.y*=factor;
  }


  ///   **************************
  ///   *        float4          *
  ///   *    implementations     *
  ///   **************************

  qualifier float4 make_float4(float x, float y, float z, float w){
    float4 f; f.x=x; f.y=y; f.z=z; f.w=w; return f;
  }

  qualifier float4 make_float4(float x, float y, float z){
    float4 f; f.x=x; f.y=y; f.z=z; return f;
  }

  qualifier float4& operator +=(float4& u, const float4& v){
    u.x+=v.x; u.y+=v.y; u.z+=v.z;
    return u;
  }

  qualifier float4& operator -=(float4& u, const float4& v){
    u.x-=v.x; u.y-=v.y; u.z-=v.z;
    return u;
  }

  qualifier float4& operator *=(float4& u, const float4& v){
    u.x*=v.x; u.y*=v.y; u.z*=v.z;
    return u;
  }

  qualifier float4& operator *=(float4& u, const float f){
    u.x*=f; u.y*=f; u.z*=f;
    return u;
  }

  qualifier float4& operator /=(float4& u, const float4& v){
    u.x/=v.x; u.y/=v.y; u.z/=v.z;
    return u;
  }

  qualifier float4& operator /=(float4& u, const float f){
    u.x/=f; u.y/=f; u.z/=f;
    return u;
  }

  qualifier float sprod(float4& u, const float4& v){
    return u.x*v.x+u.y*v.y+u.z*v.z;
  }

  qualifier void add(const float4& u, const float4& v, float4& res){
    res.x=u.x+v.x;
    res.y=u.y+v.y;
    res.z=u.z+v.z;
  }

  qualifier void sub(const float4& u, const float4& v, float4& res){
    res.x=u.x-v.x;
    res.y=u.y-v.y;
    res.z=u.z-v.z;
  }

  qualifier void lin(const float4& u, const float& b, const float4& v, float4& res){
    res.x=u.x+b*v.x;
    res.y=u.y+b*v.y;
    res.z=u.z+b*v.z;
  }

  qualifier void lin(const float& a, const float4& u, const float4& v, float4& res){
    res.x=a+u.x+v.x;
    res.y=a*u.y+v.y;
    res.z=a*u.z+v.z;
  }

  qualifier void lin(const float& a, const float4& u, const float& b, const float4& v, float4& res){
    res.x=a+u.x+b*v.x;
    res.y=a*u.y+b*v.y;
    res.z=a*u.z+b*v.z;
  }

  qualifier void normalize(float4& u){
    float factor=1.0/sqrt_(u.x*u.x+u.y*u.y+u.z*u.z);
    u.x*=factor; u.y*=factor; u.z*=factor;
  }


  ///   **************************
  ///   *        trafo4          *
  ///   *    implementations     *
  ///   **************************

  qualifier trafo4::trafo4(float a, float alpha, float q, float d){
    float cq=cos_(q);
    float sq=sin_(q);
    float ca=cos_(alpha);
    float sa=sin_(alpha);

    col[0].x=cq;    col[1].x=-sq;   col[2].x=0.0;  col[3].x=a;
    col[0].y=sq*ca; col[1].y=cq*ca; col[2].y=-sa;  col[3].y=-sa*d;
    col[0].z=sq*sa; col[1].z=cq*sa; col[2].z=ca;   col[3].z=ca*d;
  }

  qualifier void trafo4::set(float a, float alpha, float q, float d){
    float cq=cos_(q);
    float sq=sin_(q);
    float ca=cos_(alpha);
    float sa=sin_(alpha);

    col[0].x=cq;    col[1].x=-sq;   col[2].x=0.0;  col[3].x=a;
    col[0].y=sq*ca; col[1].y=cq*ca; col[2].y=-sa;  col[3].y=-sa*d;
    col[0].z=sq*sa; col[1].z=cq*sa; col[2].z=ca;   col[3].z=ca*d;
  }

  //!u=T*u
  qualifier float4& trafo4::apply(float4& u) const
  {
    float a=u.x,b=u.y,c=u.z;
    u.x=col[0].x*a+col[1].x*b+col[2].x*c+col[3].x;
    u.y=col[0].y*a+col[1].y*b+col[2].y*c+col[3].y;
    u.z=col[0].z*a+col[1].z*b+col[2].z*c+col[3].z;
    return u;
  }
  //!Tu=this*u
  qualifier float4& trafo4::apply(const float4& u, float4& Tu) const
  {
    Tu.x=col[0].x*u.x+col[1].x*u.y+col[2].x*u.z+col[3].x;
    Tu.y=col[0].y*u.x+col[1].y*u.y+col[2].y*u.z+col[3].y;
    Tu.z=col[0].z*u.x+col[1].z*u.y+col[2].z*u.z+col[3].z;
    return Tu;
  }
  //!u=R*u
  qualifier float4& trafo4::apply_rot(float4& u) const
  {
    float a=u.x,b=u.y,c=u.z;
    u.x=col[0].x*a+col[1].x*b+col[2].x*c;
    u.y=col[0].y*a+col[1].y*b+col[2].y*c;
    u.z=col[0].z*a+col[1].z*b+col[2].z*c;
    return u;
  }
  //!Ru=R*u
  qualifier float4& trafo4::apply_rot(const float4& u, float4& Ru) const
  {
    Ru.x=col[0].x*u.x+col[1].x*u.y+col[2].x*u.z;
    Ru.y=col[0].y*u.x+col[1].y*u.y+col[2].y*u.z;
    Ru.z=col[0].z*u.x+col[1].z*u.y+col[2].z*u.z;
    return Ru;
  }
  //!u=R'*u
  qualifier float4& trafo4::apply_rot_inv(float4& u) const
  {
    float a=u.x,b=u.y,c=u.z;
    u.x=col[0].x*a+col[0].y*b+col[0].z*c;
    u.y=col[1].x*a+col[1].y*b+col[1].z*c;
    u.z=col[2].x*a+col[2].y*b+col[2].z*c;
    return u;
  }
  //!Rtu=R'*u
  qualifier float4& trafo4::apply_rot_inv(const float4& u, float4& Rtu) const
  {
    Rtu.x=col[0].x*u.x+col[0].y*u.y+col[0].z*u.z;
    Rtu.y=col[1].x*u.x+col[1].y*u.y+col[1].z*u.z;
    Rtu.z=col[2].x*u.x+col[2].y*u.y+col[2].z*u.z;
    return Rtu;
  }

  //!this=this*T
  //qualifier trafo4& trafo4::apply(const trafo4& T)
  //{
  //
  //}
  //!Tres=this*T
  qualifier trafo4& trafo4::apply(const trafo4& T, trafo4& Tres) const
  {
    apply_rot(T.col[0],Tres.col[0]);
    apply_rot(T.col[1],Tres.col[1]);
    apply_rot(T.col[2],Tres.col[2]);
    apply(T.col[3],Tres.col[3]);
    return Tres;
  }
  //!this=T*this
  qualifier trafo4& trafo4::lapply(trafo4& T)
  {
    T.apply_rot(col[0]);
    T.apply_rot(col[1]);
    T.apply_rot(col[2]);
    T.apply(col[3]);
    return *this;
  }
  //!Tres=T*this (not needed)
  qualifier trafo4& trafo4::lapply(const trafo4& T, trafo4& Tres) const
  {
    T.apply_rot(col[0],Tres.col[0]);
    T.apply_rot(col[1],Tres.col[1]);
    T.apply_rot(col[2],Tres.col[2]);
    T.apply(col[3],Tres.col[3]);
    return Tres;
  }



  ///   **************************
  ///   *        output          *
  ///   *    implementations     *
  ///   **************************

  void print(const float4& v, std::ostream& out, const std::string& name=""){
    out<<std::setprecision(5);
    out<<"float4: "<<name<<"=\n"<<v.x<<"\n"<<v.y<<"\n"<<v.z<<" \n"<<std::endl;
  }

  void trafo4::print(std::ostream& out, const std::string& name){
    out<<"trafo4: "<<name<<"="<<std::endl;
    out<<std::setprecision(5);
    out<<col[0].x<<"\t"<<col[1].x<<"\t"<<col[2].x<<"\t"<<col[3].x<<std::endl;
    out<<col[0].y<<"\t"<<col[1].y<<"\t"<<col[2].y<<"\t"<<col[3].y<<std::endl;
    out<<col[0].z<<"\t"<<col[1].z<<"\t"<<col[2].z<<"\t"<<col[3].z<<std::endl<<std::endl;
  }

#define f4print(v) print(v,std::cout,#v);
#define t4print(T) T.print(std::cout,#T);


}//namespace geo4

#endif // GEO4_H
