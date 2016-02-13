#ifndef GEO4_H
#define GEO4_H

#include <math.h>
#include <fstream>
#include <iomanip>

namespace geo4{

  struct float4{
    float x,y,z,w;
  };

  class trafo4{
  public:
    float4 col[4];

    //!from DH-parameters
    trafo4(float a, float alpha, float q, float d);
    trafo4(){}
    ~trafo4(){}

    //!u=T*u
    inline float4& apply(float4& u) const;
    //!Tu=this*u
    inline float4& apply(const float4& u, float4& Tu) const;
    //!u=R*u
    inline float4& apply_rot(float4& u) const;
    //!Ru=R*u
    inline float4& apply_rot(float4& u, float4& Ru) const;
    //!u=R'*u
    inline float4& apply_rot_inv(float4& u) const;
    //!Rtu=R'*u
    inline float4& apply_rot_inv(float4& u, float4& Rtu) const;

    //!this=this*T
    //inline trafo4& apply(const trafo4& T);
    //!Tres=this*T
    inline trafo4& apply(const trafo4& T, trafo4& Tres) const;
    //!this=T*this
    inline trafo4& lapply(trafo4& T);
    //!Tres=T*this (not needed)
    inline trafo4& lapply(const trafo4& T, trafo4& Tres) const;

    inline float4& t(){return col[3];}

    inline void print(std::ostream& out, const std::string& name="");

  };



  ///   **************************
  ///   *        float4          *
  ///   *    implementations     *
  ///   **************************

  inline float4 make_float4(float x, float y, float z, float w){
    float4 f; f.x=x; f.y=y; f.z=z; f.w=w; return f;
  }

  inline float4 make_float4(float x, float y, float z){
    float4 f; f.x=x; f.y=y; f.z=z; return f;
  }

  inline float4& operator +=(float4& u, const float4& v){
    u.x+=v.x; u.y+=v.y; u.z+=v.z;
    return u;
  }

  inline float4& operator -=(float4& u, const float4& v){
    u.x-=v.x; u.y-=v.y; u.z-=v.z;
    return u;
  }

  inline float4& operator *=(float4& u, const float4& v){
    u.x*=v.x; u.y*=v.y; u.z*=v.z;
    return u;
  }

  inline float4& operator *=(float4& u, const float f){
    u.x*=f; u.y*=f; u.z*=f;
    return u;
  }

  inline float4& operator /=(float4& u, const float4& v){
    u.x/=v.x; u.y/=v.y; u.z/=v.z;
    return u;
  }

  inline float4& operator /=(float4& u, const float f){
    u.x/=f; u.y/=f; u.z/=f;
    return u;
  }

  inline const float sprod(float4& u, const float4& v){
    return u.x*v.x+u.y*v.y+u.z*v.z;
  }

  inline void add(const float4& u, const float4& v, float4& res){
    res.x=u.x+v.x;
    res.y=u.y+v.y;
    res.z=u.z+v.z;
  }

  inline void sub(const float4& u, const float4& v, float4& res){
    res.x=u.x-v.x;
    res.y=u.y-v.y;
    res.z=u.z-v.z;
  }

  inline void lin(const float4& u, const float& b, const float4& v, float4& res){
    res.x=u.x+b*v.x;
    res.y=u.y+b*v.y;
    res.z=u.z+b*v.z;
  }

  inline void lin(const float& a, const float4& u, const float4& v, float4& res){
    res.x=a+u.x+v.x;
    res.y=a*u.y+v.y;
    res.z=a*u.z+v.z;
  }

  inline void lin(const float& a, const float4& u, const float& b, const float4& v, float4& res){
    res.x=a+u.x+b*v.x;
    res.y=a*u.y+b*v.y;
    res.z=a*u.z+b*v.z;
  }



  inline void print(const float4& v, std::ostream& out, const std::string& name=""){
    out<<std::setprecision(3);
    out<<"float4: "<<name<<"=\n"<<v.x<<"\n"<<v.y<<"\n"<<v.z<<" \n"<<std::endl;
  }

#define f4print(v) print(v,std::cout,#v);


  ///   **************************
  ///   *        trafo4          *
  ///   *    implementations     *
  ///   **************************

  trafo4::trafo4(float a, float alpha, float q, float d){
    float cq=cos(q);
    float sq=sin(q);
    float ca=cos(alpha);
    float sa=sin(alpha);

    col[0].x=cq;    col[1].x=-sq;   col[2].x=0.0;  col[3].x=a;
    col[0].y=sq*ca; col[1].y=cq*ca; col[2].y=-sa;  col[3].y=-sa*d;
    col[0].z=sq*sa; col[1].z=cq*sa; col[2].z=ca;   col[3].z=ca*d;
  }

  //!u=T*u
  inline float4& trafo4::apply(float4& u) const
  {
    float a=u.x,b=u.y,c=u.z;
    u.x=col[0].x*a+col[1].x*b+col[2].x*c+col[3].x;
    u.y=col[0].y*a+col[1].y*b+col[2].y*c+col[3].y;
    u.z=col[0].z*a+col[1].z*b+col[2].z*c+col[3].z;
    return u;
  }
  //!Tu=this*u
  inline float4& trafo4::apply(const float4& u, float4& Tu) const
  {
    Tu.x=col[0].x*u.x+col[1].x*u.y+col[2].x*u.z+col[3].x;
    Tu.y=col[0].y*u.x+col[1].y*u.y+col[2].y*u.z+col[3].y;
    Tu.z=col[0].z*u.x+col[1].z*u.y+col[2].z*u.z+col[3].z;
    return Tu;
  }
  //!u=R*u
  inline float4& trafo4::apply_rot(float4& u) const
  {
    float a=u.x,b=u.y,c=u.z;
    u.x=col[0].x*a+col[1].x*b+col[2].x*c;
    u.y=col[0].y*a+col[1].y*b+col[2].y*c;
    u.z=col[0].z*a+col[1].z*b+col[2].z*c;
    return u;
  }
  //!Ru=R*u
  inline float4& trafo4::apply_rot(float4& u, float4& Ru) const
  {
    Ru.x=col[0].x*u.x+col[1].x*u.y+col[2].x*u.z;
    Ru.y=col[0].y*u.x+col[1].y*u.y+col[2].y*u.z;
    Ru.z=col[0].z*u.x+col[1].z*u.y+col[2].z*u.z;
    return Ru;
  }
  //!u=R'*u
  inline float4& trafo4::apply_rot_inv(float4& u) const
  {
    float a=u.x,b=u.y,c=u.z;
    u.x=col[0].x*a+col[0].y*b+col[0].z*c;
    u.y=col[1].x*a+col[1].y*b+col[1].z*c;
    u.z=col[2].x*a+col[2].y*b+col[2].z*c;
    return u;
  }
  //!Rtu=R'*u
  inline float4& trafo4::apply_rot_inv(float4& u, float4& Rtu) const
  {
    Rtu.x=col[0].x*u.x+col[0].y*u.y+col[0].z*u.z;
    Rtu.y=col[1].x*u.x+col[1].y*u.y+col[1].z*u.z;
    Rtu.z=col[2].x*u.x+col[2].y*u.y+col[2].z*u.z;
    return Rtu;
  }

  //!this=this*T
  //inline trafo4& trafo4::apply(const trafo4& T)
  //{
  //
  //}
  //!Tres=this*T
  inline trafo4& trafo4::apply(const trafo4& T, trafo4& Tres) const
  {
    apply(T.col[0],Tres.col[0]);
    apply(T.col[1],Tres.col[1]);
    apply(T.col[2],Tres.col[2]);
    apply(T.col[3],Tres.col[3]);
  }
  //!this=T*this
  inline trafo4& trafo4::lapply(trafo4& T)
  {
    T.apply(col[0]);
    T.apply(col[1]);
    T.apply(col[2]);
    T.apply(col[3]);
  }
  //!Tres=T*this (not needed)
  inline trafo4& trafo4::lapply(const trafo4& T, trafo4& Tres) const
  {
    T.apply(col[0],Tres.col[0]);
    T.apply(col[1],Tres.col[1]);
    T.apply(col[2],Tres.col[2]);
    T.apply(col[3],Tres.col[3]);
  }


  inline void trafo4::print(std::ostream& out, const std::string& name){
    out<<"trafo4: "<<name<<"="<<std::endl;
    out<<std::setprecision(3);
    out<<col[0].x<<"\t"<<col[1].x<<"\t"<<col[2].x<<"\t"<<col[3].x<<std::endl;
    out<<col[0].y<<"\t"<<col[1].y<<"\t"<<col[2].y<<"\t"<<col[3].y<<std::endl;
    out<<col[0].z<<"\t"<<col[1].z<<"\t"<<col[2].z<<"\t"<<col[3].z<<std::endl<<std::endl;
  }

#define t4print(T) T.print(std::cout,#T);


}//namespace geo4

#endif // GEO4_H
