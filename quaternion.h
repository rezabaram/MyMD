#ifndef QUATERNION_H
#define QUATERNION_H 

#include <stdio.h>
#include <math.h>
#include <fstream>
#include "log.h"
#include <string>
#include "common.h"

class Quaternion{
Quaternion();
protected:
public:
double u;
vec v;

Quaternion(double u,double x,double y,double z);
Quaternion(const vec &_v, double a);
Quaternion(double a, const vec &_v);

void set(double u,double x,double y,double z);
void set(const vec &axis, double angle);
void set(double angle, const vec &axis);
void setRotation(const vec &_V, double _a);
//void setRotation(double _a, const vec &_V){setRotation(_V, a);};

double abs()const;
double abs2()const;
vec getAxis(void)const{return v;};


vec rotate(const vec &v)const;
Quaternion rotate(const Quaternion &q)const;

void rotateMe(const vec &axis, double angle);

/// Convert a vector in body coordinates to a vector in world coordinates (with me as base)
vec toWorld(const vec &v)const;
/// Convert a vector in world coordinates to a vector in body coordinates (with me as base)
vec toBody(const vec &v)const;

Quaternion normalize(){ (*this)/=abs(); return *this; }
Quaternion normalized()const{ return Quaternion(*this).normalize(); }

Quaternion& operator+=(const Quaternion &);
Quaternion& operator-=(const Quaternion &);
Quaternion& operator*=(double);
Quaternion& operator/=(double);
Quaternion& operator=(const vec &);

friend Quaternion operator~(const Quaternion &);
friend Quaternion operator+(const Quaternion &,const Quaternion &);
friend Quaternion operator-(const Quaternion &,const Quaternion &);
friend Quaternion operator*(const Quaternion &,double);
friend Quaternion operator/(const Quaternion &,double);

friend Quaternion operator*(const Quaternion &,const Quaternion &);
friend Quaternion operator*(const Quaternion &,const vec &);
friend Quaternion operator*(const vec &,const Quaternion &);
friend std::ostream &operator<<(std::ostream & out, const Quaternion &q);
};


inline Quaternion::Quaternion(){
u=1;
}

inline Quaternion::Quaternion(double _u,double _x,double _y,double _z){
u=_u;  
v(0)=_x; 
v(1)=_y;
v(2)=_z;
}

inline void Quaternion::set(double _u,double _x,double _y,double _z){
u=_u;  
v(0)=_x; 
v(1)=_y;
v(2)=_z;
}


inline void Quaternion::set(double _a, const vec &_V){
u=_a;  
v=_V;
}

inline Quaternion::Quaternion(  const vec &_V, double _a){
set(_a,_V);
}

inline Quaternion::Quaternion( double _a, const vec &_V){
set(_a,_V);
}

inline void Quaternion::setRotation(const vec &_V, double _a){
u=cos(_a/2);
v=_V.normalized()*sin(_a/2.0);
}

std::ostream &operator<<(std::ostream & out, const Quaternion &q){
	out<<q.u<<"  "<<q.v;
return out;
}

inline double Quaternion::abs2()const{
return u*u+v.abs2();
}

inline double Quaternion::abs()const{
return sqrt(u*u+v.abs2());
}

inline Quaternion operator+(const Quaternion &a,const Quaternion &b){
return Quaternion(a.u+b.u, a.v+b.v);
}

inline Quaternion& Quaternion::operator+=(const Quaternion &a){
u+=a.u;
v+=a.v;
return *this;
}

inline Quaternion operator-(const Quaternion &a,const Quaternion &b){
return Quaternion(a.u-b.u, a.v-b.v);
}

inline Quaternion& Quaternion::operator-=(const Quaternion &a){
u-=a.u;
v-=a.v;
return *this;
}

inline Quaternion operator*(const Quaternion &a, double s){
return Quaternion(a.u*s, a.v*s);
}

inline Quaternion& Quaternion::operator*=(double s){
u*=s;
v*=s;
return *this;
}

inline Quaternion operator/(const Quaternion &a, double s){
return Quaternion(a.u/s, a.v/s);
}

inline Quaternion& Quaternion::operator/=(double s){
u/=s;
v/=s;
return *this;
}

inline Quaternion operator~(const Quaternion &a){
return Quaternion(a.u, -a.v);
}

inline Quaternion operator*(const Quaternion &a, const Quaternion &b){
//return Quaternion(a.u*b.u - a.x*b.x - a.y*b.y - a.z*b.z,
//a.u*b.x + a.x*b.u + a.y*b.z - a.z*b.y,
//a.u*b.y + a.y*b.u + a.z*b.x - a.x*b.z,
//a.u*b.z + a.z*b.u + a.x*b.y - a.y*b.x);
return Quaternion(a.u*b.u-a.v*b.v, a.u*b.v+b.u*a.v+cross(a.v,b.v));
}

inline Quaternion operator*(const Quaternion &q,const vec &v){
return Quaternion(-q.v*v, q.u*v+cross(q.v,v));
}

inline Quaternion operator*(const vec &v,const Quaternion &q){
return Quaternion(-q.v*v, q.u*v+cross(v,q.v));
}

inline Quaternion &Quaternion::operator=(const vec &_v){
u=0; 
v=_v;
return *this;
}

inline vec Quaternion::rotate(const vec &v)const{
Quaternion q=(*this)*v*(~(*this));
return q.v;
}

/**
* This method rotates the given quaternion over ME.
*/
inline Quaternion Quaternion::rotate(const Quaternion &q)const{
Quaternion rq=(*this)*q*(~(*this));
return rq;
}

inline void Quaternion::rotateMe(const vec &axis, double angle){
log <<"Check this!"<<std::endl; // FIXME
*this=Quaternion(axis,angle).rotate(*this);
}

/**
* Just a rotation
*/
inline vec Quaternion::toWorld(const vec &_v)const{
return ((*this)*_v*(~(*this))).getAxis();
}

/**
* Just a rotation vs the conjugate
*/
inline vec Quaternion::toBody(const vec &_v)const{
return ((~(*this))*_v*(*this)).v;
}
#endif /* QUATERNION_H */
