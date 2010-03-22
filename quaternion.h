#ifndef QUATERNION_H
#define QUATERNION_H 

#include <stdio.h>
#include <math.h>
#include <fstream>
#include "log.h"
#include <string>
#include "common.h"

class Quaternion{
protected:
public:
double u;
vec3d<double> v;

Quaternion();
Quaternion(double u,double x,double y,double z);
Quaternion(const vec3d<double> &_v, double a);
Quaternion(double a, const vec3d<double> &_v);

void set(double u,double x,double y,double z);
void set(const vec3d<double> &axis, double angle);
void set(double angle, const vec3d<double> &axis);
void setRotation(const vec3d<double> &_V, double _a);
//void setRotation(double _a, const vec3d<double> &_V){setRotation(_V, a);};

double abs();
vec3d<double> getAxis(void)const{return v;};


vec3d<double> rotate(const vec3d<double> &v)const;
Quaternion rotate(const Quaternion &q)const;

void rotateMe(const vec3d<double> &axis, double angle);

/// Convert a vector in body coordinates to a vector in world coordinates (with me as base)
vec3d<double> toWorld(const vec3d<double> &v)const;
/// Convert a vector in world coordinates to a vector in body coordinates (with me as base)
vec3d<double> toBody(const vec3d<double> &v);

Quaternion normalize(){ (*this)/=abs(); return *this; }
Quaternion normalized()const{ return Quaternion(*this).normalize(); }

Quaternion& operator+=(const Quaternion &);
Quaternion& operator-=(const Quaternion &);
Quaternion& operator*=(double);
Quaternion& operator/=(double);
Quaternion& operator=(const vec3d<double> &);

friend Quaternion operator~(const Quaternion &);
friend Quaternion operator+(const Quaternion &,const Quaternion &);
friend Quaternion operator-(const Quaternion &,const Quaternion &);
friend Quaternion operator*(const Quaternion &,double);
friend Quaternion operator/(const Quaternion &,double);

friend Quaternion operator*(const Quaternion &,const Quaternion &);
friend Quaternion operator*(const Quaternion &,const vec3d<double> &);
friend Quaternion operator*(const vec3d<double> &,const Quaternion &);
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
setRotation(v, u);
}

inline void Quaternion::set(double _u,double _x,double _y,double _z){
u=_u;  
v(0)=_x; 
v(1)=_y;
v(2)=_z;
}


inline void Quaternion::set(double _a, const vec3d<double> &_V){
u=_a;  
v=_V;
}

inline Quaternion::Quaternion(  const vec3d<double> &_V, double _a){
set(_a,_V);
}

inline Quaternion::Quaternion( double _a, const vec3d<double> &_V){
set(_a,_V);
}

inline void Quaternion::setRotation(const vec3d<double> &_V, double _a){
u=cos(_a/2);
v=_V.normalized()*sin(_a/2.0);
}

std::ostream &operator<<(std::ostream & out, const Quaternion &q){
	out<<q.u<<"  "<<q.v;
return out;
}

inline double Quaternion::abs(){
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

inline Quaternion operator*(const Quaternion &q,const vec3d<double> &v){
return Quaternion(-q.v*v, q.u*v+cross(q.v,v));
}

inline Quaternion operator*(const vec3d<double> &v,const Quaternion &q){
return Quaternion(-q.v*v, q.u*v+cross(v,q.v));
}

inline Quaternion &Quaternion::operator=(const vec3d<double> &_v){
u=0; 
v=_v;
return *this;
}

inline vec3d<double> Quaternion::rotate(const vec3d<double> &v)const{
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

inline void Quaternion::rotateMe(const vec3d<double> &axis, double angle){
log <<"Check this!"<<std::endl; // FIXME
*this=Quaternion(axis,angle).rotate(*this);
}

/**
* Just a rotation
*/
inline vec3d<double> Quaternion::toWorld(const vec3d<double> &_v)const{
return ((*this)*_v*(~(*this))).getAxis();
}

/**
* Just a rotation vs the conjugate
*/
inline vec3d<double> Quaternion::toBody(const vec3d<double> &_v){
return ((~(*this))*_v*(*this)).v;
}
#endif /* QUATERNION_H */
