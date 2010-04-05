#ifndef COMMON_H
#define COMMON_H 

#define ERROR(x)  std::cerr<<"Error: In file " __FILE__<<" line "<<__LINE__<<":  "<<x<<std::endl;
//#include"matrix.h"
#include<assert.h>
#include<memory>
#include<limits>
#include<iomanip>
#include<sstream>
#include<ostream>
#include<vector>
#include"vec3d.h"
#include"quaternion.h"
#include"matrix.h"
#include"log.h"

#include"config.h"
extern CConfig &config; // don't forget "&" or you get a vicious bug, which took me one day to find

typedef vec3d<double> vec;

double aG[]={0.0, 0, -2.0};
vec G(aG);

void define_parameters()
{
	config.add_param<vec>("Gravity", G);
	config.add_param<double>("outDt", 0.02);
	config.add_param<double>("stiffness", 5.0e+07); 
	config.add_param<double>("damping", 0.00005); 
	config.add_param<double>("density", 10000.0); 
	config.add_param<double>("particleSize", 0.05); 
	config.add_param<double>("timeStep", 0.00001); 
	config.add_param<double>("maxTime", 4.0); 
	config.add_param<int>("nParticle", 5); 
}

template <class T>
string stringify(T x, int width=15, const char ch=' ')
 {
   std::ostringstream o;
   if (!(o << setw(width)<<setfill(ch)<<x))
     cerr<<"Bad coversion to string"<<endl;
   return o.str();
 }

template<typename T>
inline
T tmax(const T &a, const T &b){
	if(a>b)return a;
	return b;
	}

using namespace math;
typedef matrix<double> Matrix;
void quaternionToMatrix(const Quaternion &q, Matrix &M){//got from wikipedia
	double Nq = q.abs2();
	static double s;
	if( Nq > 0.0) s = 2.0/Nq; else s = 0.0;
	double X = q.v(0)*s,   Y = q.v(1)*s,  Z = q.v(2)*s;

	double wX = q.u*X,    wY = q.u*Y,    wZ = q.u*Z;
	double xX = q.v(0)*X, xY = q.v(0)*Y, xZ = q.v(0)*Z;
	double yY = q.v(1)*Y, yZ = q.v(1)*Z, zZ = q.v(2)*Z;

	M(0,0)=1.0-(yY+zZ); M(1,0)=xY-wZ ;       M(2,0)= xZ+wY;
	M(0,1)= xY+wZ;      M(1,1)=1.0-(xX+zZ);  M(2,1)=yZ-wX;
	M(0,2)=xZ-wY;       M(1,2)= yZ+wX;       M(2,2)=1.0-(xX+yY);
	

return;
}

inline
vec operator *(const vec &v, const Matrix &M){

	return vec(
		v(0)*M(0,0)+v(1)*M(1,0)+v(2)*M(2,0),
		v(1)*M(0,1)+v(1)*M(1,1)+v(2)*M(2,1),
		v(2)*M(0,2)+v(1)*M(1,2)+v(2)*M(2,2)
		);
	}

inline
vec operator *(const Matrix &M, const vec &v){

	return vec(
		v(0)*M(0,0)+v(1)*M(0,1)+v(2)*M(0,2),
		v(1)*M(1,0)+v(1)*M(1,1)+v(2)*M(1,2),
		v(2)*M(2,0)+v(1)*M(2,1)+v(2)*M(2,2)
		);
	}
#endif /* COMMON_H */
