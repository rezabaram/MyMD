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

double aG[]={0.0, 0, -10.0};
vec G(aG);

void define_parameters()
{
	config.add_param<vec>("Gravity", G);
	config.add_param<double>("outDt", 0.02);
	config.add_param<double>("stiffness1", 5000000.0); 
	config.add_param<double>("stiffness2", 1000000.0); 
	config.add_param<double>("density", 10000.0); 
	config.add_param<double>("particleSize", 0.05); 
	config.add_param<double>("timeStep", 0.0001); 
	config.add_param<double>("maxTime", 1.0); 
	config.add_param<int>("nParticle", 10); 
}

template <class T>
string stringify(T x, int width=15, const char ch=' ')
 {
   std::ostringstream o;
   if (!(o << setw(width)<<setfill(ch)<<x))
     cerr<<"Bad coversion to string"<<endl;
   return o.str();
 }


using namespace math;
typedef matrix<double> Matrix;
void quaternionToMatrix(const Quaternion &q, Matrix M){
	double Nq = q.abs2();
	static double s;
	if( Nq > 0.0) s = 2.0/Nq; else s = 0.0;
	double X = q.v(0)*s,   Y = q.v(1)*s,  Z = q.v(2)*s;
	double wX = q.u*X, wY = q.u*Y, wZ = q.u*Z;
	double xX = q.v(0)*X, xY = q.v(0)*Y, xZ = q.v(0)*Z;
	double yY = q.v(1)*Y, yZ = q.v(1)*Z, zZ = q.v(2)*Z;

	M(0,0)=1.0-(yY+zZ); M(0,0)=xY-wZ ;       M(0,0)= xZ+wY;
	M(0,0)= xY+wZ;      M(0,0)=1.0-(xX+zZ);  M(0,0)=yZ-wX;
	M(0,0)=xZ-wY;       M(0,0)= yZ+wX;       M(0,0)=1.0-(xX+yY);

return;
}

#endif /* COMMON_H */
