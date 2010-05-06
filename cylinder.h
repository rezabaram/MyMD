#ifndef CYLINDER_H
#define CYLINDER_H 

#include"geombase.h"

using namespace math;

template<>
class GeomObject<tcylinder>: public GeomObjectBase{
	public:	
	GeomObject(const vec &v1,const vec &v2, double _r):GeomObjectBase((v1+v2)/2.0,tcylinder), X1(v1), X2(v2), r(_r){
		identifier=5;
		}

	void shift(const vec&){};
	void rotate(const vec& n, double alpha){ };
	void rotateTo(const Quaternion &q){};
	double vol(){return 1;};
	double I(vec n){return 1;};

	void scale(double){};
	void print(std::ostream &out)const{
		out<< identifier<< "   ";
		out<< X1<< "  "<<r<<"  ";
		out<< X2<< "  "<<r;
		};
	void parse(std::istream &in){};

	vec X1, X2;
	double r;

	};

typedef GeomObject<tcylinder> CCylinder;
#endif
