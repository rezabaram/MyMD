// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef RAY_H
#define RAY_H 
#include"vec.h"

template<typename T>
class CRay 
	{
	public:
	CRay(const T &_x1, const T &_x2):x1(_x1),x2(_x2), n((_x2-_x1)){n.normalize();};

	T operator()(double t)const{
		return x1+t*n;
		}
	const T &get_n()const{return n;}
	T &get_n(){return n;}

	void print(std::ostream &out)const{
		out<< 5 << "   ";
		out<< x1<< "  "<<0.005<<"  ";
		out<< x2<< "  "<<0.005<<endl;
		};

	T x1, x2;
	T n;
 	private:
	
	};

#endif /* RAY_H */
