#ifndef RAY_H
#define RAY_H 
#include"vec.h"

template<typename T=vec>
class CRay 
	{
	public:
	CRay(const T &_x1, const T &_x2):x1(_x1),x2(_x2), n((_x2-_x1).normalize()){};

	T operator()(double t)const{
		return x1+t*n;
		}

	void print(std::ostream &out)const{
		out<< 5 << "   ";
		out<< x1<< "  "<<0.005<<"  ";
		out<< x2<< "  "<<0.005<<endl;
		};

	T x1, x2, n;
 	private:
	
	};

#endif /* RAY_H */
