#ifndef RAY_H
#define RAY_H 

class CRay 
	{
	public:
	CRay(const vec &_x1, const vec &_x2):x1(_x1),x2(_x2), n((_x2-_x1).normalize()){};

	vec operator()(double t){
		return x1+t*n;
		}

	void print(std::ostream &out)const{
		out<< 5 << "   ";
		out<< x1<< "  "<<0.005<<"  ";
		out<< x2<< "  "<<0.005<<endl;
		};

	vec x1, x2, n;
 	private:
	
	};

#endif /* RAY_H */
