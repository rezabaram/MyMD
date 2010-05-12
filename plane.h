#ifndef PLANE_H
#define PLANE_H 
#include"geombase.h"

template<>
class GeomObject<tplane> : public GeomObjectBase
	{
	public:
		GeomObject(const vec &x0, const vec &n0):GeomObjectBase(x0, tplane), n(n0){identifier=6; n.normalize();};
		virtual ~GeomObject(){};
		void shift(const vec & v){Xc+=v;};
		void rotate(const vec &_n, double alpha){};//FIXME
		void scale(double s){Xc*=s;};//FIXME if necessary
		void print(std::ostream &out)const{
			out<< identifier<< "   ";
			out<< Xc<< "  "<<n;
			}


		void rasterprint(std::ostream &out){
			vec n2(drand48(), drand48(), drand48());
			n2.normalize();
			n2=cross(n,n2);
			vec n3=cross(n,n2);
			out<< identifier<< "   "<<Xc<<"  "<<Xc+n2<<"  "<<Xc+n3<<endl;
			}

		void parse(std::istream &in){
			in>>identifier;
			in>>Xc>>n;
			}

		vec normal_to_point(const vec & p, double shift=0)const{
			return ((p-Xc)*n-shift)*n;
			}

		vec normal_from_point(const vec & p, double shift=0)const{
			return ((Xc-p)*n-shift)*n;
			}
	double vol(){return 0;}
	double I(vec n){
		ERROR(true, "Not implemented."); //FIXME
		return 0;}
	
	vec n;//normal
 	private:
	};

typedef GeomObject<tplane> CPlane;

#endif /* PLANE_H */
