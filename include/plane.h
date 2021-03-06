// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef PLANE_H
#define PLANE_H 
#include"geombase.h"

class CPlane : public GeomObjectBase
	{
	public:
		CPlane(const vec &x0, const vec &n0):GeomObjectBase(x0, tplane), n(n0), solid(true){identifier=6; n.normalize();};
		virtual ~CPlane(){};
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
		double operator()(const vec &point)const{
			return (point-Xc)*n;
			}

		vec normal_to_point(const vec & p, double shift=0)const{
			return -normal_from_point(p,shift);
			}

		vec normal_from_point(const vec & p, double shift=0)const{
			return ((Xc-p)*n-shift)*n;
			}
	double vol()const{return 0;}
	double I(vec n){
		ERROR(true, "Not implemented."); //FIXME
		return 0;}
	
	vec n;//normal
	vec vec_to_shadow;
	bool solid;
 	private:
	};


class HomPlane
	{
	public:
	HomPlane(){}
	HomVec x, n;
 	private:
	};

#endif /* PLANE_H */
