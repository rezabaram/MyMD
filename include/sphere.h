// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef SPHERE_H
#define SPHERE_H 
#include"geombase.h"

class CSphere: public GeomObjectBase
	{
	public:	
	CSphere (const vec &v, double r):GeomObjectBase(v,tsphere){radius=r;}
	CSphere ():GeomObjectBase(vec(0,0,0),tsphere){radius=1;}
	virtual ~CSphere(){};

	void rotate(const vec&, double alpha){};
	void shift(const vec& v){Xc+=v;};
	void scale(double scale){radius*=scale; Xc*=scale;};

	void printRaster3D(std::ostream &out)const{
		print(out);
		}

	void print(std::ostream &out)const{
		out<< identifier<< "   ";
		out<< Xc<<"  "<<radius;
		}

	void parse(std::istream &in){
			in>>identifier;
			in>> Xc >>radius;
			}

	double vol()const{ 
			return 4.0/3.0*M_PI*radius*radius*radius;}
	double I(vec n){return 2.0/5.0*radius*radius;}
	private:
	};


#endif /* SPHERE_H */
