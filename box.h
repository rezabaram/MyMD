#ifndef BOX_H
#define BOX_H 
#include<limits>
#include"geombase.h"

template<>
class GeomObject<tbox>: public GeomObjectBase
	{
	public:
	virtual ~GeomObject(){};
	GeomObject(vec corner=vec(std::numeric_limits<double>::max()), vec _L=vec(0.0)):
		GeomObjectBase(corner+_L/0.5, tbox), corner(corner), L(_L), nFaces(5),
		u0(vec(1.0,0.0,0.0)), u1(vec(0.0,1.0,0.0)), u2(vec(0.0,0.0,1.0))
		 {
		identifier=6;
		face =  new GeomObject<tplane> * [nFaces];
		face[0]=new GeomObject<tplane> (corner,u0);
		face[1]=new GeomObject<tplane> (corner,u1);
		face[2]=new GeomObject<tplane> (corner,u2);

		face[3]=new GeomObject<tplane> (corner+L,-u0);
		face[4]=new GeomObject<tplane> (corner+L,-u1);
//		face[5]=new GeomObject<tplane> (corner+L,-u2);
		};

	void shift(const vec &x){corner+=x;}
	void rotate(const vec &n, double alpha){}//FIXME 
	void scale(double scale){L*=scale;};
	void print(std::ostream &out)const{//FIXME this is ad-hoc
		out<< identifier<< "   "<<corner<<"  "<<corner+vec(1,0,0)<<"  "<<corner+vec(1,1,0)<<endl;
		out<< identifier<< "   "<<corner<<"  "<<corner+vec(0,0,1)<<"  "<<corner+vec(0,1,1)<<endl;
		out<< identifier<< "   "<<corner+L<<"  "<<corner+L+vec(0,0,1)<<"  "<<corner+vec(0,1,1)<<endl;
		out<< identifier<< "   "<<corner+L<<"  "<<corner+L+vec(0,0,1)<<"  "<<corner+vec(1,0,1)<<endl;
		out<< identifier<< "   "<<corner+L<<"  "<<corner+L+vec(1,0,0)<<"  "<<corner+L+vec(1,1,0)<<endl;
		}

	double vol(){return L(0)*L(1)*L(2);}
	double I(vec n){
		ERROR(true, "Not implemented."); //FIXME
		return 0;}

	void parse(std::istream &in){
		in>>identifier;
		in>>corner>>L;
		}

	void scale(double s0, double s1, double s2){L(0)*=s0; L(1)*=s1; L(2)*=s2;}

	vec top()const{return corner+L;}
	vec corner, L;
	const size_t nFaces;
	GeomObject<tplane> ** face;
 	private:
	GeomObject();
	vec u0, u1, u2;
	};


#endif /* BOX_H */
