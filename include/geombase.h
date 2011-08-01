#ifndef GEOMBASE_H
#define GEOMBASE_H 
#include<ostream>
#include"vec.h"
#include"quaternion.h"

typedef enum {tsphere, tplane, tbox, tcomposite, tellipsoid, tcylinder} GType;

class CPlane;//forward declaration
class GeomObjectBase
	{
	public:
	GeomObjectBase(const vec &v, GType t, const Quaternion &_q=Quaternion(1,0,0,0)):Xc(v),Xc0(v),type(t), q(_q){identifier=1;};
	virtual ~GeomObjectBase(){};
	virtual void shift(const vec&)=0;//{WARNING("This function should not be called")};
	virtual void rotate(const vec& n, double alpha){//maybe overriden by derived class
		//q.rotateMe(n, alpha);
		//Xc0=q.rotate(Xc0);
		};
	virtual void rotateTo(const Quaternion &q){}
	virtual double vol()const=0;//{WARNING("This function should not be called");}
	virtual double I(vec n)=0;//{WARNING("This function should not be called");}

	virtual void scale(double)=0;//{WARNING("This function should not be called");}
	virtual void print(std::ostream &out)const=0;//{WARNING("This function should not be called");}
	virtual void printRaster3D(std::ostream &out)const{WARNING("This function should not be called");}
	virtual void print_in_euler(std::ostream &out)const{WARNING("This function should not be called");}
	virtual void parse(std::istream &in)=0;
	virtual void fixToBody(const HomVec &point){WARNING("This function should not be called");}
	virtual double operator()(const vec &point)const{return 1e+100;}
	virtual GeomObjectBase *clone()const{WARNING("This function should not be called");return NULL;}
	virtual bool doesHit(const CPlane &plane)const {WARNING("This function should not be called");return false;}

	virtual void moveto(const vec& v){//derived classes can override this. especially composite particles should!
		Xc=v;
		};

	double top()const{
		return Xc(2)+radius;
		}
	const vec& getpos()const{
		return Xc;
		};
	double distance(const GeomObjectBase *p)const{
		return (Xc-p->Xc).abs();
		};
	const vec displacement(const GeomObjectBase *p)const{
		return (Xc-p->Xc);
		};
	virtual const void print_coord_sys(ostream &out){
			ERROR(1,"Function not implemented");
			};
	vec euler()const{
		vec angles;
		//from wiki 
		//phi
		angles(0)=atan2 (2.0*(q(0)*q(1) + q(2)*q(3)), 1 - 2*(q(1)*q(1) + q(2)*q(2))) ;
		//theta
		angles(1)=asin(2.0*(q(0)*q(2) - q(3)*q(1)));
		//phi
		angles(2)=atan2(2.0*(q(0)*q(3) + q(2)*q(1)), 1 - 2*(q(2)*q(2) + q(3)*q(3))) ;
		return angles;
		}

	double radius;
	vec Xc, Xc0; //center 
	HomVec P, P0; //test point
	int identifier;
	GType type;
	Quaternion q;
	protected:
 	private:
	};
/*
template <GType T>
class GeomObject:public GeomObjectBase{ //empty class. all the shapes are created through template specialization
		private:
		GeomObject();
		virtual ~GeomObject(){};
		GeomObject(GeomObject const &);
		};
*/

#endif /* GEOMBASE_H */
