#ifndef GEOMBASE_H
#define GEOMBASE_H 
#include<ostream>
#include"vec.h"
#include"quaternion.h"

typedef enum {tsphere, tplane, tbox, tcomposite, tellipsoid, tcylinder} GType;

class GeomObjectBase
	{
	public:
	GeomObjectBase(const vec &v, GType t, const Quaternion &_q=Quaternion(1,0,0,0)):Xc(v),Xc0(v),type(t), q(_q){identifier=1;};
	virtual ~GeomObjectBase(){};
	virtual void shift(const vec&)=0;
	virtual void rotate(const vec& n, double alpha){//maybe overriden by derived class
		//q.rotateMe(n, alpha);
		//Xc0=q.rotate(Xc0);
		};
	virtual void rotateTo(const Quaternion &q){};
	virtual double vol()=0;
	virtual double I(vec n)=0;

	virtual void scale(double)=0;
	virtual void print(std::ostream &out)const=0;
	virtual void printRaster3D(std::ostream &out)const{};
	virtual void parse(std::istream &in)=0;
	virtual void fixToBody(const HomVec &point){};
	virtual double operator()(const vec &point)const{return 1e+100;};

	virtual void moveto(const vec& v){//derived classes can override this. especially composite particles should!
		Xc=v;
		};

	const vec& getpos()const{
		return Xc;
		};
	double distance(const GeomObjectBase *p)const{
		return (Xc-p->Xc).abs();
		};
	const vec displacement(const GeomObjectBase *p)const{
		return (Xc-p->Xc);
		};

	double radius;
	vec Xc, Xc0; //center 
	HomVec P, P0; //test point
	int identifier;
	GType type;
	Quaternion q;
	protected:
 	private:
	};

template <GType T>
class GeomObject:public GeomObjectBase{//empty class. all the shapes are created through template specialization
		private:
		GeomObject();
		virtual ~GeomObject(){};
		GeomObject(GeomObject const &);
		};
#endif /* GEOMBASE_H */
