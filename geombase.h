#ifndef GEOMBASE_H
#include<ostream>
#include"vec3d.h"

typedef enum {tsphere, tplane, tbox, tcomposite, tellipsoid} GType;

class GeomObjectBase
	{
	public:
	GeomObjectBase(const vec &v, GType t):Xc(v),Xc0(v),type(t){identifier=1;};
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
	virtual void parse(std::istream &in)=0;

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

	GType type;
	double radius;
	
	vec Xc, Xc0; //center 
	int identifier;
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
#define GEOMBASE_H 
#endif /* GEOMBASE_H */
