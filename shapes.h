#ifndef SHAPES_H
#define SHAPES_H 
#include"common.h"
#include"quaternion.h"
using std::cerr; 
using std::endl;
using namespace std;


typedef enum {tsphere, tplane, tbox, tcomposite} GType;


class GeomObjectBase
	{
	public:
	GeomObjectBase(vec const & v, GType t):Xc(v),Xc0(v),type(t), q(0.0, 1.0, 0.0, 0.0){identifier=1;};
	virtual ~GeomObjectBase(){};
	virtual void shift(const vec&)=0;
	virtual void rotate(const vec& n, double alpha){//maybe overriden by derived class
		//q.rotateMe(n, alpha);
		//Xc0=q.rotate(Xc0);
		};
	virtual void rotateTo(const Quaternion &q){};

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
	Quaternion q;//orientation
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

template<>
class GeomObject<tplane> : public GeomObjectBase
	{
	public:
		GeomObject(const vec &x0, const vec &n0):GeomObjectBase(x0, tplane), n(n0){n.normalize();};
		virtual ~GeomObject(){};
		void shift(const vec & v){Xc+=v;};
		void rotate(const vec &_n, double alpha){};//FIXME
		void scale(double s){Xc*=s;};//FIXME if necessary
		void print(std::ostream &out)const{
			out<< identifier<< "   ";
			out<< Xc<< "  "<<n;
			}
		void parse(std::istream &in){
			in>>identifier;
			in>>Xc>>n;
			}

		vec normal_to_point(const vec & p, double shift){
			return ((Xc-p)*n-shift)*n;
			}
	
	vec n;//normal
 	private:
	};

template<>
class GeomObject<tbox>: public GeomObjectBase
	{
	public:
	virtual ~GeomObject(){};
	GeomObject(vec corner=vec(std::numeric_limits<double>::max()), vec _L=vec(0.0)):
		GeomObjectBase(corner+_L/0.5, tbox), corner(corner), L(_L),
		u0(vec(0,1)), u1(vec(1,1.0)), u2(vec(2,1.0))
		 {
		face[0]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner,-u0));
		face[1]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner,-u1));
		face[2]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner,-u2));

		face[3]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner+L,u0));
		face[4]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner+L,u1));
		face[5]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner+L,u2));
		};

	void shift(const vec &x){corner+=x;}
	void rotate(const vec &n, double alpha){}//FIXME 
	void scale(double scale){L*=scale;};
	void print(std::ostream &out)const{
		out<< identifier<< "   ";
		out<<corner<<"  "<<L;
		}

	void parse(std::istream &in){
			in>>identifier;
			in>>corner>>L;
			}

	void scale(double s0, double s1, double s2){L(0)*=s0; L(1)*=s1; L(2)*=s2;}



	vec top()const{return corner+L;}
	vec corner, L;
	auto_ptr<GeomObject<tplane> > face[6];
 	private:
	GeomObject();
	vec u0, u1, u2;
	};

template<>
class GeomObject<tsphere>: public GeomObjectBase
	{
	public:	
	GeomObject<tsphere> (const vec &v, double r):GeomObjectBase(v,tsphere),radius(r){}
	virtual ~GeomObject(){};

	void rotate(const vec&, double alpha){};
	void shift(const vec& v){Xc+=v;};
	void scale(double scale){radius*=scale; Xc*=scale;};

	void print(std::ostream &out)const{
		out<< identifier<< "   ";
		out<< Xc<<"  "<<radius;
		}

	void parse(std::istream &in){
			in>>identifier;
			in>> Xc >>radius;
			}

	double radius;
	private:
	GeomObject<tsphere>();
	};

typedef GeomObject<tsphere> CSphere;
template<>
class GeomObject<tcomposite>: public GeomObjectBase{
	public:	
		
	GeomObject<tcomposite> (const vec &v, double r):GeomObjectBase(v,tcomposite), shell(v,2.5*r){
		N++;
		GeomObjectBase *s1=NULL;
		s1=new CSphere(vec(0,-r/4), r);
		GeomObjectBase *s2=NULL;
		s2=new CSphere(vec(0,r/4), r);
		s2->identifier=2;

		if(s1==NULL || s2==NULL){ERROR("error in memory allocation"); exit(1);}
		elems.push_back(s1);
		elems.push_back(s2);
		//s2=new CSphere(vec(2,+r/2.), r);
		//elems.push_back(s2);
		//s2=new CSphere(vec(2,-r/2.), r);
		//elems.push_back(s2);
		//s2=new CSphere(vec(1,-r/2.), r);
		//elems.push_back(s2);
		//s2=new CSphere(vec(1,+r/2.), r);
		//elems.push_back(s2);

		Xc=v;
		for(int i=0; i<elems.size(); i++){
			elems.at(i)->Xc=Xc+elems.at(i)->Xc0;
			}
		}

	~GeomObject<tcomposite>(){
		for(int i=0; i<elems.size(); i++){
			delete elems.at(i);
			elems.at(i)=NULL;
			}
		}

	void moveto(const vec &v){
		vec dx = v-Xc;
		for(int i=0; i<elems.size(); i++){
			elems.at(i)->Xc+=dx;
			}
		Xc+=dx;
		shell.Xc+=dx;
		}

	void rotateTo(const Quaternion &q){
		for(int i=0; i<elems.size(); i++){
			elems.at(i)->Xc=Xc+q.rotate(elems.at(i)->Xc0);
			}
		};

	void rotate(const vec& n , double alpha){
		q.setRotation(n, alpha);
		for(int i=0; i<elems.size(); i++){
			elems.at(i)->Xc0=q.rotate(elems.at(i)->Xc0);
			}
		};

	void shift(const vec& v){Xc+=v;};
	void scale(double scale){};//FIXME

	void print(std::ostream &out)const{
		for(int i=0; i<elems.size()-1; i++){
			elems.at(i)->print(out);
			out<< endl;
			}
			elems.back()->print(out);
		}
	
	void parse(std::istream &in)const{//FIXME
			ERROR("not implemented");
			//in>>identifier;
			}

	vector<GeomObjectBase *> elems;
	CSphere shell;
	private:
	static int N;
	GeomObject<tcomposite> (const GeomObject<tcomposite> & p);//not allow copies
	GeomObject<tcomposite> ();
	};

int GeomObject<tcomposite>::N=0;
class COverlapping{
	public:
	COverlapping();
	COverlapping(const vec &_x, const vec &_dx ):x(_x), dx(_dx){}
	static void overlaps(vector<COverlapping> &ovs, const GeomObjectBase *p1, const GeomObjectBase *p2){
		if(p1->type==tsphere && p2->type==tsphere)
			overlaps(ovs, static_cast<const GeomObject<tsphere> *>(p1), static_cast<const GeomObject<tsphere> *>(p2));
		else if(p1->type==tsphere && p2->type==tbox)
			overlaps(ovs, static_cast<const GeomObject<tsphere> *>(p1), static_cast<const GeomObject<tbox> *>(p2));
		else if(p1->type==tcomposite && p2->type==tcomposite)
			overlaps(ovs, static_cast<const GeomObject<tcomposite> *>(p1), static_cast<const GeomObject<tcomposite> *>(p2));
		else if(p1->type==tcomposite && p2->type==tbox)//FIXME
			overlaps(ovs, static_cast<const GeomObject<tcomposite> *>(p1), static_cast<const GeomObject<tbox> *>(p2));
		else ERROR("Not Implemented");
		};

	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tbox> *b);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tsphere>  * p2);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tcomposite>  * p2);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tbox>  *b);

	vec x, dx;
	};

void COverlapping::overlaps(vector<COverlapping > &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tsphere>  * p2){
	static vec v;
	static double d, dd;
	v=p2->displacement(p1);//Xc2-Xc1
	d=v.abs();
	dd=p2->radius+p1->radius-d;//FIXME can be put in the base class too

	if(dd>0) {
		v*=((p1->radius-dd/2.0)/d); //from center of p1 to contact point
		v.normalized();
		ovs.push_back(COverlapping(p1->getpos()+v, (0.5*dd)*v));
		}
	}

void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tbox> *b){
	static vec v;
	static double d, dd;
	for(int i=0; i<6; i++){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		v=b->face[i]->normal_to_point(p1->Xc, p1->radius);//vertical vector from the center of sphere to the plane
		d=v.abs();
		dd=p1->radius-d;
		if(dd>0) ovs.push_back(COverlapping(p1->getpos()+v+(0.5*dd)*b->face[i]->n, (dd/d)*v));
		}
	}

void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tcomposite>  * p2){

	vector<COverlapping> dummy;
	overlaps(dummy, &p1->shell, &p2->shell);
	if(dummy.size()==0)return;

	for(int i=0; i< (p1->elems.size()); i++){
	for(int j=0; j< (p2->elems.size()); j++){
		overlaps(ovs, p1->elems.at(i), p2->elems.at(j));
		}
		}

	for(int j=0; j< ovs.size(); j++){
//		ovs.at(j).x-=p1->Xc; // contact point with respect to center of composit particle
		}
	}

void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tbox>  * b){
	vector<COverlapping> dummy;
	overlaps(dummy, &p1->shell, b);
	if(dummy.size()==0)return;
	for(int i=0; i<p1->elems.size(); i++){
		overlaps(ovs, p1->elems.at(i), b);
		}

	for(int j=0; j< ovs.size(); j++){
		//ovs.at(j).x+=p1->Xc; // contact point with respect to center of composit particle
		}
	}

#endif /* SHAPES_H */

