#ifndef SHAPES_H
#define SHAPES_H 
#include"common.h"
#include"quaternion.h"
using std::cerr; 
using std::endl;
using namespace std;


typedef enum {tsphere, tplane, tbox, tcomposite, tellipsoid} GType;


class GeomObjectBase
	{
	public:
	GeomObjectBase(vec const & v, GType t):Xc(v),Xc0(v),type(t){identifier=1;};
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
	double vol(){return 0;}
	double I(vec n){
		ERROR("Not implemented."); //FIXME
		return 0;}
	
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
		u0(vec(1.0,0.0,0.0)), u1(vec(0.0,1.0,0.0)), u2(vec(0.0,0.0,1.0))
		 {
		face[0]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner,u0));
		face[1]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner,u1));
		face[2]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner,u2));

		face[3]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner+L,-u0));
		face[4]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner+L,-u1));
		face[5]=auto_ptr<GeomObject<tplane> >(new GeomObject<tplane> (corner+L,-u2));
		};

	void shift(const vec &x){corner+=x;}
	void rotate(const vec &n, double alpha){}//FIXME 
	void scale(double scale){L*=scale;};
	void print(std::ostream &out)const{
		out<< identifier<< "   ";
		out<<corner<<"  "<<L;
		}

	double vol(){return L(0)*L(1)*L(2);}
	double I(vec n){
		ERROR("Not implemented."); //FIXME
		return 0;}

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
	GeomObject<tsphere> (const vec &v, double r):GeomObjectBase(v,tsphere){radius=r;}
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

	double vol(){ 
			return 4.0/3.0*radius*radius*radius;}
	double I(vec n){return 2.0/5.0*vol()*radius*radius;}
	private:
	GeomObject<tsphere>();
	};

typedef GeomObject<tsphere> CSphere;
template<>
class GeomObject<tcomposite>: public GeomObjectBase{
	public:	
	GeomObject<tcomposite> (const GeomObject<tcomposite> & p):GeomObjectBase(p.Xc,tcomposite){
		for(int i=0; i<p.elems.size(); i++){
			elems.push_back(new CSphere(*(p.elems.at(i))));
			}
		radius=p.radius;
		Xc=p.Xc;
		}
		
	GeomObject<tcomposite> (const vec &v, double r):GeomObjectBase(v,tcomposite){

		radius=r;
		CSphere *s1=NULL;
		CSphere *s2=NULL, *s3;
		s1=new CSphere(vec(-2*r/3,0.0,0.0), r/3);
		s1->identifier=1;
		s2=new CSphere(vec(0.0), 2*r/3);
		s2->identifier=2;
		s3=new CSphere(vec(2*r/3, 0.0, 0.0), r/3);
		s3->identifier=1;

		if(s1==NULL || s2==NULL){ERROR("error in memory allocation"); exit(1);}
		elems.push_back(s1);
		elems.push_back(s2);
		elems.push_back(s3);

		s3=new CSphere(vec(0.0,2*r/3, 0.0), r/3);
		s3->identifier=1;
		elems.push_back(s3);

		s3=new CSphere(vec(0.0,-2*r/3, 0.0), r/3);
		s3->identifier=1;
		elems.push_back(s3);

		s3=new CSphere(vec(0.0,0.0,2*r/3), r/3);
		s3->identifier=1;
		elems.push_back(s3);

		s3=new CSphere(vec(0.0,0.0,-2*r/3), r/3);
		s3->identifier=1;
		elems.push_back(s3);
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
		}

	double vol(){
		double v=0;
		for(int i=0; i<elems.size(); i++){
			v+=elems.at(i)->vol();
			}
		return v;
		}

	double I(vec n){//FIXME this works only when the elements are on the axes
		n.normalize();
		double II=0;
		for(int i=0; i<elems.size(); i++){
			II+=elems.at(i)->I(n)+elems.at(i)->vol()*(elems.at(i)->Xc0-(elems.at(i)->Xc0*n)*n).abs2();
			}
		return II;
		}

	double min(size_t j){
		double miny=10;
		for(int i=0; i<elems.size(); i++){
			if((elems.at(i)->Xc)(j)-elems.at(i)->radius < miny) miny=((elems.at(i)->Xc)(j)-elems.at(i)->radius);
			}
		return miny;
		};
	void rotateTo(const Quaternion &q){
		for(int i=0; i<elems.size(); i++){
			elems.at(i)->Xc=Xc+q.rotate(elems.at(i)->Xc0);
			}
		};

	void rotate(const vec& n , double alpha){//FIXME
		ERROR("check this");
		//q.setRotation(n, alpha);
		for(int i=0; i<elems.size(); i++){
			//elems.at(i)->Xc0=q.rotate(elems.at(i)->Xc0);
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
	void parse(std::istream &in){ERROR("check this.");};
	
	vector<CSphere *> elems;
	private:
	GeomObject<tcomposite> ();
	};

#include"matrix.h"
using namespace math;
typedef matrix<double> Matrix;

template<>
class GeomObject<tellipsoid>: public GeomObjectBase{
	public:	
		
	GeomObject(const vec &v,double _a, double _b, double _c) :GeomObjectBase(v,tellipsoid), a(_a), b(_b), c(_c) {
		identifier=5;
		radius=tmax(a, tmax(b,c));
		setup();
		}
	~GeomObject(){}

	void setup(double beta=0){

		//Elements of the rotational matrix

		//double beta = M_PI/3.;
		  rotat_mat(2,2)= 1;
		  rotat_mat(1,2)= 0;
		  rotat_mat(2,1)= 0;
		  rotat_mat(0,2)= 0;
		  rotat_mat(2,0)= 0;

		  rotat_mat(0,0)= cos(beta);
		  rotat_mat(0,1)= -sin(beta);
		  rotat_mat(1,0)= sin(beta);
		  rotat_mat(1,1)= cos(beta);

		  //Elements of the scaling matrix

		  scale_mat(0,0)=1.0/(a*a);
		  scale_mat(1,1)=1.0/(b*b);
		  scale_mat(2,2)=1.0/(c*c);

		  ellip_mat=rotat_mat*scale_mat*~rotat_mat;

		}
	void moveto(const vec &v){
		Xc=v;
		}

	void rotateTo(const Quaternion &q){
		//static double beta=0;
		//setup(beta);
		//beta+=0.00005;
			quaternionToMatrix(q, rotat_mat);
			
		  	ellip_mat=rotat_mat*scale_mat*~rotat_mat;
		}

	double I(vec n){//FIXME only in special coordinate system
		return 2.0/5.0*vol()*radius*radius;
		ERROR("check this");
		n.normalize();
		return vol()*(n*scale_mat*n)/5.0;
		}
	double vol(){
		return 4.0/3.0*M_PI*radius*radius*radius;
		return 4.0/3.0*M_PI*a*b*c;
		}
	void rotate(const vec& n , double alpha){//FIXME
		//q.setRotation(n, alpha);
			//orientation=q.rotate(orientation);
		};

	void shift(const vec& v){Xc+=v;};
	void scale(double scale){
			a*=scale;
			b*=scale;
			c*=scale;
		  	scale_mat(0,0)=1.0/(a*a);
		  	scale_mat(1,1)=1.0/(b*b);
		  	scale_mat(2,2)=1.0/(c*c);
		  	ellip_mat=rotat_mat*scale_mat*~rotat_mat;
			radius*=scale;
			};

	void print(std::ostream &out)const{
		out<< identifier<< "   ";
		out<< Xc<< "  "<<radius+0.00001<<"  ";
		out<<ellip_mat(0,0)<< "  "<<ellip_mat(1,1)<< "  "<<ellip_mat(2,2)<< "  ";
		out<<ellip_mat(1,0)<< "  "<<ellip_mat(1,2)<< "  "<<ellip_mat(0,2)<< "  ";
		//out<<-(ellip_mat)(0,1)<< "  "<<-(Xc*ellip_mat)(1)<< "  "<<-(Xc*ellip_mat)(2)<< "  ";
		out<<0<< "  "<<0<< "  "<<0<< "  ";
		//out<<Xc*ellip_mat*Xc-1<<endl;
		out<<-1<<endl;
		}
	
	void parse(std::istream &in){//FIXME
			ERROR("not implemented");
			//in>>identifier;
			}

	Matrix rotat_mat;
	Matrix scale_mat;
	Matrix ellip_mat;
	double a,b,c;
	vec orientation;
	private:
	static int N;
	GeomObject<tellipsoid> (const GeomObject<tcomposite> & p);//not allow copies
	GeomObject<tellipsoid> ();
	};

class COverlapping{
	COverlapping();
	public:
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
		else if(p1->type==tellipsoid&& p2->type==tbox)//FIXME
			overlaps(ovs, static_cast<const GeomObject<tellipsoid> *>(p1), static_cast<const GeomObject<tbox> *>(p2));
		else if(p1->type==tellipsoid&& p2->type==tellipsoid)//FIXME
			overlaps(ovs, static_cast<const GeomObject<tellipsoid> *>(p1), static_cast<const GeomObject<tellipsoid> *>(p2));
		else ERROR("Not Implemented");
		};

	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tbox> *b);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tsphere>  * p2);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tbox> *b);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tellipsoid>  * p2);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tcomposite>  * p2);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tbox>  *b);

	vec x, dx;
	};

inline
void COverlapping::overlaps(vector<COverlapping > &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tsphere>  * p2){
	static vec v;
	static double d, dd;
	v=p2->displacement(p1);//Xc2-Xc1
	d=v.abs();
	dd=p2->radius+p1->radius-d;//FIXME can be put in the base class too

	if(dd>0) {
		v*=((p1->radius-dd/2.0)/d); //from center of p1 to contact point
		v.normalized();
		ovs.push_back(COverlapping(p1->getpos()+v, (dd)*v));
		}
	}

inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tbox> *b){
	static vec v;
	static double d, dd;
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		v=b->face[i]->normal_to_point(p1->Xc, 0);// p1->radius);//vertical vector from the center of sphere to the plane
		d=v.abs();
		dd=p1->radius-d;
		//if(dd>0) ovs.push_back( COverlapping(p1->getpos()+v+(0.5*dd)*b->face[i]->n, (dd/d)*v) );
		if(dd>0) {
				ovs.push_back( COverlapping(p1->getpos()+v*(1+0.5*dd), (dd)*v) );
				}
		}
	}
inline
void COverlapping::overlaps(vector<COverlapping > &ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tellipsoid>  * p2){
	static vec v;
	static double d, dd;
	v=p2->displacement(p1);//Xc2-Xc1
	d=v.abs();
	dd=p2->radius+p1->radius-d;//FIXME can be put in the base class too

	if(dd>0) {
		v*=((p1->radius-dd/2.0)/d); //from center of p1 to contact point
		v.normalized();
		ovs.push_back(COverlapping(p1->getpos()+v, (dd)*v));
		}
	}
inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tbox> *b){
	static vec v;
	static double d, dd;
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		v=b->face[i]->normal_to_point(p1->Xc, 0);// p1->radius);//vertical vector from the center of sphere to the plane
		d=v.abs();
		dd=p1->radius-d;
		//if(dd>0) ovs.push_back( COverlapping(p1->getpos()+v+(0.5*dd)*b->face[i]->n, (dd/d)*v) );
		if(dd>0) {
				ovs.push_back( COverlapping(p1->getpos()+v*(1+0.5*dd), (dd)*v) );
				}
		}
	}

inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tcomposite>  * p2){

	if(p1==p2){
		ERROR("A particle is checked against itself for overlapping.")
		return;
		}
	if((p1->Xc-p2->Xc).abs() > p1->radius+p2->radius)return;

	for(int i=0; i< (p1->elems.size()); ++i){
	for(int j=0; j< (p2->elems.size()); ++j){
		overlaps(ovs, p1->elems.at(i), p2->elems.at(j));
		}
		}
	}

inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tbox>  * b){
	static double d;
	bool need_to_check=false;
	for(int i=0; i<6; ++i){
		d=(b->face[i]->normal_to_point(p1->Xc,0.0)).abs2() - p1->radius*p1->radius;
		if(d<0){
			need_to_check=true;
			break;
			}
		}
//	if(!need_to_check)return;

	for(int i=0; i<p1->elems.size(); ++i){
		overlaps(ovs, p1->elems.at(i), b);
		}
	}

typedef GeomObject<tellipsoid> CEllipsoid;
vec rescale_ellipse_to_touch(const CEllipsoid &E, const CEllipsoid &E0){
//here we should if E will touch E0 or not;

static vec R, R0, R12;
R12=E0.Xc-E.Xc; 
R=E.Xc; //FIXME unnecessary, just not to modify the rest for the moment
R0=E0.Xc;


// Begin of iteration for lambda
static Matrix D(3,3), invD(3,3), F(3,3), A(3,3), B(3,3);
static Matrix invE(3,3), invE0(3,3);
invE=!E.ellip_mat;
invE0=!E0.ellip_mat;
int itermax=1000;
double lambda=0.5;
for(int loop=1;loop<=itermax;loop++){

	D=invE0+lambda*invE;	
	invD=!D;
//these two lines can be optimized
	double d=lambda*lambda*invE.Det()+lambda*((invE+invE0).Det()-invE.Det()-invE0.Det())+invE0.Det();
	double dd=2*lambda*invE.Det()+((invE+invE0).Det()-invE.Det()-invE0.Det());
	
	F=invE*invE0*invD;
	A=invD*invE0*invD;
	
	B=2.0/d*(F-dd*A);



	double dx= (R*A*R)/(R*B*R);

	double threshold=1e-4;
	static Matrix E12(3,3);
	static vec rc;
	static vec rE;
	if(fabs(dx) < threshold){
		//CVector rc;
		E12=!(E.ellip_mat-(lambda*E0.ellip_mat));
		rc=R0-E12*(E.ellip_mat*R12);
		double s=lambda*sqrt(R*(E0.ellip_mat*E12*E.ellip_mat*E12*E0.ellip_mat)*R);
		rE=(1-1/s)*rc+(1/s)*R;
		cerr<< s <<endl;
		return rE;
		}
	lambda -= dx;
	}   

  return E.Xc;   //NON Convergence.

}

#endif /* SHAPES_H */
