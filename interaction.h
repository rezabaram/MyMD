#ifndef OVERLAPPING_H
#define OVERLAPPING_H 
#include"shapes.h"
#include"eigen.h"

class Contact{
	public:
	Contact(const vec &_x, const vec &_dx ):x(_x), dx(_dx){}
	vec x, dx;
	};

template<typename T1=double, typename T2=double>
class PairContact : public Contact {
	public:
	PairContact(const vec &_x, const vec &_dx ):Contact(_x,_dx){}
 	private:
	};

class ShapeContactHolder : public vector<Contact*>{
	public:
	ShapeContactHolder():vector<Contact*>(){}
	void add(const Contact &c){
		push_back(new Contact(c));
		}
	~ShapeContactHolder(){
		ShapeContactHolder::iterator it;
		for(it=this->begin(); it!=this->end(); ++it)
			delete (*it);
		}
	};

class CInteraction{
	public:
	CInteraction();

	static void overlaps(ShapeContactHolder* ovs, GeomObjectBase *p1, GeomObjectBase *p2);

	static void append(ShapeContactHolder&v, ShapeContactHolder&v2);
	//every new kind of particle needs to define two functions
	static void overlaps(ShapeContactHolder* ovs, const GeomObject <tsphere>     *p1, const GeomObject <tbox>        *b );
	static void overlaps(ShapeContactHolder* ovs, const GeomObject <tsphere>     *p1, const GeomObject <tsphere>     *p2);
	static void overlaps(ShapeContactHolder* ovs, const GeomObject <tellipsoid>  *p1, const GeomObject <tbox>        *b );
	static void overlaps(ShapeContactHolder* ovs, GeomObject <tellipsoid>  *p1, GeomObject <tellipsoid>  *p2);
	static void overlaps(ShapeContactHolder* ovs, const GeomObject <tcomposite>  *p1, const GeomObject <tcomposite>  *p2);
	static void overlaps(ShapeContactHolder* ovs, const GeomObject<tellipsoid>   *p1, const GeomObject <tplane>   *plane);
	static void overlaps(ShapeContactHolder* ovs, const GeomObject <tcomposite>  *p1, const GeomObject <tbox>        *b );


	};

inline
void CInteraction::overlaps(ShapeContactHolder* ovs, const GeomObject<tsphere>  *p1, const GeomObject<tsphere>  * p2){
	static vec v;
	static double d, dd;
	v=p2->displacement(p1);//Xc2-Xc1
	d=v.abs();
	dd=p2->radius+p1->radius-d;//FIXME can be put in the base class too

	if(dd>0) {
		v*=((p1->radius-dd/2.0)/d); //from center of p1 to contact point
		v.normalized();
		ovs->add(Contact(p1->getpos()+v, (dd)*v));
		}
	}

inline
void CInteraction::overlaps(ShapeContactHolder* ovs, GeomObjectBase *p1, GeomObjectBase *p2){
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
			overlaps(ovs, static_cast<GeomObject<tellipsoid> *>(p1), static_cast<GeomObject<tellipsoid> *>(p2));
		else ERROR(true, "Not Implemented");
		};
inline
void CInteraction::overlaps(ShapeContactHolder* ovs, const GeomObject<tsphere>  *p1, const GeomObject<tbox> *b){
	static vec v;
	static double d, dd;
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		v=b->face[i]->normal_from_point(p1->Xc, 0);// p1->radius);//vertical vector from the center of sphere to the plane
		d=v.abs();
		dd=p1->radius-d;
		//if(dd>0) ovs->push_back( CInteraction(p1->getpos()+v+(0.5*dd)*b->face[i]->n, (dd/d)*v) );
		if(dd>0) {
			ovs->add(Contact(p1->getpos()+v*(1+0.5*dd), (dd)*v) );
			}
		}
	}


inline
void CInteraction::overlaps(ShapeContactHolder* ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tplane> *plane){
	static vec v, vp;
	if(plane->normal_from_point(p1->Xc).abs()-p1->radius > 0) return;
	vp=p1->point_to_plane(*(plane));
	v=plane->normal_from_point(vp, 0);
	if(v*plane->n <0)return;
	ovs->add(  Contact(vp, -v) );
	return;
	}

inline
void CInteraction::overlaps(ShapeContactHolder* ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tbox> *b){
	try{
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		overlaps(ovs, p1, b->face[i]);
		}
	}catch(...){
		ERROR(1,"Some error in the interaction with wall");
		}
	}

inline
void CInteraction::overlaps(ShapeContactHolder* ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tcomposite>  * p2){

	ERROR(p1==p2, "A particle is checked against itself for overlapping.")
	if((p1->Xc-p2->Xc).abs() > p1->radius+p2->radius)return;

	for(indexType i=0; i< (p1->elems.size()); ++i){
		for(indexType j=0; j< (p2->elems.size()); ++j){
			overlaps(ovs, p1->elems.at(i), p2->elems.at(j));
			}
		}
	}

inline
void CInteraction::overlaps(ShapeContactHolder* ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tbox>  * b){
	static double d;
	bool need_to_check=false;
	for(indexType i=0; i<6; ++i){
		d=(b->face[i]->normal_from_point(p1->Xc,0.0)).abs2() - p1->radius*p1->radius;
		if(d<0){
			need_to_check=true;
			break;
			}
		}
//	if(!need_to_check)return;

	for(indexType i=0; i<p1->elems.size(); ++i){
		overlaps(ovs, p1->elems.at(i), b);
		}
	}

	bool separatingPlane(CPlane &plane,  CEllipsoid  &E1, CEllipsoid  &E2){
	Matrix M=(-(!E1.ellip_mat)*E2.ellip_mat);
	//CQuartic q=characteristicPolynom(M);

	vector<complex<double> > eigenvals;
	vector<HomVec> eigenvecs;
	eigens(M, eigenvals, eigenvecs);


	if(fabs(eigenvals.at(3).imag() ) < epsilon){
		CRay<HomVec> ray(eigenvecs.at(3),eigenvecs.at(2));
		CQuadratic q1(intersect(ray, E1));
		CQuadratic q2(intersect(ray, E2));
		//the roots are sorted ascending
		HomVec X1= ray(q1.root(1).real()); //on the surface of E1
		HomVec X2= ray(q2.root(0).real()); //on the surface of E2
		E1.fixToBody(X1);
		E2.fixToBody(X2);

		vec n1=HomVec(E1.ellip_mat*X1).project();
		vec n2=HomVec(E2.ellip_mat*X2).project();
/*
		CPlane p1(X1.project(),n1);
		CPlane p2(X2.project(),n2);
		CSphere S((X1.project()+X2.project())*.5, 0.01);
		CRay<vec> nr(S.Xc, S.Xc+0.05*(n1-n2));
		if(gout!=NULL and gout->is_open()){
			S.print(*gout);
			(*gout)<<endl;
			
			nr.print(*gout);
			}
*/

		plane=CPlane((X1.project()+X2.project())*.5, (n1-n2)*.5);
		return true;
		}
	return false;
	}

void adjust2(CEllipsoid &E1, CEllipsoid &E2, long nIter=1){
TRY
	for(long i=0; i<nIter; i++){
        	HomVec X=(E1.P+E2.P)/2;
        	HomVec X1, X2;
	 	CRay<HomVec> ray1(HomVec(E1.Xc,1), HomVec(X(0), X(1), X(2), X(3))); 
	 	CRay<HomVec> ray2(HomVec(E2.Xc,1), HomVec(X(0), X(1), X(2), X(3))); 
		X1=ray1((intersect(ray1, E1)).root(1).real());
		X2=ray2((intersect(ray2, E2)).root(1).real());
		E1.fixToBody(X1);
		E2.fixToBody(X2);
		}
	 
	return;
CATCH
	}
void adjust(CEllipsoid &E1, const CEllipsoid &E2, long nIter=20){
TRY
	double dx=(E1.P-E2.P).abs();	
        HomVec X=E1.P, Xp, Xpp;
	bool b=E2(X)<0?true:false;
	for(long i=0; i<nIter; i++){
	 	CRay<HomVec> raytest(HomVec(E1.Xc,1), HomVec(X(0)+dx*drand48(), X(1)+dx*drand48(), X(2)+dx*drand48(),X(3))); 
		Xp=raytest((intersect(raytest, E1)).root(1).real());
		Xpp=raytest((intersect(raytest, E2)).root(0).real());
		//cerr<< " here " <<E2(Xp) << "  " <<E1(Xpp) <<endl;
		//cerr<< " here " <<E2(Xp) * E1(Xpp) <<endl;
		if(E2(Xp)<E2(X)){
			cerr<<setprecision(14)<< (Xp*X)/Xp.abs()/X.abs() <<endl;
			//cerr<< (Xp-X).abs() <<endl;
			X=Xp;
			}
		}
	if(b and E2(X)>0)cerr<< "errorooooooooooooo in "<<__FILE__ <<endl;
	E1.fixToBody(X);
	 
	return;
CATCH
	}

void CInteraction::append(ShapeContactHolder &v, ShapeContactHolder &v2){
	for(size_t i=0; i<v2.size(); i++){
		v.push_back(v2.at(i));	
		}

	return;
	}
#define XOR(p, q) ( ((p) || (q)) && !((p) && (q)) ) 
inline
void CInteraction::overlaps(ShapeContactHolder* ovs, CEllipsoid  *E1, CEllipsoid  *E2){
TRY
	static CPlane plane(vec(0,0,0), vec(0,0,1));
        static HomVec X1=E1->P, X2=E2->P;
	
	//bool b=(*E2)(X1)>0?true:false;
	if(!separatingPlane(plane, *E1, *E2)){
		//adjust2(*E1, *E2);
		adjust(*E1, *E2);
		adjust(*E2, *E1);
	//if(b and (*E2)(X1)>0)cerr<< "errorooooooooooooo "<<++i<<endl;
		//ovs->push_back( Contact( (E1->P.project()+E2->P.project())/2, (E1->P.project()-E2->P.project())) );
		double dd=(E1->P.project()-E2->P.project()).abs();
		vec x=(E1->P.project()+E2->P.project())/2;
		vec dx=(E1->gradient(x)-E2->gradient(x));
		dx.normalize();
		ovs->add(Contact(x,dd*dx) );
		}
/*
	overlaps(ovtest1, E, p);
	overlaps(ovtest2, E0, p);
	if(ovtest1.size()>0 and ovtest2.size()>0){
		//append(ovs, ovtest1);
		throw 1;
		collide=true;
		}

	*p= separatingPlane(*E, *E0);//update the plane
	if(XOR(ovtest1.size()>0 , ovtest2.size()>0) and !collide){
		 *p= separatingPlane(*E, *E0);//update the plane
		}
	
	//ERROR(1, "stopped");
	//ERROR(1,"Stopped here");
*/
CATCH
	}
/*
inline
void overlapsTemp(vector<Contact> ovs, const CEllipsoid  *E, const CEllipsoid  *E0){
	//here we should if E will touch E0 or not;
// ------------------------------------------
	static vec v;
	static double d;
	v=E0->displacement(E);//Xc2-Xc1
	d=v.abs();
	double dp=E0->radius+E->radius-d;//FIXME can be put in the base class too

	if(dp<0) {//return;
		v*=((E->radius-dp/2.0)/d); //from center of p1 to contact point
		v.normalized();
		ovs->push_back(Contact(E->getpos()+v, (dp)*v));
		}
//------------------------------------
	static vec R, R0, R12;
	R12=E0->Xc-E->Xc; 
	R=E->Xc; //FIXME unnecessary, just not to modify the rest for the moment
	R0=E0->Xc;


	// Begin of iteration for lambda
	static Matrix D(3,3), invD(3,3), F(3,3), A(3,3), B(3,3);
	static Matrix invE(3,3), invE0(3,3);
	invE=E->inv();
	invE0=E0->inv();

	int itermax=1000;
	double lambda=0.0;
	double dx;
	for(int loop=1;loop<=itermax;loop++){

		for(lambda=-8; lambda<8; lambda+=0.001){

		D=invE0+lambda*invE;	
		invD=!D;
	//these two lines can be optimized
		double d=lambda*lambda*(invE.Det())+lambda*((invE+invE0).Det()-invE.Det()-invE0.Det())+invE0.Det();
		double dd=2*lambda*(invE.Det())+((invE+invE0).Det()-invE.Det()-invE0.Det());
		
		F=invE*invE0*invD;
		A=invD*invE0*invD;
		
		B=2.0/d*(F-dd*A);

		dx= (R12*A*R12 - 1.0)/(R12*B*R12);
		cout<<lambda<<"   "<<(R12*A*R12 - 1.0)<<" vector ("<< invD*R12<<")   "<<invD.Det()<<endl;
		}
		exit(0);

		double threshold=1e-6;
		static Matrix E12(3,3);
		static vec rc;
		static vec rE;
		if(fabs(dx) < threshold){
			//CVector rc;
			E12=!(E->ellip_mat+(lambda*E0->ellip_mat));
			rc=R0-E12*(E->ellip_mat*R12);
			double s2=lambda*lambda*(R12*(E0->ellip_mat*E12*E->ellip_mat*E12*E0->ellip_mat)*R12);
			double s=sqrt(s2);
			cerr<< s2 <<endl;
			if(s>1){
				break;//not overlapping
				}
			rE=(1-1/s)*rc+(1/s)*R;
			//cerr<< s <<endl;
			
			//return CContact(rc-R, rE-R);
			cerr<< dx<<"  "<<lambda <<endl;
			cerr<< setprecision(10)<<R <<"   x "<<rc<<"  dx "<<rE <<endl;
			//ovs->push_back(Contact(rc, R-rE));
			return;
			}
		lambda += dx;
		}   
//	cerr<< dx <<endl;

	  return ;

}
*/
#endif /* OVERLAPPING_H */
