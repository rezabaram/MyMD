#ifndef OVERLAPPING_H
#define OVERLAPPING_H 
#include "shapes.h"
#include "eigen.h"
#include"shapecontact.h"


class CInteraction{
	public:
	CInteraction();

	static void overlaps(ShapeContact* ovs, GeomObjectBase *p1, GeomObjectBase *p2);

	//every new kind of particle needs to define two functions
	static void overlaps(ShapeContact* ovs, const GeomObject <tsphere>     *p1, const GeomObject <tbox>        *b );
	static void overlaps(ShapeContact* ovs, const GeomObject <tsphere>     *p1, const GeomObject <tsphere>     *p2);
	static void overlaps(ShapeContact* ovs, const GeomObject <tellipsoid>  *p1, const GeomObject <tbox>        *b );
	static void overlaps(ShapeContact* ovs, GeomObject <tellipsoid>  *p1, GeomObject <tellipsoid>  *p2);
	static void overlaps(ShapeContact* ovs, const GeomObject <tcomposite>  *p1, const GeomObject <tcomposite>  *p2);
	static void overlaps(ShapeContact* ovs, const GeomObject<tellipsoid>   *p1, const GeomObject <tplane>   *plane);
	static void overlaps(ShapeContact* ovs, const GeomObject <tcomposite>  *p1, const GeomObject <tbox>        *b );

	static void append(ShapeContact&v, ShapeContact&v2);

	};

inline
void CInteraction::overlaps(ShapeContact* ovs, const GeomObject<tsphere>  *p1, const GeomObject<tsphere>  * p2){
	static vec v;
	static double d, dd;
	v=p2->displacement(p1);//Xc2-Xc1
	d=v.abs();
	dd=p2->radius+p1->radius-d;//FIXME can be put in the base class too

	if(dd>0) {
		v*=((p1->radius-dd/2.0)/d); //from center of p1 to contact point
		v.normalized();
		ovs->add(Contact(p1->getpos()+v, v, dd));
		}
	}

inline
void CInteraction::overlaps(ShapeContact* ovs, GeomObjectBase *p1, GeomObjectBase *p2){
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
void CInteraction::overlaps(ShapeContact* ovs, const GeomObject<tsphere>  *p1, const GeomObject<tbox> *b){
	static vec v;
	static double d, dd;
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		v=b->face[i]->normal_from_point(p1->Xc, 0);// p1->radius);//vertical vector from the center of sphere to the plane
		d=v.abs();
		dd=p1->radius-d;
		//if(dd>0) ovs->push_back( CInteraction(p1->getpos()+v+(0.5*dd)*b->face[i]->n, (dd/d)*v) );
		if(dd>0) {
			ovs->add(Contact(p1->getpos()+v*(1+0.5*dd), v, dd) );
			}
		}
	}



inline
void CInteraction::overlaps(ShapeContact* ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tcomposite>  * p2){
TRY

	ERROR(p1==p2, "A particle is checked against itself for overlapping.")
	if((p1->Xc-p2->Xc).abs() > p1->radius+p2->radius)return;

	for(indexType i=0; i< (p1->elems.size()); ++i){
		for(indexType j=0; j< (p2->elems.size()); ++j){
			overlaps(ovs, p1->elems.at(i), p2->elems.at(j));
			}
		}
CATCH
	}

inline
void CInteraction::overlaps(ShapeContact* ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tbox>  * b){
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

inline
void CInteraction::overlaps(ShapeContact* ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tplane> *plane){
	static vec v, vp;
	static double dx;
	if(plane->normal_from_point(p1->Xc).abs()-p1->radius > 0) return;
	vp=p1->point_to_plane(*(plane));
	v=plane->normal_from_point(vp, 0);
	dx=v.abs();
	v/=dx;
	if(v*plane->n <0)return;
	ovs->add(Contact(vp, -v, dx) );
	return;
	}

inline
void CInteraction::overlaps(ShapeContact* ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tbox> *b){
TRY
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		overlaps(ovs, p1, b->face[i]);
		}
CATCH
	}
void adjust3(ShapeContact &ovs, CEllipsoid &E1, CEllipsoid &E2, long nIter=0){
TRY
		vec x=(ovs.x1.project()+ovs.x2.project())/2;
		vec dx=(E1.gradient(x)-E2.gradient(x));
		dx.normalize();
	for(long i=0; i<nIter; i++){
        	HomVec X=(ovs.x1+ovs.x2)/2;
        	HomVec X1, X2;
	 	CRay<HomVec> ray1(HomVec(x,1), HomVec(x+dx,1)); 
	 	CRay<HomVec> ray2(HomVec(x,1), HomVec(x-dx,1)); 
		X1=ray1((intersect(ray1, E1)).root(1).real());
		X2=ray2((intersect(ray2, E2)).root(1).real());
		E1.P=X1;
		E2.P=X2;
		ovs.x1=X1;
		ovs.x2=X2;
		ovs.x01=E1.toBody(X1);
		ovs.x02=E2.toBody(X2);
		}
	 
		double dd=(ovs.x1-ovs.x2).abs();
		ovs.add(Contact(x,dx,dd) );
	return;
CATCH
	}
void adjust1(ShapeContact &ovs, CEllipsoid &E1, const CEllipsoid &E2, long nIter=0){
TRY
        HomVec X=ovs.x1, Xp, Xpp;
	double dx=(X-ovs.x2).abs();	
	bool b=E2(X)<0?true:false;
	
	for(long i=0; i<nIter; i++){
	 	CRay<HomVec> raytest(HomVec(E1.Xc,1), HomVec(X(0)+dx*drand48(), X(1)+dx*drand48(), X(2)+dx*drand48(),X(3))); 
		Xp=raytest((intersect(raytest, E1)).root(1).real());
		Xpp=raytest((intersect(raytest, E2)).root(0).real());
		if(E2(Xp)<E2(X)){
			//cerr<<setprecision(14)<< (Xp*X)/Xp.abs()/X.abs() <<endl;
			X=Xp;
			}
		}
	if(b and E2(X)>0)cerr<< "errorooooooooooooo in "<<__FILE__ <<endl;
	ovs.x1=X;
	ovs.x01=E1.toBody(X);
	 
	return;
CATCH
	}
void adjust2(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2, long nIter=0){
TRY
        HomVec X=ovs.x2, Xp, Xpp;
	double dx=(X-ovs.x1).abs();	
	bool b=E1(X)<0?true:false;
	
	for(long i=0; i<nIter; i++){
	 	CRay<HomVec> raytest(HomVec(E2.Xc,1), HomVec(X(0)+dx*drand48(), X(1)+dx*drand48(), X(2)+dx*drand48(),X(3))); 
		Xp=raytest((intersect(raytest, E2)).root(1).real());
		Xpp=raytest((intersect(raytest, E1)).root(0).real());
		if(E1(Xp)<E1(X)){
	//		cerr<<setprecision(14)<< (Xp*X)/Xp.abs()/X.abs() <<endl;
			X=Xp;
			}
		}
	if(b and E1(X)>0)cerr<< "errorooooooooooooo in "<<__FILE__ <<endl;
	ovs.x2=X;
	ovs.x02=E2.toBody(X);
	 
	return;
CATCH
	}

bool hit_plane(const CEllipsoid  *p1, const CPlane *plane){
TRY
	static vec v, vp;
	//if(plane->normal_from_point(p1->Xc).abs()-p1->radius > 0) return false;
	vp=p1->point_to_plane(*(plane));
	v=plane->normal_from_point(vp, 0);
	if(v*plane->n <0)return false;
	return true;
CATCH
	}
bool separatingPlane(ShapeContact &ovs,  CEllipsoid  &E1, CEllipsoid  &E2){
TRY
	//if(ovs.set)if( !hit_plane(&E1, &ovs.plane) ) return true;
	//if(ovs.set)if( !hit_plane(&E1, &ovs.plane) and !hit_plane(&E2, &ovs.plane)) return false;
	 //ovs.set=true;
	
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

		ovs.x1=X1;
		ovs.x2=X2;
		ovs.x01=E1.toBody(X1);
		ovs.x02=E2.toBody(X2);

		vec n1=HomVec(E1.ellip_mat*X1).project();
		vec n2=HomVec(E2.ellip_mat*X2).project();

		ovs.plane=CPlane((X1.project()+X2.project())*.5, (n1-n2)*.5);
		return true;
		}
	return false;
CATCH
	}


void CInteraction::append(ShapeContact &v, ShapeContact &v2){
	for(size_t i=0; i<v2.size(); i++){
		v.push_back(v2.at(i));	
		}

	return;
	}
#define XOR(p, q) ( ((p) || (q)) && !((p) && (q)) ) 
inline
void CInteraction::overlaps(ShapeContact* ovs, CEllipsoid  *E1, CEllipsoid  *E2){
TRY
	
	ERROR(E1==E2, "A particle is checked for overlapping against itself for overlapping.")
	if((E1->Xc-E2->Xc).abs()>1.01*(E1->radius+E2->radius))return;

	if(!separatingPlane(*ovs, *E1, *E2)){
		//adjust2(*E1, *E2);
		
		ovs->x1=E1->toWorld(ovs->x01);
		ovs->x2=E2->toWorld(ovs->x02);

		adjust1(*ovs, *E1, *E2, 10);
		adjust2(*ovs, *E1, *E2, 10);
		adjust3(*ovs, *E1, *E2, 1);
		}
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
