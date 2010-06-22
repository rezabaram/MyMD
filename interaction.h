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
	static void overlaps(ShapeContact* ovs, GeomObject<tellipsoid>   *p1, const GeomObject <tplane>   *plane);
	static void overlaps(ShapeContact* ovs, GeomObject <tellipsoid>  *p1, const GeomObject <tbox>        *b );
	static void overlaps(ShapeContact* ovs, GeomObject <tellipsoid>  *p1, GeomObject <tellipsoid>  *p2);
	static void overlaps(ShapeContact* ovs, const GeomObject <tcomposite>  *p1, const GeomObject <tcomposite>  *p2);
	static void overlaps(ShapeContact* ovs, const GeomObject <tcomposite>  *p1, const GeomObject <tbox>        *b );

	static void append(ShapeContact&v, ShapeContact&v2);

	};

inline
void CInteraction::overlaps(ShapeContact* ovs, const GeomObject<tsphere>  *p1, const GeomObject<tsphere>  * p2){
	static vec v;
	static double d, dd;
	v=p2->Xc-p1->Xc;
	d=v.abs();
	dd=p2->radius+p1->radius-d;//FIXME can be put in the base class too

	if(dd>0) {
		v*=((p1->radius-dd/2)/d); //from center of p1 to contact point
		ovs->add(Contact(p1->getpos()+v, v.normalized(), dd));
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
			overlaps(ovs, static_cast<GeomObject<tellipsoid> *>(p1), static_cast<const GeomObject<tbox> *>(p2));
		else if(p1->type==tellipsoid&& p2->type==tellipsoid)//FIXME
			overlaps(ovs, static_cast<GeomObject<tellipsoid> *>(p1), static_cast<GeomObject<tellipsoid> *>(p2));
		else ERROR(true, "Not Implemented");
		};
inline
void CInteraction::overlaps(ShapeContact* ovs, const GeomObject<tsphere>  *p1, const GeomObject<tbox> *b){
	static vec v;
	static double d, dd;
	for(size_t i=0; i<b->nFaces; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		v=b->face[i]->normal_from_point(p1->Xc, 0);// p1->radius);//vertical vector from the center of sphere to the plane
		d=v.abs();
		dd=p1->radius-d;
		//if(dd>0) ovs->push_back( CInteraction(p1->getpos()+v+(0.5*dd)*b->face[i]->n, (dd/d)*v) );
		if(dd>0) {
			ovs->add( Contact(p1->getpos()+v+1.0*dd*v.normalized(), v.normalized(), dd) );
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
	for(indexType i=0; i<b->nFaces; ++i){
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
void CInteraction::overlaps(ShapeContact* ovs, GeomObject<tellipsoid>  *p1, const GeomObject<tplane> *plane){
	static vec v, vp;
	static double dx;
	if(plane->normal_from_point(p1->Xc).abs()-p1->radius > 0) return;
	vp=p1->point_to_plane(*(plane));
	v=plane->normal_to_point(vp, 0);
	dx=v.abs();
	v/=dx;
	if(v*plane->n >0)return;
	ovs->add(Contact(vp, v, dx) );
	p1->fixToBody(HomVec(vp,1));
	return;
	}

inline
void CInteraction::overlaps(ShapeContact* ovs, GeomObject<tellipsoid>  *p1, const GeomObject<tbox> *b){
TRY
	for(size_t i=0; i<b->nFaces; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		overlaps(ovs, p1, b->face[i]);
		}
CATCH
	}

void fixcontact(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2){
TRY
	ovs.x01=E1.toBody(ovs.x1);
	ovs.x02=E2.toBody(ovs.x2);
CATCH
	}
void updatecontact(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2){
TRY
	ovs.x1=E1.toWorld(ovs.x01);
	ovs.x2=E2.toWorld(ovs.x02);
CATCH
	}

void updateplane(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2){
TRY
	ovs.plane.Xc=(ovs.x1.project()+ovs.x2.project())/2;
	ovs.plane.n=(E1.gradient(ovs.x1.project())-E2.gradient(ovs.x2.project())).normalized();
CATCH
	}

void correctpoints(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2){
TRY
	ovs.x1=HomVec(E1.point_to_plane((ovs.plane)),1);
	ovs.x2=HomVec(E2.point_to_plane((ovs.plane)),1);
	ovs.x01=E1.toBody(ovs.x1);
	ovs.x02=E2.toBody(ovs.x2);
	//static int i=0;
	//if( E2(ovs.x1) > epsilon or E1(ovs.x2) >epsilon ) cout<<"  "<< ++i<<"  "<<E2(ovs.x1) <<" "<<E1(HomVec(ovs.plane.Xc,1))<<" "<<E2(HomVec(ovs.plane.Xc,1))<<"  "<<E1(ovs.x2)<<endl;
CATCH
	}

void setcontact(ShapeContact &ovs,CEllipsoid  &E1, CEllipsoid  &E2){
TRY
	//ovs.add(Contact(ovs.plane.Xc, ovs.plane.n, fabs((ovs.x1.project()-ovs.x2.project())*ovs.plane.n)));
	ovs.x01=E1.toBody(ovs.x1);
	ovs.x02=E2.toBody(ovs.x2);
	ovs.add(Contact((ovs.x1.project()+ovs.x2.project())/2, (ovs.x1.project()-ovs.x2.project()).normalized(), (ovs.x1-ovs.x2).abs()));
	//ovs.add(Contact((ovs.x1.project()+ovs.x2.project())/2, (ovs.x1.project()-ovs.x2.project()).normalized(), (ovs.x1-ovs.x2).abs()));
CATCH
	}


bool separatingPlane(ShapeContact &ovs,  CEllipsoid  &E1, CEllipsoid  &E2){
TRY

	if(0)if(ovs.has_sep_plane){
		if(!(E1.doesHit(ovs.plane) or E2.doesHit(ovs.plane))) {
			return true;
			}
		}
	
	Matrix M=(-(!E1.ellip_mat)*E2.ellip_mat);
	//CQuartic q=characteristicPolynom(M);

	vector<complex<double> > eigenvals;
	vector<HomVec> eigenvecs;
	eigens(M, eigenvals, eigenvecs);


	if(fabs(eigenvals.at(3).imag() ) < epsilon){
		ERROR(eigenvals.at(2).imag()>epsilon,"one eigenvalue complex one not.");

		CRay<HomVec> ray(eigenvecs.at(3).project4d(),eigenvecs.at(2).project4d());
		CQuadratic q1(intersect(ray, E1));
		CQuadratic q2(intersect(ray, E2));
		//the roots are sorted ascending
		HomVec X1= ray(q1.root(1).real()).project4d(); //on the surface of E1
		HomVec X2= ray(q2.root(0).real()).project4d(); //on the surface of E2

		ovs.x1=X1;
		ovs.x2=X2;
		ovs.x01=E1.toBody(X1);
		ovs.x02=E2.toBody(X2);


		//calculating the separating plane

		vec n1=HomVec(E1.ellip_mat*X1).project();
		vec n2=HomVec(E2.ellip_mat*X2).project();
		HomVec mp=(X1+X2)*.5;
		double alpha=-(mp.project()-X1.project())*n1 /( (mp.project()-X2.project())*n2);
		ovs.plane=CPlane(mp.project(), n1+alpha*n2);
		cerr<< ovs.plane.n*(E2.Xc-E1.Xc).normalize() <<endl<<endl;;

		bool hit=(E1.doesHit(ovs.plane) or E2.doesHit(ovs.plane));
		if(hit){
			cerr<< E1(ovs.plane.Xc) <<"\t"<< E2(ovs.plane.Xc) <<endl;
			}
		ERROR(hit, "the separation plane set incorrectly");
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

// find min of x on E1, in the potentional of E2
void findMin(HomVec &x,  CEllipsoid  &E1, CEllipsoid  &E2, long nIter=1){
TRY
	double lambda;
	vec xp0, xp=x.project();
	Matrix Em1(3,3), Em2(3,3);
	for(size_t i=0; i<3; i++){
	for(size_t j=0; j<3; j++){
		Em1(i,j)=E1.ellip_mat(i,j);
		Em2(i,j)=E2.ellip_mat(i,j);
		}
		}
	
	long iter=0;
	bool converged=false;
	do{
		++iter;
		xp0=xp;
		lambda=fabs((xp-E1.Xc)*E2.ellip_mat*(xp-E2.Xc));
		xp=(!(Em2+lambda*Em1))*(Em2*E2.Xc+lambda*Em1*E1.Xc);
		converged= (xp-xp0).abs()<1e-14;
		}
	while(iter<nIter and !converged);
	//	cout<<setprecision(14)<<iter<<"  lambda="<< lambda<<"   Xp="<<xp <<"  C1="<< E1.Xc<<"  C2="<<E2.Xc<<endl;
	

	if(!converged)WARNING("minimization not converged: "<<E2(xp));
	if(converged and E2(xp) > 0)WARNING("Coverged to maximum instread of minimum: "<<E2(xp));
	x(0)=xp(0);
	x(1)=xp(1);
	x(2)=xp(2);
CATCH
}

#define XOR(p, q) ( ((p) || (q)) && !((p) && (q)) ) 
inline
void CInteraction::overlaps(ShapeContact* ovs, CEllipsoid  *E1, CEllipsoid  *E2){
TRY
	
	ERROR(E1==E2, "A particle is checked for overlapping against itself for overlapping.")
	if((E1->Xc-E2->Xc).abs()>1.001*(E1->radius+E2->radius))return;

	//if(ovs->set)if( !E1->doesHit(ovs->plane) and !E2->doesHit(ovs->plane)) return;
	 //ovs->set=true;

	ovs->has_sep_plane=separatingPlane(*ovs, *E1, *E2);
	if(!ovs->has_sep_plane){
		
		//cerr<< E1->doesHit(ovs->plane) <<"\t"<< E2->doesHit(ovs->plane) <<endl;
		//ERROR(( !E1->doesHit(ovs->plane) and !E2->doesHit(ovs->plane)), " ");

		updatecontact(*ovs, *E1, *E2);
		findMin(ovs->x1, *E1, *E2, 50);
		findMin(ovs->x2, *E2, *E1, 50);
		//updateplane(*ovs, *E1, *E2);
		//correctpoints(*ovs, *E1, *E2);
		setcontact(*ovs, *E1, *E2);

		}
	else{
		//dont collid
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
