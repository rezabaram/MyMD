#ifndef OVERLAPPING_H
#define OVERLAPPING_H 
#include"shapes.h"
#include"eigen.h"


class COverlapping{
	COverlapping();
	public:
	COverlapping(const vec &_x, const vec &_dx ):x(_x), dx(_dx){}

	static void overlaps(vector<COverlapping> &ovs, const GeomObjectBase *p1, const GeomObjectBase *p2);


	//every new kind of particle needs to define two functions
	static void overlaps(vector<COverlapping> &ovs, const GeomObject <tsphere>     *p1, const GeomObject <tbox>        *b );
	static void overlaps(vector<COverlapping> &ovs, const GeomObject <tsphere>     *p1, const GeomObject <tsphere>     *p2);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject <tellipsoid>  *p1, const GeomObject <tbox>        *b );
	static void overlaps(vector<COverlapping> &ovs, const GeomObject <tellipsoid>  *p1, const GeomObject <tellipsoid>  *p2);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject <tcomposite>  *p1, const GeomObject <tcomposite>  *p2);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject<tellipsoid>   *p1, const GeomObject <tplane>   *plane);
	static void overlaps(vector<COverlapping> &ovs, const GeomObject <tcomposite>  *p1, const GeomObject <tbox>        *b );

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
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObjectBase *p1, const GeomObjectBase *p2){
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
		else ERROR(true, "Not Implemented");
		};
inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tsphere>  *p1, const GeomObject<tbox> *b){
	static vec v;
	static double d, dd;
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		v=b->face[i]->normal_from_point(p1->Xc, 0);// p1->radius);//vertical vector from the center of sphere to the plane
		d=v.abs();
		dd=p1->radius-d;
		//if(dd>0) ovs.push_back( COverlapping(p1->getpos()+v+(0.5*dd)*b->face[i]->n, (dd/d)*v) );
		if(dd>0) {
			ovs.push_back( COverlapping(p1->getpos()+v*(1+0.5*dd), (dd)*v) );
			}
		}
	}


inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tplane> *plane){
	static vec v, vp;
	if(plane->normal_from_point(p1->Xc).abs()-p1->radius > 0) return;
	vp=p1->point_to_plane(*(plane));
	v=plane->normal_from_point(vp, 0);
	if(v*plane->n <0)return;
	ovs.push_back( COverlapping(vp, -v) );
	return;
	}

inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tbox> *b){
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		overlaps(ovs, p1, b->face[i]);
		}
	}

inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tcomposite>  * p2){

	ERROR(p1==p2, "A particle is checked against itself for overlapping.")
	if((p1->Xc-p2->Xc).abs() > p1->radius+p2->radius)return;

	for(indexType i=0; i< (p1->elems.size()); ++i){
		for(indexType j=0; j< (p2->elems.size()); ++j){
			overlaps(ovs, p1->elems.at(i), p2->elems.at(j));
			}
		}
	}

inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tcomposite>  *p1, const GeomObject<tbox>  * b){
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

CPlane separatingPlane(const CEllipsoid  &E1, const CEllipsoid  E2){
	Matrix M=(-(!E1.ellip_mat)*E2.ellip_mat);
	//CQuartic q=characteristicPolynom(M);
	static int ii=0;

	vector<complex<double> > eigenvals;
	vector<HomVec> eigenvecs;
	eigens(M, eigenvals, eigenvecs);

		cout<<ii++<<"  "<< (E1.Xc-E2.Xc).abs()-E1.c-E2.c<<"  ";
		for(int i=0; i<4; i++){
			cout<< eigenvals.at(i).real()<<"  "<<eigenvals.at(i).imag()<<"  ";
		}
		cout<<endl;
	if(eigenvals.size()==4){
		CRay<HomVec> ray(eigenvecs.at(3),eigenvecs.at(2));
		CQuadratic q1(intersect(ray, E1));
		CQuadratic q2(intersect(ray, E2));
		//the roots are sorted ascending
		HomVec X1= ray(q1.root(1).real()); //on the surface of E1
		HomVec X2= ray(q2.root(0).real()); //on the surface of E2
		vec n1=HomVec(E1.ellip_mat*X1).get3d();
		vec n2=HomVec(E2.ellip_mat*X2).get3d();
		CPlane p1(X1.get3d(),n1);
		CPlane p2(X2.get3d(),n2);
		//CSphere S((X1.get3d()+X2.get3d())*.5, 0.001);
		//S.print(*gout);
		//(*gout)<<endl;
		//cerr<< (X1-X2).abs() <<endl;

		return CPlane((X1.get3d()+X2.get3d())*.5, (n1-n2)*.5);
		}
	
	

	ERROR(1, "no separating plane");
	}

void append(vector<COverlapping> &v, vector<COverlapping> &v2){
	for(size_t i=0; i<v2.size(); i++){
		v.push_back(v2.at(i));	
		}

	return;
	}
#define XOR(p, q) ( ((p) || (q)) && !((p) && (q)) ) 
inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const CEllipsoid  *E, const CEllipsoid  *E0){
	static CPlane *p=NULL;
	static bool collide=false;
	static vector<COverlapping> ovtest1;
	static vector<COverlapping> ovtest2;
	ovtest1.clear();
	ovtest2.clear();

	if(p==NULL and !collide){
		p= new CPlane(separatingPlane(*E, *E0));
		}
	overlaps(ovtest1, E, p);
	overlaps(ovtest2, E0, p);
	if(ovtest1.size()>0 and ovtest2.size()>0){
		//append(ovs, ovtest1);
		throw 1;
		collide=true;
		p->print(cerr);
		cerr<<endl;
		}
	if(XOR(ovtest1.size()>0 , ovtest2.size()>0) and !collide){
		 *p= separatingPlane(*E, *E0);//update the plane
		}
	
	//ERROR(1, "stopped");
	//ERROR(1,"Stopped here");
	}
inline
void overlapsTemp(vector<COverlapping> &ovs, const CEllipsoid  *E, const CEllipsoid  *E0){
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
		ovs.push_back(COverlapping(E->getpos()+v, (dp)*v));
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
			//ovs.push_back(COverlapping(rc, R-rE));
			return;
			}
		lambda += dx;
		}   
//	cerr<< dx <<endl;

	  return ;

}
#endif /* OVERLAPPING_H */
