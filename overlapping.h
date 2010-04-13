#ifndef OVERLAPPING_H
#define OVERLAPPING_H 
#include"shapes.h"


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
		else ERROR("Not Implemented");
		};
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
void COverlapping::overlaps(vector<COverlapping> &ovs, const GeomObject<tellipsoid>  *p1, const GeomObject<tbox> *b){
	static vec v, vp;
	static double d, dd2;
	for(int i=0; i<6; ++i){//FIXME to generalize Box to any polygon, 6 should be the number of faces
		if(b->face[i]->normal_to_point(p1->Xc).abs()-p1->radius > 0) continue;
		vp=p1->point_to_plane(*(b->face[i]));
		v=b->face[i]->normal_to_point(vp, 0);
		if(v*b->face[i]->n <0)continue;
		
		ovs.push_back( COverlapping(vp, -v) );
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


inline
void COverlapping::overlaps(vector<COverlapping> &ovs, const CEllipsoid  *E, const CEllipsoid  *E0){
	//here we should if E will touch E0 or not;

// ------------------------------------------
	static vec v;
	static double d, dd;
	v=E0->displacement(E);//Xc2-Xc1
	d=v.abs();
	dd=E0->radius+E->radius-d;//FIXME can be put in the base class too

	if(dd<0) {return;
		v*=((E->radius-dd/2.0)/d); //from center of p1 to contact point
		v.normalized();
		ovs.push_back(COverlapping(E->getpos()+v, (dd)*v));
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

	int itermax=10000;
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

		double threshold=1e-3;
		static Matrix E12(3,3);
		static vec rc;
		static vec rE;
		if(fabs(dx) < threshold){
			//CVector rc;
			E12=!(E->ellip_mat-(lambda*E0->ellip_mat));
			rc=R0-E12*(E->ellip_mat*R12);
			double s2=lambda*lambda*(R*(E0->ellip_mat*E12*E->ellip_mat*E12*E0->ellip_mat)*R);
			double s=sqrt(s2);
			if(s>1){
				cerr<< s <<endl;
				break;//not overlapping
				}
			rE=(1-1/s)*rc+(1/s)*R;
			//cerr<< s <<endl;
			
			//return CContact(rc-R, rE-R);
			cerr<< E->radius<<"   x "<<(rc-R).abs()<<"  dx "<<rE-R <<endl;
			ovs.push_back(COverlapping(rc, R-rE));
			return;
			}
		lambda -= dx;
		}   

	  return ;

}
#endif /* OVERLAPPING_H */
