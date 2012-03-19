// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef ELLIPS_CONTACT_H
#define ELLIPS_CONTACT_H 
#include "exception.h"
#include "eigen.h"
#include"multicontact.h"
#include"ellipsoid.h"

void fixcontact(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2)
	{
TRY
	ovs.x01=E1.toBody(ovs.x1);
	ovs.x02=E2.toBody(ovs.x2);
CATCH
	}

void updatecontact(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2)
	{
TRY
	ovs.x1=E1.toWorld(ovs.x01);
	ovs.x2=E2.toWorld(ovs.x02);
CATCH
	}

void updateplane(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2)
	{
TRY
	ovs.plane.Xc=(ovs.x1.project()+ovs.x2.project())/2;
	ovs.plane.n=(E1.gradient(ovs.x1.project())-E2.gradient(ovs.x2.project())).normalized();
CATCH
	}

void correctpoints(ShapeContact &ovs, const CEllipsoid &E1, CEllipsoid &E2)
	{
TRY
	ovs.x1=HomVec(E1.point_to_plane((ovs.plane)),1);
	ovs.x2=HomVec(E2.point_to_plane((ovs.plane)),1);
	ovs.x01=E1.toBody(ovs.x1);
	ovs.x02=E2.toBody(ovs.x2);
	//static int i=0;
	//if( E2(ovs.x1) > epsilon or E1(ovs.x2) >epsilon ) cout<<"  "<< ++i<<"  "<<E2(ovs.x1) <<" "<<E1(HomVec(ovs.plane.Xc,1))<<" "<<E2(HomVec(ovs.plane.Xc,1))<<"  "<<E1(ovs.x2)<<endl;
CATCH
	}


//this function takes two intersection ellipsoids, and gives back the intersection of a 
//ray through midpoint mp (which should be inside both of them) in direction of the average gradient 
// at point mp. the intersecting points are X1 on E1 and X2 on E2
void intersect(HomVec &X1, HomVec &X2,  const CEllipsoid &E1, const CEllipsoid &E2)
	{
TRY
	HomVec mp=(X1+X2)/2;
	if(E1(mp)>1e-12 or E2(mp) > 1e-12)
		{
		WARNING("contact point is not inside both ellipsoids: "<<E1(mp)<<"   "<<E2(mp));
		return;
		}
	
	HomVec g1=HomVec(E1.gradient(mp.project() ),0);
	HomVec g2=HomVec(E2.gradient(mp.project() ), 0);
	HomVec g=g2-g1; g.normalize();
	//if(g1*g2<0)g*=-1;

	CRay<HomVec> ray(mp, mp+g);
	CQuadratic q1(intersect(ray, E1));
	CQuadratic q2(intersect(ray, E2));

	ERROR(fabs(q1.root(0).imag()) > epsilon, "the intersection of line with ellipsoid is complex.");
	ERROR(fabs(q2.root(1).imag()) > epsilon, "the intersection of line with ellipsoid is complex.");

	X1= ray(q1.root(0).real());//on the surface of E1
	X2= ray(q2.root(1).real());//on the surface of E2

CATCH
	}

void setcontact(ShapeContact &ovs,CEllipsoid  &E1, CEllipsoid  &E2)
	{
TRY
	//ovs.add(Contact(ovs.plane.Xc, ovs.plane.n, fabs((ovs.x1.project()-ovs.x2.project())*ovs.plane.n)));

	vec diff=(ovs.x1.project()-ovs.x2.project());
	vec mp=((ovs.x1.project()+ovs.x2.project())/2.0);
	vec g1=E1.gradient(mp);
	vec g2=E2.gradient(mp);
	double dx=diff.abs();
	diff.normalize();

	//if the contact points are not along the normal direction, correct them
	if(fabs(g1*diff/g1.abs()) <0.999 or fabs(g2*diff/g1.abs())<0.999 )
		{
		intersect(ovs.x1, ovs.x2, E1, E2);
		diff=(ovs.x1.project()-ovs.x2.project());
		mp=((ovs.x1.project()+ovs.x2.project())/2.0);
		dx=diff.abs();
		diff.normalize();
		}

	ovs.add(Contact(mp, diff, dx));

	ovs.x01=E1.toBody(ovs.x1);
	ovs.x02=E2.toBody(ovs.x2);
	//ovs.add(Contact((ovs.x1.project()+ovs.x2.project())/2, (ovs.x1.project()-ovs.x2.project()).normalized(), (ovs.x1-ovs.x2).abs()));
CATCH
	}

/*
void charpolynom(const CEllipsoid &A, const CEllipsoid &B){

	double a=1/A.a/A.a;
	double b=1/A.b/A.b;
	double c=1/A.c/A.c;
	
	u=

	}
*/

bool doOverlap(ShapeContact &ovs,  CEllipsoid  &E1, CEllipsoid  &E2){
TRY

	if(0)if(ovs.has_sep_plane){
		if(!(E1.doesHit(ovs.plane) or E2.doesHit(ovs.plane))) {
			return false;
			}
		}
	
	Matrix M=(-(!E1.ellip_mat)*E2.ellip_mat);
	//CQuartic q=characteristicPolynom(M);

	vector<complex<double> > eigenvals;
	vector<HomVec> eigenvecs;
	eigens(M, eigenvals, eigenvecs);


	if(fabs(eigenvals.at(3).imag() ) < epsilon){
		ERROR(eigenvals.at(2).imag()>epsilon,"one eigenvalue complex one not.");

		// eigenvec 3 is inside E1, and 2 inside E2
		if(E1(eigenvecs.at(3))>0 or E2(eigenvecs.at(2))>0){
			WARNING("error in calculation of poles.");
			return false;
			}

		//line through two poles (from E1 to E2)
		CRay<HomVec> ray(eigenvecs.at(3).project4d(),eigenvecs.at(2).project4d());
		//CRay<HomVec> ray(HomVec(E1.Xc,1),HomVec(E2.Xc,1));
		
		CQuadratic q1(intersect(ray, E1));
		CQuadratic q2(intersect(ray, E2));
		//the roots are sorted ascending
		ERROR(fabs(q1.root(0).imag()) > epsilon, "the intersection of line with ellipsoid is complex.");
		ERROR(fabs(q2.root(1).imag()) > epsilon, "the intersection of line with ellipsoid is complex.");

		HomVec X1= ray(q1.root(1).real());//on the surface of E1
		HomVec X2= ray(q2.root(0).real());//on the surface of E2

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

		//cerr<< ovs.plane.n*(E2.Xc-E1.Xc).normalize() <<endl;
		//cerr<< (X2.project()-X1.project()).normalize()*(E2.Xc-E1.Xc).normalize() <<endl;
		//cerr <<endl;

		//FIXME there are some problems with calculation of separation plane. there for
		//it is deactivated for now
		/*
		bool hit=(E1.doesHit(ovs.plane) or E2.doesHit(ovs.plane));
		if(hit){
			WARNING("The separation plane hits the ellipsoids: E1(plane.x)="<< E1(ovs.plane.Xc) <<"\t E2(plane.x)="<< E2(ovs.plane.Xc) );
			return false;
			}
		*/

		return false ;
		}
	return true;
CATCH
	}

// find min of x on E1, in the potentional of E2 (iteretively)
// for fast convergence x should be initially on E1 and near to minimum 
void findMin(HomVec &x,  CEllipsoid  &E1, CEllipsoid  &E2, long nIter=1){
TRY
	double lambda, lambda0;
	vec xp0, xp=x.project();
	static Matrix Em1(3,3), Em2(3,3);
	for(size_t i=0; i<3; i++){
	for(size_t j=0; j<3; j++){
		Em1(i,j)=E1.ellip_mat(i,j);
		Em2(i,j)=E2.ellip_mat(i,j);
		}
		}
	
	long iter=0;
	bool converged=false;
	lambda=fabs((xp-E1.Xc)*E2.ellip_mat*(xp-E2.Xc));
	do{
		++iter;
		xp0=xp;
		lambda0=lambda;
		//using the fact that the gradients of the potentials of the 
		//the ellipsoids are in opposite directions at the minimum
		lambda=fabs((xp-E1.Xc)*E2.ellip_mat*(xp-E2.Xc));
		xp=(!(Em2+lambda*Em1))*(Em2*E2.Xc+lambda*Em1*E1.Xc);
		converged= (xp-xp0).abs()<1e-13 and fabs(lambda0-lambda)<1e-10;
		}
	while(iter<nIter and !converged);
	//	cout<<setprecision(14)<<iter<<"  lambda="<< lambda<<"   Xp="<<xp <<"  C1="<< E1.Xc<<"  C2="<<E2.Xc<<endl;
	

	if(!converged)WARNING("minimization not converged: "<<(xp-xp0).abs()<<"    "<<fabs(lambda0-lambda));
	if(converged and E2(xp) > 0)WARNING("A minimum point is not inside the corresponding ellipse: "<<xp<<". E(x)= "<<E2(xp)<<endl<<E1<<endl<<E2);
	x(0)=xp(0);
	x(1)=xp(1);
	x(2)=xp(2);
CATCH
}
#endif /* ELLIPS_CONTACT_H */
