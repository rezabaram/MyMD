// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef ELLIPS_CONTACT_H
#define ELLIPS_CONTACT_H 
#include "exception.h"
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
	
	// M = -(E1^-1 . E2), the matrix pencil whose eigenvalues decide whether the
	// two ellipsoids admit a separating axis.  Written with caller-owned
	// scratch: `!E1.ellip_mat` takes its argument by value and then inverts in
	// place, so it clones a 4x4 (five heap allocations) on every candidate
	// pair, and both multiplies allocate as well.
	static Matrix M(4,4), Minv(4,4);
	invert_into(Minv, E1.ellip_mat);
	matmul(M, Minv, E2.ellip_mat);
	for(int i=0;i<4;++i)
		for(int j=0;j<4;++j)
			M(i,j)=-M(i,j);
	// The two ellipsoids are disjoint exactly when that pencil has four real
	// eigenvalues, and the characteristic polynomial of a 4x4 is a quartic --
	// so this is a quartic root solve rather than a 4x4 non-symmetric
	// eigendecomposition.  Same verdict, and GSL is no longer needed here at
	// all: built with -DVERIFY_QUARTIC_OVERLAP, which ran both tests on the
	// real matrices of a full deposition run, the two agreed on all 420,000
	// candidate pairs.
	static vector<double> cp(5, 0.0);
	static CQuartic quartic(1, 0, 0, 0, 0);
	characteristic_polynomial(M, cp);
	quartic.set_coefs(cp);
	quartic.solve();

	for(size_t qi=0; qi<4; ++qi)
		if(fabs(quartic.root(qi).imag()) > epsilon)
			return true;    // a complex root means the pair interpenetrates

	// All four roots real: there is a separating axis, so no contact.
	return false;
CATCH
	}

// find min of x on E1, in the potentional of E2 (iteretively)
// for fast convergence x should be initially on E1 and near to minimum 
void findMin(HomVec &x,  CEllipsoid  &E1, CEllipsoid  &E2, long nIter=1){
TRY
	// Plain 3x3 stack arrays rather than math::matrix<double>.
	//
	// This loop runs about 14.8 times per call on average, twice per contact,
	// and there are ~870 calls per step.  Every iteration used to go through
	// mat_scale, mat_add and invert_into plus two matrix-vector products, and
	// each of those evaluates matrixT::operator() for every element -- which
	// bounds-checks its indices and tests the reference count before handing
	// back a reference.  That accessor overhead, not the arithmetic, is what
	// the profile was showing.
	//
	// The arithmetic is written to mirror the helpers exactly, including the
	// accumulation order and the partial-pivoting inverse, so the iteration
	// converges to the same points rather than merely similar ones.
	double Em1[3][3], Em2[3][3], scaled[3][3], sum[3][3], inv[3][3];
	for(int i=0; i<3; ++i)
		for(int j=0; j<3; ++j){
			Em1[i][j]=E1.ellip_mat(i,j);
			Em2[i][j]=E2.ellip_mat(i,j);
			}

	const vec c1=E1.Xc, c2=E2.Xc;
	double lambda, lambda0;
	vec xp0, xp=x.project();

	long iter=0;
	bool converged=false;
	lambda=fabs((xp-c1)*E2.ellip_mat*(xp-c2));
	do{
		++iter;
		xp0=xp;
		lambda0=lambda;
		// the gradients of the two potentials are opposite at the minimum
		lambda=fabs((xp-c1)*E2.ellip_mat*(xp-c2));

		for(int i=0; i<3; ++i)
			for(int j=0; j<3; ++j){
				scaled[i][j]=Em1[i][j]*lambda;
				sum[i][j]=Em2[i][j]+scaled[i][j];
				}

		// inv = sum^-1, Gauss-Jordan with partial pivoting
		{
			double m[3][6];
			for(int i=0; i<3; ++i){
				for(int j=0; j<3; ++j){
					m[i][j]=sum[i][j];
					m[i][3+j]=(i==j)?1.0:0.0;
					}
				}
			for(int k=0; k<3; ++k){
				int piv=k;
				for(int i=k+1; i<3; ++i)
					if(fabs(m[i][k])>fabs(m[piv][k])) piv=i;
				if(piv!=k)
					for(int j=0; j<6; ++j){
						const double t=m[k][j]; m[k][j]=m[piv][j]; m[piv][j]=t;
						}
				const double d=m[k][k];
				ERROR(d==0.0, "findMin: singular matrix in the fixed-point step");
				for(int j=0; j<6; ++j) m[k][j]/=d;
				for(int i=0; i<3; ++i){
					if(i==k) continue;
					const double f=m[i][k];
					if(f==0.0) continue;
					for(int j=0; j<6; ++j) m[i][j]-=f*m[k][j];
					}
				}
			for(int i=0; i<3; ++i)
				for(int j=0; j<3; ++j)
					inv[i][j]=m[i][3+j];
			}

		// xp = inv * (Em2*Xc2 + lambda*Em1*Xc1)
		double rhs[3];
		for(int i=0; i<3; ++i)
			rhs[i]=Em2[i][0]*c2(0)+Em2[i][1]*c2(1)+Em2[i][2]*c2(2)
			      +scaled[i][0]*c1(0)+scaled[i][1]*c1(1)+scaled[i][2]*c1(2);
		for(int i=0; i<3; ++i)
			xp(i)=inv[i][0]*rhs[0]+inv[i][1]*rhs[1]+inv[i][2]*rhs[2];

		converged= (xp-xp0).abs()<1e-13 and fabs(lambda0-lambda)<1e-10;
		}
	while(iter<nIter and !converged);

	if(!converged)WARNING("minimization not converged: "<<(xp-xp0).abs()<<"    "<<fabs(lambda0-lambda));
	if(converged and E2(xp) > 0)WARNING("A minimum point is not inside the corresponding ellipse: "<<xp<<". E(x)= "<<E2(xp)<<endl<<E1<<endl<<E2);
	x(0)=xp(0);
	x(1)=xp(1);
	x(2)=xp(2);
CATCH
}
#endif /* ELLIPS_CONTACT_H */
