#include<iostream>
#include"../ellipsoid.h"
#include"../MersenneTwister.h"
#include"/home/reza/workstation/mysrc/CStat.h"

using namespace std;

vec randomVec(const vec &x1, const vec &x2){
	return vec( x1(0)+(x2(0)-x1(0))*rgen(), x1(1)+(x2(1)-x1(1))*rgen(), x1(0)+(x2(2)-x1(2))*rgen());
	}

vec randomVecSurfaceSphere(const vec &x, const double r){
	
	vec X( rgen.randNorm(0, 0.5),rgen.randNorm(0, 0.5), rgen.randNorm(0, 0.5));
	X.normalize();
	X*=r;
	X+=x;
	
	return X;
	}

long RNGSeed;
extern MTRand rgen;

void Initialize(){
	rgen.seed(RNGSeed);
	cerr<< "RNG Seed: "<< RNGSeed <<endl;
	//define_parameters();
	//config.parse("config");
	}

void Run(){

	for(double ee=0.01; ee<100; ee*=5){
	double r=1.0;
	double a =r/pow(ee,1./3.);
	double b =a;//*rgen();
	double c =ee*a;//*rgen();
	GeomObject<tellipsoid> E(vec(0,0,0), a,b,c);


	CStat stat;
	for(int i=0; i<10000000; i++){
	vec xs=randomVecSurfaceSphere(E.Xc, r);
	CRay<HomVec> ray(HomVec(0,0,0,1), HomVec(xs,1));
	CQuadratic q=intersect(ray, E);
	assert(q.root(1).imag()<1e-10);
	vec xe= (ray(q.root(1).real())).project();
	stat<<(xe-xs).abs();
	}
	double sign=1;
	if(ee<1)sign=-1;
	cout<<ee<<"  "<< sign*stat.mean()<<"  "<<stat.standard_deviation() <<endl;
	}

}

int main(int pi, char **params){
	if(pi==1)
		RNGSeed=0;
	else
		RNGSeed=313*atoi(params[1])+1;
	
	try {
	Initialize();
	Run();
	//Shutdown();
	return 0;
	} catch(CException e)
	{
	e.Report();
	return 1;
	}
return 0;
}
