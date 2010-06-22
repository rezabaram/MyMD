#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;

void particles_on_grid(CSys &sys){ 

double size=sys.maxr;
double margin=2.3*size;
vec x(0.0, 0.0, .0);

double k=0;
double t=0;
double Dt=config.get_param<double>("timeStep");
while(sys.particles.size()<sys.maxNParticle){
	k=sys.maxh+margin;
	if(k==1+10*margin)break;
for(double i=margin/2; i<1-margin/2; i+=margin){
for(double j=1-margin/2; j>margin/2; j-=margin){
	if(sys.particles.size()==sys.maxNParticle)break;

	x(1)=i+size*drand48()/10; 
	x(0)=j+size*drand48()/10;
	x(2)=k+size*drand48()/10; 
	double alpha=drand48()*M_PI;
	Quaternion q=Quaternion(cos(alpha),sin(alpha),0,0)*Quaternion(cos(alpha),0,0,sin(alpha) );
	//CParticle *p = new CParticle(GeomObject<tsphere>(x,size*(1-0.0*drand48())));
	//GeomObject<tellipsoid> E(x, 1-0.0*drand48(), 1-0.0*drand48(),1-0.0*drand48(), size*(1+0.0*drand48()));
	//CParticle *p = new CParticle(E);
	double r=size*(1-0.1*drand48());
	GeomObject<tsphere> E1(x,r);
	//GeomObject<tellipsoid> E2(x, 1, 1, 1, size, q);

	double ee=0.0;
	double a =1;
	double b =1-ee*drand48();
	double c =1-ee*drand48();
	GeomObject<tellipsoid> E2(x, a,b,c, r);
	CParticle *p = new CParticle(E2);
	//p->x(1)(0)=1-2*drand48();
	//p->x(1)(2)=1-2*drand48();
	p->w(1)(2)=10*(1-2*drand48());
	sys.add(p);
	
	}
	}
	t+=0.2;
	if(t>config.get_param<double>("maxTime"))exit(0);
	sys.setup_verlet();
	sys.solve(t, Dt);
	}
}

long RNGSeed;
void Initialize(){
	RNGSeed=0;
	define_parameters();
	config.parse("config");


	}

void Run(){

	CSys sys(config.get_param<size_t>("nParticle"));
	sys.outDt=config.get_param<double>("outDt");
	sys.maxr=config.get_param<double>("particleSize");

	sys.G=config.get_param<vec>("Gravity");
	particles_on_grid(sys);
	cerr<< "Number of Particles: "<<sys.particles.size() <<endl;
	double Dt=config.get_param<double>("timeStep");
	sys.solve(config.get_param<double>("maxTime"), Dt);
}

int main(){
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
