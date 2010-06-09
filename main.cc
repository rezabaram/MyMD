#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;

void particles_on_grid(CSys &sys){ 

double size=sys.maxr;
double margin=2.3*size;
vec x(0.0, 0.0, .0);

for(double i=margin/2; i<1-margin/2; i+=margin){
for(double j=1-margin/2; j>margin/2; j-=margin){
for(double k=1-margin/2; k>margin/2; k-=margin){

	if(sys.particles.size()==sys.maxNParticle) {
		break;
		}
	x(1)=i+size*drand48()/10; 
	x(0)=j+size*drand48()/10;
	x(2)=k+size*drand48()/10; 
	double alpha=drand48()*M_PI;
	Quaternion q=Quaternion(cos(alpha),sin(alpha),0,0)*Quaternion(cos(alpha),0,0,sin(alpha) );
	CParticle *p = new CParticle(GeomObject<tsphere>(x,size*(1-0.0*drand48())));
	//GeomObject<tellipsoid> E(x, 1-0.0*drand48(), 1-0.0*drand48(),1-0.0*drand48(), size*(1+0.0*drand48()));
	//CParticle *p = new CParticle(E);
	sys.add(p);
	}
	}
	}
}

void Run(){
	define_parameters();
	config.parse("config");


	CSys sys(config.get_param<size_t>("nParticle"));
	sys.outDt=config.get_param<double>("outDt");
	double Dt=config.get_param<double>("timeStep");
	sys.maxr=config.get_param<double>("particleSize");

	sys.G=config.get_param<vec>("Gravity");
	particles_on_grid(sys);

	cerr<< "Number of Particles: "<<sys.particles.size() <<endl;
	sys.solve(config.get_param<double>("maxTime"), Dt);
}

int main(){
	try {
	//Initialize();
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
