#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;


void Run(){
define_parameters();
config.parse("config");


CSys sys(config.get_param<size_t>("nParticle"));
sys.outDt=config.get_param<double>("outDt");
double Dt=config.get_param<double>("timeStep");


sys.G=config.get_param<vec>("Gravity");
//sys.read_packing3("coord2.dat");
//sys.write_packing("test.dat");
//sys.solve(2, 0.0001);
//return 0;


double size=config.get_param<double>("particleSize");

double margin=2.3*size;
vec x(0.0, 0.0, .0);

double ratio=0.3;
for(double i=margin/2; i<1-margin/2; i+=margin){
for(double j=1-margin/2; j>margin/2; j-=margin){
for(double k=1-margin/2; k>margin/2; k-=margin){

	if(sys.particles.size()==sys.maxNParticle) {
		break;
		}
	x(1)=i+size*drand48()/10; 
	x(0)=j+size*drand48()/10;
	x(2)=k+size*drand48()/10; 
	//CParticle *p = new CParticle(GeomObject<tsphere>(x,size*(1+0.2*drand48())));
	double alpha=drand48()*M_PI;
	Quaternion q=Quaternion(cos(alpha),sin(alpha),0,0)*Quaternion(cos(alpha),0,0,sin(alpha) );
	GeomObject<tellipsoid> E(x, 1-0.5*drand48(), 1-0.5*drand48(),1-0.5*drand48(), size*(1+0.0*drand48()));
	ratio=ratio*-1.0;
	//GeomObject<tellipsoid> E(x, 0.7+ratio, 0.7-ratio,0.7-ratio, size*(1+0.0*drand48()), q);
	CParticle *p = new CParticle(E);
	vec axis(0.0);
	axis(1)=30.0*drand48();
	//axis(0)=axis(2);
	//axis(0)=2.0*drand48();
	//p->w(1)=axis;
	//p->x(1)(2)=1;
	//p->rotate(vec(1.0), drand48()*3.14);
	sys.add(p);
	}
	}
	}
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
