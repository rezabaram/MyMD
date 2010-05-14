#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;


void Run(){
define_parameters();
config.parse("config");


CSys sys;
sys.outDt=config.get_param<double>("outDt");
double Dt=config.get_param<double>("timeStep");


sys.G=config.get_param<vec>("Gravity");
//sys.read_packing3("coord2.dat");
//sys.write_packing("test.dat");
//sys.solve(2, 0.0001);
//return 0;


double size=config.get_param<double>("particleSize");

double margin=3.0*size;
vec x(0.0, 0.0, .0);
for(double i=margin; i<1-margin; i+=margin){
for(double j=1-margin; j>margin; j-=margin){
for(double k=1-margin; k>margin; k-=margin){

	if(sys.particles.size()==config.get_param<size_t>("nParticle")) {
		break;
		}
	x(1)=0.5;//i+size*drand48()/10; 
	x(0)=j+size*drand48()/10;
	x(2)=k+size*drand48()/10; 
	//CParticle *p = new CParticle(GeomObject<tsphere>(x,size*(1+0.2*drand48())));
	GeomObject<tellipsoid> E(x, 1, 1-0.4*drand48(), 0.6);
	E.scale(size*(1+0.2*drand48()));
	CParticle *p = new CParticle(E);
	vec axis(0.0);
	axis(1)=30.0*drand48();
	//axis(0)=axis(2);
	//axis(0)=2.0*drand48();
	//p->q=Quaternion(cos(M_PI/8.),sin(M_PI/8.),0,0 )*Quaternion(cos(M_PI/8.),0,0,sin(M_PI/8.) );
	//p->w(1)=axis;
	p->x(1)(2)=1;
	//p->rotate(vec(1.0), drand48()*3.14);
	sys.add(p);
	}
	}
	}
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
