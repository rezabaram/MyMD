#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;


int main(){
define_parameters();
config.parse("config");

double size=config.get_param<double>("particleSize");

CSys sys;
sys.outDt=config.get_param<double>("outDt");
double Dt=config.get_param<double>("timeStep");


double shift[3]={0};
//cerr<< sys.read_packing("Mini9.dat",vec3d(shift), 1./2100.0)<<endl;
//cerr<< sys.read_packing("Mini9.dat")<<endl;
//sys.overlappings();

sys.G=config.get_param<vec>("Gravity");
//sys.read_packing3("coord2.dat");
//sys.write_packing("test.dat");
//sys.solve(2, 0.0001);
//return 0;

double time=0;
vec x;
double margin=3.1*size;
for(double i=margin; i<1-3.1*margin; i+=3.1*size){
for(double j=1-3.1*size; j>margin; j-=3.1*size){
for(double k=1-3.1*size; k>margin; k-=3.1*size){
	x(1)=i+size*drand48()/10; 
	x(0)= 0.5+size*drand48()/10;
	x(2)=k+size*drand48()/10; 

	CParticle *p = new CParticle(x,size*(1+0.2*drand48()));
	double ran=drand48();
	if(ran<0.1)p->material.color="1 0 0";
	else if(ran<0.45)p->material.color="0 1 0";
	else p->material.color="0 0 1";
	double col=drand48();
	p->identifier=1;
	vec axis(0.0);
	axis(1)=20.0*drand48();
	p->w(1)=axis;
	//axis(0)=2.0*drand48();
	p->q=Quaternion(cos(M_PI/8.),sin(M_PI/8.),0,0 )*Quaternion(cos(M_PI/8.),0,0,sin(M_PI/8.) );
	//p->q=Quaternion(cos(M_PI/8.),0,sin(M_PI/8.),0 );
//	p->x(1)=axis;
	//p->rotate(vec(1.0), drand48()*3.14);
	if(sys.particles.size()<config.get_param<int>("nParticle")) sys.add(p);
	}
	}
	}
	sys.solve(config.get_param<double>("maxTime"), Dt);

return 0;
}
