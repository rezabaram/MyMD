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

double v[]={0,0.0, -.0};
vec G(0.0);
G(2)=-10;
sys.G=G;
//sys.read_packing3("coord2.dat");
//sys.write_packing("test.dat");
//sys.solve(2, 0.0001);
//return 0;

double time=0;
vec x;
double margin=2*size;
for(double i=margin; i<1-margin; i+=4.3*size){
for(double j=margin; j<1-margin; j+=4.3*size){
for(double k=margin; k<1-2*size; k+=4.3*size){
	x(0)=i+size*drand48()/10; 
	x(1)=j+size*drand48()/10;
	x(2)=k+size*drand48()/10; 
	CParticle *p = new CParticle(x,size);
	double ran=drand48();
	if(ran<0.1)p->material.color="1 0 0";
	else if(ran<0.45)p->material.color="0 1 0";
	else p->material.color="0 0 1";
	double col=drand48();
	p->identifier=1;
	p->rotate(vec(1.0), drand48()*3.14);
	sys.add(p);
	if(sys.particles.size()>15)break;
	}
	}
	}
	sys.solve(2, Dt);

return 0;
}
