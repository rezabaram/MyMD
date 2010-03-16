#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;

int main(){
define_parameters();
config.parse("config");

double size=config.get_param<double>("particleSize");
cerr<< paramsDouble("stiffness2") <<endl;

CSys sys;
sys.outDt=config.get_param<double>("outDt");


double shift[3]={0};
//cerr<< sys.read_packing("Mini9.dat",vec3d(shift), 1./2100.0)<<endl;
//cerr<< sys.read_packing("Mini9.dat")<<endl;
//sys.overlappings();

double v[]={0,0.0, -.0};
vec G(0.0);
sys.G=G;
//sys.solve(1, 0.005);
//sys.write_packing("coord1.dat");
G(2)=-10;
sys.G=G;
//sys.solve(3., 0.005);
//sys.write_packing("coord2.dat");
G(2)=-10;
sys.G=G;
//sys.solve(5., 0.005);
//sys.write_packing("coord3.dat");
//return 0;
sys.read_packing2("coord2.dat");
sys.write_packing("coord3.dat");
return 0;
double time=0;
vec x;
while(sys.particles.size()<1200){
x(0)=drand48()*0.9; x(1)=0.9*drand48(); x(2)=0.9;
CParticle *p = new CParticle(x,size);
time+=0.03;
double ran=drand48();
if(ran<0.1)p->material.color="1 0 0";
else if(ran<0.2)p->material.color="0 1 0";
else p->material.color="0 0 1";
sys.add(p);
sys.solve(time, 0.0001);
}
sys.solve(100, 0.0001);
//sys.write_packing("coord3.dat");
return 0;
double margin=2*size;
for(double i=margin; i<1-margin; i+=2.3*size)
for(double j=margin; j<1-margin; j+=2.3*size)
for(double k=margin; k<1-2.0*size; k+=2.3*size){
	x(0)=i+size*drand48()/10; x(1)=j+size*drand48()/10; x(2)=k+size*drand48()/10;
	CParticle *p = new CParticle(x,size);
	double ran=drand48();
	if(ran<0.1)p->material.color="1 0 0";
	else if(ran<0.45)p->material.color="0 1 0";
	else p->material.color="0 0 1";
	double col=drand48();
	sys.add(p);
	}
	sys.solve(3, 0.0001);

return 0;
}
