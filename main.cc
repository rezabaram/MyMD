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

vec x;
double margin=8*size;
for(double i=margin; i<1-margin; i+=3*size)
for(double j=margin; j<1-margin; j+=3*size)
for(double k=margin; k<1-2.0*size; k+=3*size){
	x(0)=i+size*drand48()/2; x(1)=j+size*drand48()/2; x(2)=k+size*drand48()/2;
	CParticle *p = new CParticle(x,size);
	double col=drand48();
	sys.add(p);
	}
	sys.solve(3, 0.0001);

return 0;
}
