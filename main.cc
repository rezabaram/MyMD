#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;

int main(){
define_parameters();



double v[]={0,0.0, -.0};

CSys sys;
sys.outDt=0.02;

double shift[3]={0};
//cerr<< sys.read_packing("Mini9.dat",vec3d(shift), 1./2100.0)<<endl;
//cerr<< sys.read_packing("Mini9.dat")<<endl;
//sys.overlappings();

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
for(double xx=1.0; xx < 7.2; xx+=.20){
	//p.set_size(0.03+0.01*(1.0-2*drand48()));
	x(0)=(0.3+0.2*(drand48()));
	x(1)=0.5;//(0.08+0.8*(drand48()));
	x(2)=0.5;
	//x*=0.5; //FIXME some where it is multipled by two
	CParticle *p = new CParticle(x,0.07);
	double col=drand48();
	//if(col<0.6)p.material.color=" 1";
	//else if(col<0.8)p.material.color=" 2";
	//else  p.material.color=" 3";
	sys.add(p);
	sys.solve(xx, 0.0001);
	}

return 0;
}
