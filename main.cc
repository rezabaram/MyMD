#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;


int main(){
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



double time=0;
vec x(0.0, 0.0, .0);

GeomObject<tellipsoid> E(vec(0.5, 0.5, 0.3), 1, 1, 0.5);
E.scale(0.2);
CParticle *p = new CParticle(E);
//p->q=Quaternion(cos(M_PI/7.),sin(M_PI/7.),0,0 )*Quaternion(cos(M_PI/15.),0,0,sin(M_PI/15.) );
//sys.add(p);

GeomObject<tellipsoid> E2(vec(0.5, 0.5, 0.7), 1, 1, 2);
E2.scale(0.2);
CParticle *p2 = new CParticle(E2);
p->q=Quaternion(cos(M_PI/18.),sin(M_PI/18.),0,0 )*Quaternion(cos(M_PI/13.),0,0,sin(M_PI/13.) );
//E2.rotateTo(p->q);
//E2.moveto(vec(1, 2,4.35));
Matrix M=(-(!E2.ellip_mat)*E.ellip_mat);
cerr<< M <<endl;
cerr<< "Det= "<<M.Det() <<endl;
CQuartic quart=characteristicPolynom(M);
//CQuartic quart=CQuartic(1., - 468.562 ,  - 939.125, - 468.563  ,1);
quart.print(cerr);
quart.plot(cout,-3, 500, 0.01);
quart.print_roots(cerr);
//sys.add(p2);

//sys.solve(config.get_param<double>("maxTime"), Dt);
return 0;

double size=config.get_param<double>("particleSize");

double margin=2.0*size;
for(double i=margin; i<1-margin; i+=margin){
for(double j=1-margin; j>margin; j-=margin){
for(double k=1-margin; k>margin; k-=margin){

	if(sys.particles.size()==config.get_param<int>("nParticle")) {
				break;
				}
	x(1)=0.5;//i+size*drand48()/10; 
	x(0)=j+size*drand48()/10;
	x(2)=k+size*drand48()/10; 
	//CParticle *p = new CParticle(GeomObject<tsphere>(x,size*(1+0.2*drand48())));
	GeomObject<tellipsoid> E(x, 1, 0.7, 0.5);
	E.scale(size*(1+0.2*drand48()));
	CParticle *p = new CParticle(E);
	vec axis(0.0);
	axis(1)=30.0*drand48();
	//axis(0)=axis(2);
	//axis(0)=2.0*drand48();
	p->q=Quaternion(cos(M_PI/8.),sin(M_PI/8.),0,0 )*Quaternion(cos(M_PI/8.),0,0,sin(M_PI/8.) );
	//p->q=Quaternion(cos(M_PI/22.0),0,sin(M_PI/22.0),0 );
	//p->w(1)=axis;
	//p->x(1)(0)=2;
	//p->rotate(vec(1.0), drand48()*3.14);
	sys.add(p);
	}
	}
	}
	sys.solve(config.get_param<double>("maxTime"), Dt);

return 0;
}
