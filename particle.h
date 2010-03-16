#ifndef PARTICLE_H
#define PARTICLE_H 
#include<sstream>
#include<fstream>
#include<vector>
#include<string>
#include<iomanip>
#include"common.h"
#include"dfreedom.h"
#include"shapes.h"

using namespace std;
class CProperty
	{
	public:
	CProperty():
	 	stiffness1(paramsDouble("stiffness1")), stiffness2(paramsDouble("stiffness2")), color(" 1"), density(paramsDouble("density")){
		};
	~CProperty(){}
	double stiffness1, stiffness2;
	double density;
	string color;
 	private:
	};

//#define particleType tcomposite
#define particleType tsphere
class CParticle : public GeomObject<particleType>{
	public:
	CParticle(const vec & _x0, double r):GeomObject<particleType>(_x0, r), id(-1), forces(vec(0.0)), frozen(false){
		mass=material.density*4.0/3.0*M_PI*r*r*r;
		x(0)=_x0;
		w(0)=0.0;
		identifier=-1;
		};

	double kEnergy(){
		return 0.5*mass*x(1).abs();
		}

	double get_mass()const{return mass;}

	void addforce(const vec force){
		static vec prev(0.0);
		forces+=force;
		avgforces=(force);
		//avgforces=(prev+force)/2.0;
		prev=force;
		};


	void calPos(double dt);
	void calVel(double dt);

	vec forces, avgforces;
	CProperty material;
	CDFreedom<6> x, x0;//TranslationalDFreedom;
	CDFreedom<6> w, w0;//Rotational;
	int id, identifier;
	bool frozen;
	protected:
	double mass;
 	private:
	CDFreedom<5> RotationalDFreedom;
	};

ostream &operator <<(ostream &out, const CParticle &p){
		out<< p.identifier<< "   ";
		p.print(out);
		out<< "  "<<p.material.color;
	}

double friction=1;

void CParticle::calPos(double dt){
	static const double c=1./6.0;
	x(0) += x(1)*dt + x(2)*(dt*dt*4.0*c) - x0(2)*(dt*dt*c);
	x(1) += x(2)*(dt*5.0*c) - x0(2)*(dt*c);

	static vec n(1); //FIXME this temporary
	static double t=0;
	t+=dt;

	rotate(vec(1.0), 10*dt);
	moveto(x(0));
	}

void CParticle::calVel(double dt){
	static const double c=1./6.0;
	x0(2)=x(2);
	x(2)=forces/mass;
	x(1)+= x(2)*(dt*2*c);
	}

#endif /* PARTICLE_H */
