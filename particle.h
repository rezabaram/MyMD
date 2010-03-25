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
		x(0)=_x0;
		w(0)=0.0;
		for(int i=1; i<6; ++i){
			x(i)=0.0;
			w(i)=0.0;
			}

		
		mass=material.density*4.0/3.0*M_PI*r*r*r;
		if(type==tcomposite){
			mass=mass*(1.5);//FIXME for a specific composite particle
			Iyy=2*(2*mass*r*r/5.0+mass*r*r/128.0);
			Ixx=2*mass*r*r/5.0+2*mass*r*r/32.0/5.0;
			Izz=Iyy;
			}
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

	void addtorque(const vec torque){
		torques+=torque;
		};

	void parse(std::istream &in){
			GeomObject<particleType>::parse(in);
			x(0)=Xc;
			mass=material.density*4.0/3.0*M_PI*radius*radius*radius;
			}

	void calPos(double dt);
	void calVel(double dt);

	vec forces, avgforces;
	vec torques, avgtorque;
	CProperty material;
	CDFreedom<6> x, x0;//TranslationalDFreedom;
	CDFreedom<6> w, w0;//Rotational;
	int id;
	bool frozen;
	vec test;
	protected:
	double mass, Ixx, Iyy, Izz;
 	private:
	CDFreedom<5> RotationalDFreedom;
	};

ostream &operator <<(ostream &out, const CParticle &p){
	p.print(out);
	//out<<"  "<<p.x(1);
	out<<"  "<<drand48()<<" "<<drand48()<<" "<<drand48();
	}

double friction=1;

//using beeman method
void CParticle::calPos(double dt){
	static const double c=1./6.0;
	//translational degree
	x(0) += x(1)*dt + x(2)*(dt*dt*4.0*c) - x0(2)*(dt*dt*c);
	x(1) += x(2)*(dt*5.0*c) - x0(2)*(dt*c);


	//rotational degree
	static vec wp;
	static Quaternion dq;

	wp=q.toBody(w(1));
	w(1) += w(2)*(dt*5.0*c) - w0(2)*(dt*c);
	dq.u =    -q.v(0)*w(1)(0) - q.v(1)*w(1)(1) - q.v(2)*w(1)(2);
	dq.v(0) =  q.u  * w(1)(0) - q.v(2)*w(1)(1) + q.v(1)*w(1)(2);
	dq.v(1) =  q.v(2)*w(1)(0) + q.u *  w(1)(1) - q.v(0)*w(1)(2);
	dq.v(2) = -q.v(1)*w(1)(0) + q.v(0)*w(1)(1) + q.u *  w(1)(2);

	q+=dq*dt*0.5;
	//cout<< q.abs() <<endl;
	q.normalize();
	
	rotateTo(q);
	moveto(x(0));
	}

void CParticle::calVel(double dt){
	static const double c=1./6.0;
	x0(2)=x(2);
	x(2)=forces/mass;
	x(1)+= x(2)*(dt*2*c);
	
	w0(2)=w(2);
	static vec wp, wwp, torquep;
	torquep=q.toBody(torques);
	wp=q.toBody(w(1));
	wwp(0)=(torquep(0)+w(1)(1)*w(1)(2)*(Iyy-Izz))/Ixx;
	wwp(1)=(torquep(1)+w(1)(0)*w(1)(2)*(Izz-Ixx))/Iyy;
	wwp(2)=(torquep(2)+w(1)(0)*w(1)(1)*(Ixx-Iyy))/Izz;
	w(2)=q.toWorld(wwp);

	w(1)+=w(2)*(dt*2*c);
	}


#endif /* PARTICLE_H */
