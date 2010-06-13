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


typedef enum {frozen, onhold, rejected, ready_to_go} tState;

using namespace std;
class CProperty
	{
	public:
	CProperty():
	 	stiffness(paramsDouble("stiffness")), damping(paramsDouble("damping")), friction(paramsDouble("friction")), cohesion(paramsDouble("cohesion")), density(paramsDouble("density")), color(" 1"){
		};
	~CProperty(){}
	double stiffness, damping, friction, cohesion;
	double density;
	string color;
 	private:
	};


class CParticle{
	CParticle(const CParticle&);
	public:
	//template<GType shapeType>
	//CParticle(const vec & _x0, double r):shape(new GeomObject<shapeType>(_x0, r)), q(1.0, 0.0, 0.0, 0.0), id(-1), forces(vec(0.0)), frozen(false){init();}
	template<GType shapeType>
	explicit CParticle(const GeomObject<shapeType> &_shape):shape(new GeomObject<shapeType>(_shape)),  id(-1),  forces(vec(0.0)), state(ready_to_go){init();}
	~CParticle(){
		delete shape;
		}

	void init(){
		x(0)=shape->Xc;
		w(0)=0.0;
		for(int i=1; i<6; ++i){
			x(i)=0.0;
			w(i)=0.0;
			}
		
		mass=material.density*shape->vol();
		Ixx=mass*shape->I(vec(1,0,0));
		Iyy=mass*shape->I(vec(0,1,0));
		Izz=mass*shape->I(vec(0,0,1));

		};

	double kEnergy(){
		return 0.5*mass*x(1).abs2();
		}

	double pEnergy(const vec &g){
		return -mass*(g*x(0));
		}
	double rEnergy(){
		vec wp=shape->q.toBody(w(1));
		return 0.5*(Ixx*wp(0)*wp(0)+Iyy*wp(1)*wp(1)+Izz*wp(2)*wp(2));
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
			shape->parse(in);
			x(0)=shape->Xc;
			//mass=material.density*4.0/3.0*M_PI*radius*radius*radius;
			}

	void calPos(double dt);
	void calVel(double dt);

	GeomObjectBase *shape;
	CProperty material;
	CDFreedom<3> x, x0, x_p;//TranslationalDFreedom;
	CDFreedom<3> w, w0, w_p;//Rotational;
	long id;
	vec test;
	//Quaternion q;//orientation
	vec forces, avgforces;
	vec torques, avgtorque;
	tState state;
	protected:
	double mass, Ixx, Iyy, Izz;
 	private:
	CDFreedom<5> RotationalDFreedom;
	};

ostream &operator <<(ostream &out, const CParticle &p){
	p.shape->print(out);
	//out<<"  "<<p.x(1);
	return out;
	}

double friction=1;

//using beeman method
void CParticle::calPos(double dt){
TRY
	static const double c=1./6.0;
	//translational degree
	x(0) += x(1)*dt + x(2)*(dt*dt*4.0*c) - x0(2)*(dt*dt*c);
	x(1) += x(2)*(dt*5.0*c) - x0(2)*(dt*c);

	//rotational degree
	static vec wp;
	static Quaternion dq(0,0,0,0);

	w(1) += w(2)*(dt*5.0*c) - w0(2)*(dt*c);
	wp=shape->q.toBody(w(1));//FIXME make sure which should be used
	//wp=w(1);			or this
	dq.u =    -shape->q.v(0)*wp(0) - shape->q.v(1)*wp(1) - shape->q.v(2)*wp(2);
	dq.v(0) =  shape->q.u  * wp(0) - shape->q.v(2)*wp(1) + shape->q.v(1)*wp(2);
	dq.v(1) =  shape->q.v(2)*wp(0) + shape->q.u *  wp(1) - shape->q.v(0)*wp(2);
	dq.v(2) = -shape->q.v(1)*wp(0) + shape->q.v(0)*wp(1) + shape->q.u *  wp(2);

	shape->q+=dq*dt*0.5;
	shape->q.normalize();
	
	
	//cerr<< shape->q <<endl;
	shape->rotateTo(shape->q);
	shape->moveto(x(0));
CATCH
	}

void CParticle::calVel(double dt){
	static const double c=1./6.0;
	x0(2)=x(2);
	x(2)=forces/mass;
	x(1)+= x(2)*(dt*2*c);
	
	w0(2)=w(2);
	static vec wp, wwp, torquep;
	torquep=shape->q.toBody(torques);

	wp=shape->q.toBody(w(1));
	wwp(0)=(torquep(0)+wp(1)*wp(2)*(Iyy-Izz))/Ixx;
	wwp(1)=(torquep(1)+wp(0)*wp(2)*(Izz-Ixx))/Iyy;
	wwp(2)=(torquep(2)+wp(0)*wp(1)*(Ixx-Iyy))/Izz;
	w(2)=shape->q.toWorld(wwp);


	w(1)+=w(2)*(dt*2*c);
	}


#endif /* PARTICLE_H */
