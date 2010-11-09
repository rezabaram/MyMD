#ifndef PARTICLE_H
#define PARTICLE_H 

#include<sstream>
#include<fstream>
#include<vector>
#include<set>
#include<string>
#include<algorithm>
#include<iomanip>

#include"common.h"
#include"params.h"
#include"define_params.h"
#include"dfreedom.h"
#include"shapes.h"
#include"verlet.h"
#include"grid.h"


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


class CParticle 
	{
	CParticle(); 
	protected:
	//CParticle(const CParticle&); // for making shallow copies for shadow particles
	public:
	GeomObjectBase *shape;
	template<class T>
	explicit CParticle(const T &_shape, bool sh=false)
	:shape(_shape.clone()),  id(-1),  state(ready_to_go),vlist(this),vlistold(this),isshadow(sh), orig_p(this)
	{
	TRY
		init();
	CATCH
	}

	virtual ~CParticle()
	{
	TRY
		shape->destroy();
		shape=NULL;
		clear();
	CATCH
		}

	virtual void clear(){
	TRY
		if(!isshadow){
			delete forces;
			delete torques;
			}
	CATCH
		}
	virtual long orig_id()const{
		return id;
		}

	virtual bool expire(){
		return false;
		}

	list<void *> shadows;
	CParticle *Shadow(void *plane);

	virtual void reset_forces(const vec &v=vec(0.0)){
		*forces=v;
		}
	virtual void reset_torques(const vec &v=vec(0.0)){
		*torques=v;
		}

	virtual void init(){
		ERROR(isshadow,"This is not allowed");

		forces= (new vec(0.0));
		torques=(new vec(0.0)); 

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

	virtual void  print(ostream &out){
		out<<id<<"  ";
		this->shape->print(out);
		out<<endl;
		}

	double get_mass()const{return mass;}

	void addforce(const vec force){
		static vec prev(0.0);
		*forces+=force;
		avgforces=(force);
		//avgforces=(prev+force)/2.0;
		prev=force;
		};

	void addtorque(const vec torque){
		*torques+=torque;
		};

	virtual void parse(std::istream &in){
			in>>id;
			shape->parse(in);
			x(0)=shape->Xc;
			//mass=material.density*4.0/3.0*M_PI*radius*radius*radius;
			}

	virtual void calPos(double dt);
	virtual void calVel(double dt);
	void get_grid_neighbours(set<CParticle *> &neigh)const;

	CProperty material;
	CDFreedom<3> x, x0, x_p;//TranslationalDFreedom;
	CDFreedom<3> w, w0, w_p;//Rotational;
	long id;
	vec test;
	//Quaternion q;//orientation
	vec *forces, avgforces;
	vec *torques, avgtorque;
	tState state;
	CVerlet<CParticle> vlist;
	CVerlet<CParticle> vlistold;
	vector< CNode3D<CParticle> *> grid_nodes;

	//to hold neighbours on the grid
	//set is chosen to avoid repeatition
	set<CParticle *> neighbours;
	double mass, Ixx, Iyy, Izz;
	bool isshadow;

	CParticle *orig_p;//original particle (if this is a shawdow)
	protected:
 	private:
	CDFreedom<5> RotationalDFreedom;
	};

void CParticle::get_grid_neighbours(set<CParticle *> &neigh)const{
	CNode3D<CParticle>::iterator it;
	for(size_t i=0; i<grid_nodes.size(); i++){
	for(it=grid_nodes.at(i)->begin(); it!= grid_nodes.at(i)->end(); it++){
		if(this!=(*it))neigh.insert(*it);
		}
		}
	
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
	
	
	shape->rotateTo(shape->q);
	shape->moveto(x(0));
CATCH
	}

void CParticle::calVel(double dt){
	static const double c=1./6.0;
	x0(2)=x(2);
	x(2)=*forces/mass;
	x(1)+= x(2)*(dt*2*c);
	
	w0(2)=w(2);
	static vec wp, wwp, torquep;
	torquep=shape->q.toBody(*torques);

	wp=shape->q.toBody(w(1));
	wwp(0)=(torquep(0)+wp(1)*wp(2)*(Iyy-Izz))/Ixx;
	wwp(1)=(torquep(1)+wp(0)*wp(2)*(Izz-Ixx))/Iyy;
	wwp(2)=(torquep(2)+wp(0)*wp(1)*(Ixx-Iyy))/Izz;
	w(2)=shape->q.toWorld(wwp);


	w(1)+=w(2)*(dt*2*c);
	}

class ShadowParticle : public CParticle
	{
	public:
	explicit ShadowParticle(CParticle *_p, CPlane *_plane) 
	:CParticle(*_p), plane(_plane)
	{
		
		ERROR(shape->type!=tellipsoid, "Incomplete construction of shadow particle.");
		orig_p=_p;
		init();
		}
	virtual void init(){
	TRY 	
		assert(orig_p);
		shape=orig_p->shape->clone();
		isshadow=true;
		shift=plane->vec_to_shadow;
		assert(shift.abs()>1e-10);
		shape->moveto(orig_p->shape->Xc+shift);
		forces=orig_p->forces;
		torques=orig_p->torques;
		vlist.clear();
		vlistold.clear();
		vlist.self_p=this;
		vlistold.self_p=this;
	CATCH
	}

	virtual ~ShadowParticle(){
		WARNING("out shadow");
		clear();
		}
	virtual void clear(){
		}
	virtual bool expire()
	{
	TRY
		return !(orig_p->shape->doesHit(*plane));
	CATCH
	}


	virtual long orig_id(){
		return orig_p->id;
		}

	virtual void  print(ostream &out){
		out<<orig_p->id<<"  ";
		this->shape->print(out);
		out<<endl;
		}

	virtual void reset_forces(const vec &v=vec(0.0)){
		//do nothing
		}
	virtual void reset_torques(const vec &v=vec(0.0)){
		//do nothing
		}
	virtual void parse(std::istream &in){
		WARNING("A shadow particle may not be parsed in directly");
			}

	virtual void calPos(double dt)
	{
	TRY
		x=orig_p->x;
		x(0)=orig_p->x(0)+shift;
		shape->q=orig_p->shape->q;
		shape->rotateTo(shape->q);
		shape->moveto(orig_p->shape->Xc+shift);
	CATCH
	};

	virtual void calVel(double dt)
	{
	TRY
		x(1)=orig_p->x(1);
		x(2)=orig_p->x(2);
		w=orig_p->w;
	CATCH
	};

	//mechanism for periodic boundary
	CPlane  *plane;//the plane being crossed
	vec shift;//the shift vector
 	private:
	};

	CParticle* CParticle::Shadow(void *plane)
	{
	TRY
		ERROR(isshadow, "the particle is itself a shadow");
		ERROR(!plane, "invalid pointer to plane");

		if( !((CPlane*)plane)->has_shadow )return NULL;
		if(std::find( shadows.begin(), shadows.end(), plane ) != shadows.end() )return NULL;//already has that shadow

		shadows.push_back(plane);
		return (new ShadowParticle(this, (CPlane*)plane));
		
	CATCH
	}

#endif /* PARTICLE_H */
