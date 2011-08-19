#ifndef INTERACTION_FORCE_H
#define INTERACTION_FORCE_H 
#include"particle.h"

class Test{
	public:
	static vec static_friction(Contact &c, double fn, vec vt, CProperty &m){
	TRY
		const GeomObjectBase *p1=static_cast<const GeomObjectBase*>(c.p1);
		const GeomObjectBase *p2=static_cast<const GeomObjectBase*>(c.p2);
		if(!c.static_friction_on){
			c.x1=p1->q.toBody(c.x);
			c.x2=p2->q.toBody(c.x); 
			c.static_friction_on=true;
			return vec(0,0,0);
			}
		vec x1=p1->q.toWorld(c.x1);
		vec x2=p2->q.toWorld(c.x2);
		vt.normalize();
		double project=vt*(c.x1-c.x2);
		double dx2=project*project;
		double f=m.static_friction*dx2;
		
		cerr<< f <<endl;
		if(f>m.friction_threshold){
			c.x1=p1->q.toBody(c.x);
			c.x2=p2->q.toBody(c.x); 
			c.static_friction_on=false;
			return vec(0,0,0);
			}
		return f*vt;
	CATCH
		}
	static vec contactForce(Contact &c, const vec &dv, CProperty &m, double tmp=1){
	TRY
		double proj=(dv*c.n);
		double ksi=c.dx_n;

		ksi=(m.stiffness*ksi+tmp*m.damping*proj)*sqrt(ksi); 
		
		if(ksi<0)ksi=0; //to eliminate artifical attractions
		vec fn=-ksi*c.n; //normal force
		vec vt=(dv - proj*c.n);
		vec ft=-tmp*m.friction*(fn.abs())*vt;//dynamic frictions
		vec fs=static_friction(c,fn.abs(), vt, m);

		return fn+ft+fs; //visco-elastic Hertz law
	CATCH
		}
	static bool interact(CParticle *p1,CParticle *p2){
TRY
	//CParticle *p1=particles.at(i);
	//CParticle *p2=particles.at(j);
	ShapeContact &overlaps=p1->vlist[p2];
	//ERROR(p1->vlist.size()>25, "You need to clean this list");
	//ShapeContact overlaps;

	overlaps.clear();
	CInteraction::overlaps(&overlaps, p1->shape, p2->shape);
	static vec r1, r2, v1, v2, dv, force, torque(0.0);
	//static double proj, ksi;
	if(overlaps.size()==0)return false;
	for(size_t ii=0; ii<overlaps.size(); ii++){
		r1=overlaps(ii).x-p1->x(0);
		r2=overlaps(ii).x-p2->x(0);
		v1=p1->x(1)+cross(p1->w(1), r1);
		v2=p2->x(1)+cross(p2->w(1), r2);
		dv=v1-v2;

		force=Test::contactForce(overlaps(ii), dv, p1->material);
		p1->shape->fixToBody(HomVec(overlaps(ii).x,1));
		//cerr<< (p2->x(0)-p1->x(0)).normalized() <<endl;


		p1->addforce(force);
		torque=cross(r1, force);
		p1->addtorque(torque);
			
		force*=-1.0;
		p2->addforce(force);
		torque=cross(r2, force);
		p2->addtorque(torque);

		}

	return true;
CATCH
	}
};
#endif /* INTERACTION_FORCE_H */
