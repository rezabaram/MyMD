#ifndef INTERACTION_FORCE_H
#define INTERACTION_FORCE_H 
#include"particle.h"

class Test{
	public:
	static vec contactForce(const Contact &c, const vec &dv, CMaterial &m, double tmp=1){
	TRY
		double proj=(dv*c.n);
		double ksi=c.dx_n;

		ksi=(m.stiffness*ksi+tmp*m.damping*proj)*sqrt(ksi); 
		
		if(ksi<0)ksi=0; //to eliminate artifical attractions
		vec fn=-ksi*c.n; //normal force
		vec ft=-tmp*m.friction*(fn.abs())*((dv - proj*c.n));//dynamic frictions

		return fn+ft; //visco-elastic Hertz law
	CATCH
		}
	static bool interact(CParticle *p1,CParticle *p2){
TRY
	//CParticle *p1=particles.at(i);
	//CParticle *p2=particles.at(j);
	ShapeContact &overlaps=p1->vlist[p2];
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
