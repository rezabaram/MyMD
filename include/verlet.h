#ifndef VERLET_H
#define VERLET_H 
#include<map>
#include"exception.h"
#include"vec.h"
#include"multicontact.h"
#include"particlecontact.h"
#include"packing.h"

using namespace std;
template<class particleT>
//class CVerletList: public list<particleT*> 
class CVerletList: public map<particleT*, ParticleContactHolder<particleT> > 
	{
	public:

	CVerletList(particleT *p):self_p(p), set(false)
		{
		ERROR(p==NULL,"Improper initialization of verlet list");
		}

	void add(particleT *p)
		{
		ERROR(p==self_p, "A particle cannot be added to its own verlet list");
		//push_back(p);
		insert(pair<particleT*, ParticleContactHolder<particleT> > (p, ParticleContactHolder<particleT>(this->self_p, p)));
		}

	vec x; //position when the list was updated

	particleT * self_p;
	bool set;
 	private:
	};

template<class particleT>
class CVerletManager{
	public:
	CVerletManager(CPacking<particleT> *p):need_update (true), packing(p){}

	bool add_particle(particleT *p){
		packing->push_back(p);
		setup(packing->back());
		return true;
		}
	void setup(particleT *p);
	void print(ostream &out);
	void update();

	void set_distance(double d){
		distance=d;
		min_distance=distance/5;
		max_distance=distance*5;
		}

	double distance, min_distance, max_distance;
	bool need_update;
	private:
	CPacking<particleT> *packing;
	};

template<class particleT>
void CVerletManager<particleT>::print(ostream &out){
	typename CPacking<particleT>::iterator it;
	typename CVerletList<particleT >::iterator neigh;
	for(it=packing->begin(); it!=packing->end(); ++it){
		out<< (*it)->id <<": ";
		for(neigh=(*it)->vlist.begin(); neigh!=(*it)->vlist.end(); ++neigh ){
			out<< (*neigh).first->id <<" ";
			}
			out<< endl;
		}
	}
//construct the verlet list of all particles
template<class particleT>
void CVerletManager<particleT>::update(){
	distance*=0.99;
	if(distance<min_distance)distance=min_distance;

	typename CPacking<particleT>::iterator it;
	for(it=packing->begin(); it!=packing->end(); ++it){
		if(need_update==false and ((*it)->x(0)-(*it)->vlist.x).abs() > distance/2.0-epsilon) need_update=true;
		}

	if(!need_update) return;

	for(it=packing->begin(); it!=packing->end(); it++){
		(*it)->vlistold=(*it)->vlist;
		(*it)->vlist.clear();
		setup(*it);
		}

	distance*=1.5;
	if(distance>max_distance)distance=max_distance;

	need_update=false;
	//cerr<< "Verlet updated at: "<<t <<"\t verlet distance: "<<verlet_distance<<endl;
	}

//construct the verlet list of one particle 
template<class particleT>
void CVerletManager<particleT>::setup(particleT *p){
TRY
	typename CPacking<particleT>::iterator it;
	typename CVerletList<particleT>::iterator v_it_old;
	for(it=packing->begin(); (*it)!=p; it++){ //checking particles before in the list
	//for(it=particles.begin(); it!=particles.end(); it++){ //checking particles before in the list
		if(p->min_distance(*it) < distance){
			v_it_old=p->vlistold.find(*it);
			if(v_it_old!= p->vlistold.end())
				p->vlist.insert(*v_it_old);
				else
				p->vlist.add((*it));
				}
		}
		assert((*it)->id == p->id);
	p->vlist.x=p->x(0);//save the position at which the list has been updated
	p->vlist.set=true;
CATCH
	}
#endif /* VERLET_H */
