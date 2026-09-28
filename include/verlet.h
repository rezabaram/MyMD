// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef VERLET_H
#define VERLET_H 
#include<map>
#include"exception.h"
#include"vec.h"
#include"multicontact.h"
#include"particlecontact.h"
#include"packing.h"

using namespace std;

// Despite the name, and the file it lives in, this is NOT a Verlet list.  It is
// the per-particle contact cache: CParticle::vlist maps a neighbour particle to
// the contact data computed for that pair this step, which lets Test::interact
// hand the overlap geometry back and forth.  The actual Verlet neighbour list
// this file was written for (CVerletManager) was never compiled and has been
// removed -- see ROADMAP.md.
//
// TODO(design): rename to ContactCache and move out of verlet.h.
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

#endif /* VERLET_H */
