#ifndef VERLET_H
#define VERLET_H 
#include<map>
#include"exception.h"
#include"vec.h"
#include"shapecontact.h"
#include"particlecontact.h"

using namespace std;
template<class T>
//class CVerlet: public list<T*> 
class CVerlet: public map<T*, ParticleContactHolder<T> > 
	{
	public:

	CVerlet(T *p):self_p(p), set(false)
		{
		ERROR(p==NULL,"Improper initialization of verlet list");
		}

	void add(T *p)
		{
		ERROR(p==self_p, "A particle cannot be added to its own verlet list");
		//push_back(p);
		insert(pair<T*, ParticleContactHolder<T> > (p, ParticleContactHolder<T>(this->self_p, p)));
		}

	vec x; //position when the list was updated

	T * self_p;
	bool set;
 	private:
	};

#endif /* VERLET_H */
