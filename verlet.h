#ifndef VERLET_H
#define VERLET_H 
#include<set>
#include<map>
#include"exception.h"
#include"vec.h"
#include"shapecontact.h"

using namespace std;
template<class T>
class CVerlet: public list<T*> 
	{
	public:

	CVerlet(T *p, double r1=0, double r2=0):self_p(p)
		{
		ERROR(p==NULL,"Improper initialization of verlet list");
		}

	void add(T *p)
		{
		ERROR(p==self_p, "A particle cannot be added to its own verlet list");
		push_back(p);
		}

	vec x; //position when the list was updated
	
	

 	private:
	T * self_p;
	};

#endif /* VERLET_H */