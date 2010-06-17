#ifndef VERLET_H
#define VERLET_H 
#include<set>
#include<list>
#include"exception.h"
#include"vec.h"

using namespace std;
template<class T>
class CVerlet: public list<T*> {
	public:
	CVerlet(T *p, double r1=0, double r2=0):self_p(p){
		ERROR(p==NULL,"Improper initialization of verlet list");
		}
	void add(T *p){
		ERROR(p==self_p, "A particle cannot be added to its own verlet list");
		//typename list<T*>::iterator it=this->find(p);
		//if(!it){
			//WARNING("The particle is already in the verlet list");
			//return;
			//}
		push_back(p);
		}
	//position when the list was updated
	vec x;
 	private:
	T * self_p;
	};

#endif /* VERLET_H */
