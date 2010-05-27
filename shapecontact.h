#ifndef SHAPECONTACT_H
#define SHAPECONTACT_H 
#include<vector>
#include"contact.h"
#include"exception.h"

class ShapeContact : public std::vector<Contact>
	{
	public:
	ShapeContact(GeomObjectBase *_p1=NULL, GeomObjectBase *_p2=NULL, const CPlane &p=CPlane(vec(0,0,0), vec(0,0,1)) )
		:std::vector<Contact >(), p1(_p1), p2(_p2), plane(p), x1(HomVec(0,0,0,1)), x2(HomVec(0,0,0,1)), set(false)
		{
		N++;
		}

	void add(const Contact &contact){
		if(this->size()>0)this->at(0)=contact;
		else push_back(contact);
		}
	Contact &operator()(size_t i){
		return (this->at(i));
		}

	~ShapeContact(){
		N--;
		}
	GeomObjectBase  *const p1, * const p2;
	CPlane plane;
	HomVec x1, x2;
	HomVec x01, x02;
	bool set;
 	private:
	static long N;
	};
	long ShapeContact::N=0;


#endif /* SHAPECONTACT_H */
