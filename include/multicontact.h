// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef SHAPECONTACT_H
#define SHAPECONTACT_H 
#include<vector>
#include"exception.h"
#include"contact.h"


//This class is used to be able to work with shapes for which more than one 
//contact per pair is possible. For spheres and ellipsoids this is not the case.

template <class T>
class MultiContact : public std::vector<Contact>
	{
	public:
	explicit MultiContact(T *_p1=NULL, T *_p2=NULL, const CPlane &p=CPlane(vec(0,0,0), vec(0,0,1)) )
		:std::vector<Contact >(), p1(_p1), p2(_p2), plane(p), x1(HomVec(0,0,0,1)), x2(HomVec(0,0,0,1)), has_sep_plane(false), set(false)
		{
		N++;
		}

	void add(const Contact &contact)
		{
		//if(this->size()>0)this->at(0)=contact;
		//else 
		push_back(contact);
		}

	Contact &operator()(size_t i)
		{
		return (this->at(i));
		}

	~MultiContact()
		{
		N--;
		}

	T *const p1, * const p2;
	CPlane plane;
	HomVec x1, x2;
	HomVec x01, x02;
	bool has_sep_plane;
	bool set;
 	private:
	static long N;
	};
	template<class T>
	long MultiContact<T>::N=0;

typedef MultiContact<GeomObjectBase> ShapeContact ;

#endif /* SHAPECONTACT_H */
