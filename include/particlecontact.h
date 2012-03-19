// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef PARTICLECONTACT_H
#define PARTICLECONTACT_H 

template<class particleT>
class ParticleContactHolder : public ShapeContact{
	public:
	ParticleContactHolder(particleT &_p1, particleT& _p2):ShapeContact(), p1(&_p1), p2(&_p2), in_contact(false){}
	ParticleContactHolder(particleT *_p1, particleT * _p2):ShapeContact(), p1(_p1), p2(_p2), in_contact(false){}
	ParticleContactHolder():ShapeContact(), p1(NULL), p2(NULL), in_contact(false){}

	particleT * const p1, * const p2;
	bool in_contact;
	double col_time;
 	private:
	};

#endif /* PARTICLECONTACT_H */
