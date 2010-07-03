#ifndef PARTICLECONTACT_H
#define PARTICLECONTACT_H 

template<class particleT>
class ParticleContactHolder : public ShapeContact{
	public:
	ParticleContactHolder(particleT &_p1, particleT& _p2):ShapeContact(), p1(&_p1), p2(&_p2){}
	ParticleContactHolder(particleT *_p1, particleT * _p2):ShapeContact(), p1(_p1), p2(_p2){}
	ParticleContactHolder():ShapeContact(), p1(NULL), p2(NULL){}

	particleT * const p1, * const p2;
 	private:
	};

#endif /* PARTICLECONTACT_H */
