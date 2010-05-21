#ifndef MDSYS_H
#define MDSYS_H 
#include "common.h"
//#include"grid.h"
#include"particle.h"
#include"interaction.h"

class ParticleContactHolder : public ShapeContactHolder{
	public:
	ParticleContactHolder(const CParticle &_p1, const CParticle & _p2, bool b=false):ShapeContactHolder(), p1(&_p1), p2(&_p2),persist(b){}
	ParticleContactHolder(const CParticle *_p1=NULL, const CParticle * _p2=NULL, bool b=false):ShapeContactHolder(), p1(_p1), p2(_p2),persist(b){}

	const CParticle *p1, *p2;
	
	bool persist;
 	private:
	};

typedef vector<CParticle *> ParticleContainer;

typedef GeomObjectBase * BasePtr;
class CSys{
	public:
	CSys(unsigned long maxnparticle):t(0), outDt(0.01), box(vec(0.0), vec(1.0)), maxr(0),  G(vec(0.0)), maxNParticle(maxnparticle),epsFreeze(1.0e-12){
		pairs=new ParticleContactHolder[maxNParticle*maxNParticle];
		};
	~CSys();

	void initialize();
	void solve(double tMax, double dt);
	void forward(double dt);
	void calForces();
	void interactions();
	void cleaninteractions();
	inline bool interact(unsigned int i, unsigned int j)const; //force from p2 on p1
	inline bool interact(CParticle *p1, GeomObject<tbox> *p2)const;

	int read_packing(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing2(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing3(string infilename, const vec &shift=vec(), double scale=1);
	void write_packing(string infilename);
	ParticleContactHolder & pair(int i, int j){
		return pairs[j+i*particles.size()];
		}
	//void setup_grid(double d);

	bool add(CParticle *p);
	inline bool exist(int i);

	double ke;//kinetic energy
	double t, outDt;
	ParticleContainer particles;
	GeomObject<tbox> box;
	CPlane *sp;
	//CRecGrid *grid;
	double maxr;
	vec G;
	const unsigned maxNParticle;
 	private:
	ParticleContactHolder *pairs;
	double epsFreeze;
	};

CSys::~CSys(){
	ParticleContainer::iterator it;
	for(it=particles.begin(); it!=particles.end(); ++it){
		delete (*it);
		}
	delete [] pairs;
	}

void CSys::initialize(){
	//create for each pair an object for the contact
	for(unsigned int i=0; i<maxNParticle; i++)
		for(unsigned int j=0; j<maxNParticle; j++){
			pair(i,j).clear();
			}
	}

//void CSys::setup_grid(double _d){
		//grid=new CRecGrid(box.corner, box.L, _d*maxr);
		//}

bool CSys::add(CParticle *p){
TRY
	ERROR(particles.size()==maxNParticle,"Reached max number of particles.");
	vec force(0.0);
	//if(interact(p, box, force)){
		//ERROR("The particle is intially within the box: "<<p->x(0));
		//return false;
		//}
	if(maxr<p->shape->radius)maxr=p->shape->radius;
	particles.push_back(p);
	//assert(grid);
	//grid->add(p);
	return true;
CATCH
	}

void CSys::cleaninteractions(){
TRY
	//create for each pair an object for the contact
	for(unsigned int i=0; i<maxNParticle; i++)
		for(unsigned int j=0; j<maxNParticle; j++){
			//if(! pair(i,j)->persist) pair(i,j).clear();
			pair(i,j).clear();
			}
	return;
CATCH
	}

void CSys::calForces(){
TRY
	ParticleContainer::iterator it1, it2;
	//reset forces
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		(*it1)->forces=G*((*it1)->get_mass());
		(*it1)->torques=0.0;
		}

	//interactions
	cleaninteractions();
	for(unsigned int i=0; i<particles.size(); i++){
		for(unsigned int j=i+1; j<particles.size(); j++){
		//ittemp=it1;++ittemp;
		//for(it2=ittemp; it2!=particles.end(); ++it2){
			//if((*it1)->frozen && (*it2)->frozen)continue;
			if(interact(i, j)){ }
				//rescale_ellipse_to_touch(*(static_cast<CEllipsoid*> ((*it1)->shape)), *(static_cast<CEllipsoid*> ((*it2)->shape)));
				}
		//the walls
		if(interact(particles.at(i), &box)){  }
		}
CATCH
};

inline
vec contactForce(const vec dx, const vec dv, double stiff, double damp){
	double proj=(dv*dx.normalized());
	double ksi=dx.abs();
	ksi=(stiff*ksi+damp*proj);//*sqrt(ksi); //to eliminate artifical attractions
	if(ksi<0)ksi=0;
	return -ksi*dx;//visco-elastic Hertz law
	}

inline bool CSys::interact(unsigned int i,unsigned int  j)const{
TRY
	CParticle *p1=particles.at(i);
	CParticle *p2=particles.at(j);
	//if(pairs[i+particles.size()*j]==NULL) pairs[i+particles.size()*j]=new ParticleContactHolder*[maxNParticle*maxNParticle];
	ShapeContactHolder &overlaps=pairs[i+particles.size()*j];
	CInteraction::overlaps(&overlaps, p1->shape, p2->shape);
	static vec r1, r2, v1, v2, dv, force, torque(0.0);
	//static double proj, ksi;
	if(overlaps.size()==0)return false;
	for(size_t i=0; i<overlaps.size(); i++){
		r1=overlaps.at(i).x-p1->x(0);
		r2=overlaps.at(i).x-p2->x(0);
		v1=p1->x(1)+cross(r1, p1->w(1));
		v2=p2->x(1)+cross(r2, p2->w(1));
		dv=v1-v2;

		force=contactForce(overlaps.at(i).dx, dv, p1->material.stiffness, p1->material.damping);

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


void CSys::forward(double dt){
TRY
	static int count=0, outN=0,outPutN=outDt/dt;
	static ofstream out;

	ParticleContainer::iterator it;
	if(count%outPutN==0){
			stringstream outname;
			outname<<"out"<<setw(5)<<setfill('0')<<outN;
			out.open(outname.str().c_str());
			box.print(out);
			gout=&out;
			for(it=particles.begin(); it!=particles.end(); ++it){
				out<<**it<<endl;
				}
			count=0;
			outN++;
			}
	//double energy=0.0;
	//bool allforwarded=false;
	for(it=particles.begin(); it!=particles.end(); ++it){
	//	if(!(*it)->frozen) 
		(*it)->calPos(dt);
		//if(!it->frozen) it->x.gear_predict<4>(dt);
		//energy+=(*it)->rEnergy()+(*it)->kEnergy()+(*it)->pEnergy(G);
		}
	calForces();
//	if(!allforwarded)foward(dt/2.0, 2);

	ke=0;
	for(it=particles.begin(); it!=particles.end(); ++it){
	//	if(!(*it)->frozen) 
		(*it)->calVel(dt);
		//vec dA=(it->forces/it->get_mass()-it->x(2));
		//if(!it->frozen) it->x.gear_correct<5>(dt,dA);
		if(0)if((*it)->x(1).abs()< epsFreeze  && (*it)->avgforces.abs()< epsFreeze ) {
			//(*it)->frozen=true;
			(*it)->material.color="0.5 0.5 0.5";
			(*it)->x(1)=0.0;
			}
		ke+=(*it)->kEnergy();
		//cout<< it->x <<"  "<<it->size<< " cir"<<endl;
		}
//output 
	count++;
	if(out.is_open())out.close();
CATCH
	}

void CSys::solve(double tMax, double dt){
	try{
	bool stop=false;
	initialize();
	while(true){
		if(t+dt>tMax){//to stop exactly at tMax
			dt=tMax-t;
			stop=true;
			}

		calForces();
		forward(dt);
		t+=dt;
		if(stop)break;
		}
	}catch(CException e){
		ERROR(1,"Some error in the solver at t= "+ stringify(t)+"\n\tfrom "+e.where());
		}
	}

void CSys::write_packing(string outfilename){
	ParticleContainer::iterator it;
	ofstream out(outfilename.c_str());
	assert(!out);
	for(it=particles.begin(); it!=particles.end(); ++it){
		out<<**it<<endl;
		}
	}
/*
int CSys::read_packing3(string infilename, const vec &shift, double scale){
        cerr<< "Reading contacts ..." <<endl;
        ifstream inputFile(infilename.c_str());
        if(!inputFile.good())
        {
        std::cerr << "Unable to open input file: " << infilename << std::endl;
        return 0;
        }

        string line;
        string vname;

	double a1[]={1020.0, 1020.0};
	vec center(a1);
	vec dist(0.0);
        //Parse the line
        while(getline(inputFile,line))
                {
                stringstream ss(line);
                CParticle *p=new CParticle(vec(0.0), 0);
		p->parse(ss);
		p->frozen=false;//FIXME if it is not frozen it doesnt work, why?
		
		p->x(0)+=shift;
		p->x(0)*=scale;
		double ran=drand48();
                add(p);
                }

        cerr<< "done" <<endl;
        inputFile.close();
	return particles.size();
        }

int CSys::read_packing2(string infilename, const vec &shift, double scale){
        cerr<< "Reading contacts ..." <<endl;
        ifstream inputFile(infilename.c_str());
        if(!inputFile.good())
        {
        std::cerr << "Unable to open input file: " << infilename << std::endl;
        return 0;
        }

        string line;
        string vname;

	double a1[]={1020.0, 1020.0};
	vec center(a1);
	vec dist(0.0);
        //Parse the line
        while(getline(inputFile,line))
                {
                stringstream ss(line);
                CParticle *p=new CParticle(vec(0.0), 0);
                ss>>p->x(0)(0)>>p->x(0)(1)>>p->x(0)(2)>>p->shape->radius;
		//if(exist(p->id))continue;
		p->frozen=true;//FIXME if it is not frozen it doesnt work, why?
		p->shape->identifier=0;
		if(particles.size()>1000)break;
		
		p->x(0)+=shift;
		p->x(0)*=scale;
		double ran=drand48();
		if(ran<0.05)p->shape->identifier=1;
		else if(ran<0.2)p->shape->identifier=2;
		else p->shape->identifier=3;
		p->material.color=stringify(drand48())+stringify(drand48())+stringify(drand48());
		p->shape->Xc=p->x(0);
		p->x(1)=0.0;
		p->x(2)=0.0;
		p->x0(1)=0.0;
		p->x0(2)=0.0;
                add(p);
                }

        cerr<< "done" <<endl;
        inputFile.close();
	return particles.size();
        }

//int CSys::read_packing(string infilename, const vec &shift, double scale){
int CSys::read_packing(string infilename, const vec &shift, double scale){
        cerr<< "Reading contacts ..." <<endl;
        ifstream inputFile(infilename.c_str());
        if(!inputFile.good())
        {
        std::cerr << "Unable to open input file: " << infilename << std::endl;
        return 0;
        }

        string line;
        string vname;

// 413 1026 1638  
//406 1021 1635 
// 191 1845
	double a1[]={1020.0, 1020.0};
	vec center(a1);
	vec dist(0.0);
        //Parse the line
        while(getline(inputFile,line))
                {
                stringstream ss(line);
                CParticle *p=new CParticle(vec(0.0), 0);
                ss>>p->id>>p->x(0)(0)>>p->x(0)(1)>>p->x(0)(2)>>p->shape->radius;
		//if(exist(p->id))continue;
		if(p->x(0)(2)<270 || p->x(0)(2)>1725)p->frozen=true;
		dist(0)=p->x(0)(0); dist(1)=p->x(0)(1);
		if((dist-center).abs()> 525)p->frozen=true;
		
		p->x(0)+=shift;
		p->x(0)*=scale;
		//p.shape->radius=44;
		p->shape->radius*=scale*1.04;
		p->x0(0)=p->x(0);
                add(p);
                }
        cerr<< "done" <<endl;
        inputFile.close();
	return particles.size();
        }

*/
bool CSys::exist(int i){
	ParticleContainer::iterator it1;
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		if((*it1)->id==i)return true;
		}
		return false;
		}

void CSys::interactions(){
	ParticleContainer::iterator it1, it2, ittemp;
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
	ittemp=it1;++ittemp;
	for(it2=ittemp; it2!=particles.end(); ++it2){
		double d=((*it1)->x(0)-(*it2)->x(0)).abs()- (*it1)->shape->radius - (*it2)->shape->radius;
		if(d<0)cerr<< d<<" "<< (*it1)->id<<"  "<< (*it2)->id <<endl;
		}
		}

		}


inline bool CSys::interact(CParticle *p1, GeomObject<tbox> *p2)const{
TRY
	ShapeContactHolder overlaps;
	CInteraction::overlaps(&overlaps, p1->shape, (GeomObjectBase*)p2);

	static vec dv, r1, force, torque, vt, vn;
	if(overlaps.size()==0)return false;
	for(size_t i=0; i<overlaps.size(); i++){
		r1=overlaps.at(i).x-p1->x(0);
		dv=p1->x(1)+cross(r1, p1->w(1));

		force=contactForce(overlaps.at(i).dx, dv, p1->material.stiffness, p1->material.damping);
		p1->addforce(force);

		torque=cross(r1, force);
		p1->addtorque(torque);

		//tangent=(1.0-fabs(proj)/dv.abs()/overlaps.at(i)->dx.abs())*dv;

		//cerr<< tangent <<endl;
		//force=200*tangent;
		//p1->addforce(force);
		//torque=cross(r1, force);
		//p1->addtorque(torque);
		}

	return true;
CATCH
	}

#endif /* MDSYS_H */
