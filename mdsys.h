#ifndef MDSYS_H
#define MDSYS_H 
#include "common.h"
#include"grid.h"
#include"particle.h"
#include"interaction.h"

class ParticleContactHolder : public ShapeContact{
	public:
	ParticleContactHolder(CParticle &_p1, CParticle & _p2):ShapeContact(), p1(&_p1), p2(&_p2){}
	ParticleContactHolder(CParticle *_p1, CParticle * _p2):ShapeContact(), p1(_p1), p2(_p2){}
	ParticleContactHolder():ShapeContact(), p1(NULL), p2(NULL){}

	CParticle * const p1, * const p2;
 	private:
	};

typedef vector<CParticle *> ParticleContainer;

typedef GeomObjectBase * BasePtr;
class CSys{
	CSys();
	public:
	CSys(unsigned long maxnparticle):t(0), outDt(0.01), box(vec(0.0), vec(1.0)), maxr(0),  G(vec(0.0)), maxNParticle(maxnparticle),epsFreeze(1.0e-12), outEnergy("log_energy"){
	TRY
		pairs=new ParticleContactHolder*[maxNParticle*maxNParticle];
		for(unsigned int i=0; i<maxNParticle*maxNParticle; i++) pairs[i]=NULL;
	CATCH
		};
	~CSys();

	void initialize();
	void solve(double tMax, double dt);
	void forward(double dt);
	void calForces();
	void interactions();
	inline bool interact(unsigned int i, unsigned int j)const; //force from p2 on p1
	inline bool interact(CParticle *p1, GeomObject<tbox> *p2)const;
	vec contactForce(const Contact &c, const vec &dv, CProperty &m)const;

	int read_packing(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing2(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing3(string infilename, const vec &shift=vec(), double scale=1);
	void write_packing(string infilename);
	//void setup_grid(double d);

	bool add(CParticle *p);
	inline bool exist(int i);

	double t, outDt;
	ParticleContainer particles;

	
	ParticleContactHolder *  &pair(size_t i, size_t j)const{
		return pairs[j+i*particles.size()];
			}
	ParticleContactHolder **pairs;
	GeomObject<tbox> box;
	CPlane *sp;
	//CRecGrid *grid;
	double maxr;
	vec G;
	const unsigned maxNParticle;
 	private:
	double epsFreeze;
	ofstream outEnergy;
	};

CSys::~CSys(){
TRY
	for(size_t i=0; i<particles.size(); i++)
	for(size_t j=0; j<particles.size(); j++){
		delete pair(i,j);
		}
	delete [] pairs;
	ParticleContainer::iterator it;
	for(it=particles.begin(); it!=particles.end(); ++it){
		delete (*it);
		}
CATCH
	}

void CSys::initialize(){
TRY
	//create for each pair an object for the contact
	for(unsigned int i=0; i<particles.size(); i++){
		for(unsigned int j=0; j<particles.size(); j++){
			pair(i,j)=new ParticleContactHolder(particles.at(i), particles.at(j));
		//	pair(i,j).set(particles.at(i), particles.at(j));
			}
		}
CATCH
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

void CSys::calForces(){
TRY
	ParticleContainer::iterator it1, it2;
	//reset forces
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		(*it1)->forces=G*((*it1)->get_mass());
		(*it1)->torques=0.0;
		}

	//interactions
	for(unsigned int i=0; i<particles.size(); i++){
		for(unsigned int j=i+1; j<particles.size(); j++){
			if(interact(i, j)){ }
			}
		//the walls
		if(interact(particles.at(i), &box)){  }
		}
CATCH
};

inline
vec CSys::contactForce(const Contact &c, const vec &dv, CProperty &m)const{
TRY
	double proj=(dv*c.n);
	double ksi=c.dx_n;
	

	ksi=(m.stiffness*ksi+m.damping*proj)*sqrt(ksi); 
	if(ksi<0)ksi=0;//to eliminate artifical attractions
	vec fn=-ksi*c.n;//normal force
	vec ft=m.friction*(fn.abs())*((dv - proj*c.n));//dynamic frictions
	//cerr<< ft.abs()/(fn+ft).abs()<<"   "<<ft.abs()<<"  "<<fn.abs()<<endl;

	//static int i2=0;
	//i2++;
	//if(i2==160)exit(0);
	//cerr<< t<<"  "<<setprecision(14)<<c.x<<"  "<<c.n<<"  "<< c.dx_n << "  dv="<<dv<<"  Fn="<<fn<<endl;

	return fn+ft;//visco-elastic Hertz law
CATCH
	}

inline bool CSys::interact(unsigned int i,unsigned int  j)const{
TRY
	CParticle *p1=particles.at(i);
	CParticle *p2=particles.at(j);
	ShapeContact &overlaps=*pair(i,j);
	overlaps.clear();
	CInteraction::overlaps(&overlaps, p1->shape, p2->shape);
	static vec r1, r2, v1, v2, dv, force, torque(0.0);
	//static double proj, ksi;
	if(overlaps.size()==0)return false;
	for(size_t ii=0; ii<overlaps.size(); ii++){
		r1=overlaps(ii).x-p1->x(0);
		r2=overlaps(ii).x-p2->x(0);
		v1=p1->x(1)-cross(p1->w(1), r1);
		v2=p2->x(1)-cross(p2->w(1), r2);
		dv=v1-v2;

		force=contactForce(overlaps(ii), dv, p1->material);
		p1->shape->fixToBody(HomVec(overlaps(ii).x,1));


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
	static double Energy=0.0, rEnergy=0, pEnergy=0, kEnergy=0;

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
			Energy=rEnergy+kEnergy+pEnergy;
			outEnergy<<setprecision(14)<<t<<"  "<<Energy<<"  "<<kEnergy<<"  "<<pEnergy<<"  "<<rEnergy <<endl;
			rEnergy=0; pEnergy=0; kEnergy=0; Energy=0;
			}
	//bool allforwarded=false;
	for(it=particles.begin(); it!=particles.end(); ++it){
	//	if(!(*it)->frozen) 
		(*it)->calPos(dt);
		//if(!it->frozen) it->x.gear_predict<4>(dt);
		}
	calForces();
//	if(!allforwarded)foward(dt/2.0, 2);

	Energy=0.0, rEnergy=0, pEnergy=0, kEnergy=0;
	for(it=particles.begin(); it!=particles.end(); ++it){
	//	if(!(*it)->frozen) 
		(*it)->calVel(dt);
		if(0)if((*it)->x(1).abs()< epsFreeze  && (*it)->avgforces.abs()< epsFreeze ) {
			//(*it)->frozen=true;
			(*it)->material.color="0.5 0.5 0.5";
			(*it)->x(1)=0.0;
			}
		rEnergy+=(*it)->rEnergy();
		pEnergy+=(*it)->pEnergy(G);
		kEnergy+=(*it)->kEnergy();
		//cout<< it->x <<"  "<<it->size<< " cir"<<endl;
		}
	count++;
	if(out.is_open())out.close();
CATCH
	}

void CSys::solve(double tMax, double dt){
	try{
	initialize();
	bool stop=false;
	while(true){
		if(t+dt>tMax){//to stop exactly at tMax
			dt=tMax-t;
			stop=true;
			}

		//calForces();
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
	static ShapeContact overlaps;
	overlaps.clear();
	CInteraction::overlaps(&overlaps, p1->shape, (GeomObjectBase*)p2);

	static vec dv, r1, force, torque, vt, vn;
	if(overlaps.size()==0)return false;
	for(size_t i=0; i<overlaps.size(); i++){
		r1=overlaps(i).x-p1->x(0);
		dv=p1->x(1)+cross(p1->w(1), r1);

		force=contactForce(overlaps(i), -dv, p1->material);
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
