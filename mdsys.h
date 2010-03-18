#ifndef MDSYS_H
#define MDSYS_H 
#include"common.h"
#include"grid.h"
#include"particle.h"

typedef GeomObjectBase * BasePtr;
class CSys{
	public:
	CSys():t(0), outDt(0.01),epsFreeze(1.0e-12), box(vec(0.0), vec(1.0)), maxr(0),  G(vec(0.0)){};
	~CSys();

	bool add(CParticle *p);
	void solve(double tMax, double dt);
	void forward(double dt);
	void calForces();
	vector<CParticle *> particles;
	int read_packing(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing2(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing3(string infilename, const vec &shift=vec(), double scale=1);
	int write_packing(string infilename);
	//void setup_grid(double d);

	inline bool interact(CParticle *p1,CParticle *p2, vec &force, vec &torque)const; //force from p2 on p1
	inline bool interact(CParticle *p1, GeomObject<tbox> *p2, vec &force, vec &torque)const;
	bool interact(const CParticle *p1, const GeomObject<tbox> &b, vec &force)const;//force from box on p1

	void overlappings();
	bool exist(int i);
	double t;
	double ke;//kinetic energy
	double outDt;
	GeomObject<tbox> box;
	//CRecGrid *grid;
	double maxr;
	vec G;
 	private:
	double epsFreeze;
	};

CSys::~CSys(){
	for(int i=0; i<particles.size(); i++){
		delete particles.at(i);
		}
	//delete grid;
	}

//void CSys::setup_grid(double _d){
		//grid=new CRecGrid(box.corner, box.L, _d*maxr);
		//}

bool CSys::add(CParticle *p){
	vec force(0.0);
	//if(interact(p, box, force)){
		//ERROR("The particle is intially within the box: "<<p->x(0));
		//return false;
		//}
	if(maxr<p->radius)maxr=p->radius;
	particles.push_back(p);
	//assert(grid);
	//grid->add(p);
	return true;
	}

void CSys::calForces(){
	vector<CParticle *>::iterator it1, it2;
	//reset forces
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		(*it1)->forces=G*((*it1)->get_mass());
		}
	//calculation of interaactions
	vec force;
	vec torque;
	int temp=0;
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		temp++;
		for(it2=it1+1; it2!=particles.end(); ++it2){
			if((*it1)->frozen && (*it2)->frozen)continue;
			force=0;
			if(interact(*it1, *it2, force, torque)){
				(*it1)->addforce(force);
				(*it1)->addtorque(force);
				
				(*it2)->addforce(force*-1.0);
				(*it2)->addtorque(force*-1.0);
				}
			}
		//the walls
		force=0;
		if(interact(*it1, &box, force, torque)){
			(*it1)->addforce(force);
			(*it1)->addforce(torque);
			}
		}
	
	};

void CSys::forward(double dt){
	static int count=0, outN=0,outPutN=outDt/dt;
	
	vector<CParticle *>::iterator it;
	for(it=particles.begin(); it!=particles.end(); ++it){
		if(!(*it)->frozen) (*it)->calPos(dt);
		//if(!it->frozen) it->x.gear_predict<4>(dt);
		}
	calForces();
	ke=0;
	for(it=particles.begin(); it!=particles.end(); ++it){
		if(!(*it)->frozen) (*it)->calVel(dt);
		//vec dA=(it->forces/it->get_mass()-it->x(2));
		//if(!it->frozen) it->x.gear_correct<5>(dt,dA);
		if(0)if((*it)->x(1).abs()< epsFreeze  && (*it)->avgforces.abs()< epsFreeze ) {
			(*it)->frozen=true;
			(*it)->material.color="0.5 0.5 0.5";
			(*it)->x(1)=0.0;
			}
		ke+=(*it)->kEnergy();
		//cout<< it->x <<"  "<<it->size<< " cir"<<endl;
		}
	if(count%outPutN==0){
			stringstream outname;
			outname<<"out"<<setw(5)<<setfill('0')<<outN;
			ofstream out(outname.str().c_str());
			for(it=particles.begin(); it!=particles.end(); ++it){
				//if(!(*it)->frozen)out<<**it<<endl;
				out<<**it<<endl;
				}
			count=0;
			outN++;
			}
	count++;
	}

void CSys::solve(double tMax, double dt){
	bool stop=false;
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
	}

int CSys::write_packing(string outfilename){
	vector<CParticle*>::iterator it;
	ofstream out(outfilename.c_str());
	for(it=particles.begin(); it!=particles.end(); ++it){
		out<<**it<<endl;
		}
	}
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
                ss>>p->x(0)(0)>>p->x(0)(1)>>p->x(0)(2)>>p->radius;
		//if(exist(p->id))continue;
		p->frozen=true;//FIXME if it is not frozen it doesnt work, why?
		p->identifier=0;
		if(particles.size()>1000)break;
		
		p->x(0)+=shift;
		p->x(0)*=scale;
		double ran=drand48();
		if(ran<0.05)p->identifier=1;
		else if(ran<0.2)p->identifier=2;
		else p->identifier=3;
		p->material.color=stringify(drand48())+stringify(drand48())+stringify(drand48());
		p->Xc=p->x(0);
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
                ss>>p->id>>p->x(0)(0)>>p->x(0)(1)>>p->x(0)(2)>>p->radius;
		//if(exist(p->id))continue;
		if(p->x(0)(2)<270 || p->x(0)(2)>1725)p->frozen=true;
		dist(0)=p->x(0)(0); dist(1)=p->x(0)(1);
		if((dist-center).abs()> 525)p->frozen=true;
		
		p->x(0)+=shift;
		p->x(0)*=scale;
		//p.radius=44;
		p->radius*=scale*1.04;
		p->x0(0)=p->x(0);
                add(p);
                }
        cerr<< "done" <<endl;
        inputFile.close();
	return particles.size();
        }

bool CSys::exist(int i){
	vector<CParticle *>::iterator it1;
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		if((*it1)->id==i)return true;
		}
		return false;
		}

void CSys::overlappings(){
	vector<CParticle *>::iterator it1, it2;
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
	for(it2=it1+1; it2!=particles.end(); ++it2){
		double d=((*it1)->x(0)-(*it2)->x(0)).abs()- (*it1)->radius - (*it2)->radius;
		if(d<0)cerr<< d<<" "<< (*it1)->id<<"  "<< (*it2)->id <<endl;
		}
		}

		}

inline bool CSys::interact(CParticle *p1, CParticle *p2, vec &force, vec &torque)const{

	vector<COverlapping> overlaps;
	COverlapping::overlaps(overlaps, (GeomObjectBase*)p1, (GeomObjectBase*)p2);
	//cerr<< overlaps.size()<<endl;
	vec dv=p1->x(1)-p2->x(1);
	double proj;
	if(overlaps.size()==0)return false;
	for(int i=0; i<overlaps.size(); i++){
		proj=dv*overlaps.at(i).dx;
		if(proj>0){
			force-=p1->material.stiffness1*overlaps.at(i).dx;
			torque+=cross(force, overlaps.at(i).x);
			}
		else {
			force-=p1->material.stiffness2*overlaps.at(i).dx;
			torque+=cross(force, overlaps.at(i).x);
			}
		}
/// ?????????????????? this is just a test
	//force-=friction*dv;//FIXME
/////////////////////////////
	return true;
	}

inline bool CSys::interact(CParticle *p1, GeomObject<tbox> *p2, vec &force, vec &torque)const{

	vector<COverlapping> overlaps;
	COverlapping::overlaps(overlaps, (GeomObjectBase*)p1, (GeomObjectBase*)p2);
	//cerr<< overlaps.size()<<endl;
	vec dv=p1->x(1);
	double proj;
	if(overlaps.size()==0)return false;
	for(int i=0; i<overlaps.size(); i++){
		proj=dv*overlaps.at(i).dx;
		if(proj>0){
			force-=p1->material.stiffness1*overlaps.at(i).dx;
			torque+=cross(force, overlaps.at(i).x);
			}
		else {
			force-=p1->material.stiffness2*overlaps.at(i).dx;
			torque+=cross(force, overlaps.at(i).x);
			}

		}
/// ?????????????????? this is just a test
	//force-=friction*dv;
/////////////////////////////
	return true;
	}
inline bool CSys::interact(const CParticle *p1, const GeomObject<tbox> &b, vec &force)const{
	bool overlap=false;
	double d=b.top()(0)-p1->x(0)(0)-p1->radius;

	if(d<0){
		if(p1->x(1)(0)>=0)force(0)+= (p1->material.stiffness1*d);
		else force(0)+= (p1->material.stiffness2*d);
		overlap=true;

/// ?????????????????? this is just a test
	force(0)-=friction*p1->x(1)(0);
/////////////////////////////
		}

	d=-b.corner(0)+p1->x(0)(0)-p1->radius;
	if(d<0){
		if(p1->x(1)(0)<=0)force(0)-= (p1->material.stiffness1*d);
		else force(0)-= (p1->material.stiffness2*d);
		overlap=true;
/// ?????????????????? this is just a test
	force(0)-=friction*p1->x(1)(0);
/////////////////////////////
		}

	d=b.top()(1)-p1->x(0)(1)-p1->radius;
	if(d<0){
		if(p1->x(1)(1)>=0)force(1)+= (p1->material.stiffness1*d);
		else force(1)+= (p1->material.stiffness2*d);
		overlap=true;
/// ?????????????????? this is just a test
	force(1)-=friction*p1->x(1)(1);
/////////////////////////////
		}

	d=-b.corner(1)+p1->x(0)(1)-p1->radius;
	if(d<0){
		if(p1->x(1)(1)<=0)force(1)-= (p1->material.stiffness1*d);
		else force(1)-= (p1->material.stiffness2*d);
		overlap=true;
/// ?????????????????? this is just a test
	force(1)-=friction*p1->x(1)(1);
/////////////////////////////
		}

	//for 3d
	d=b.top()(2)-p1->x(0)(2)-p1->radius;
	if(d<0){
		if(p1->x(1)(2)>=0)force(2)+= (p1->material.stiffness1*d);
		else force(2)+= (p1->material.stiffness2*d);
		overlap=true;
/// ?????????????????? this is just a test
	force(2)-=friction*p1->x(1)(2);
/////////////////////////////
		}

	d=-b.corner(2)+p1->x(0)(2)-p1->radius;
	if(d<0){
		if(p1->x(1)(2)<=0)force(2)-= (p1->material.stiffness1*d);
		else force(2)-= (p1->material.stiffness2*d);
		overlap=true;
/// ?????????????????? this is just a test
	force(2)-=friction*p1->x(1)(2);
/////////////////////////////
		}

	return overlap;
	}

#endif /* MDSYS_H */
