#ifndef MDSYS_H
#define MDSYS_H 
#include "common.h"
#include"grid.h"
#include"particle.h"
#include"interaction.h"
#include"particlecontact.h"

#define VERLET

typedef vector<CParticle *> ParticleContainer;

typedef GeomObjectBase * BasePtr;
class CSys{
	CSys();
	public:
	CSys(unsigned long maxnparticle):t(0), outDt(0.01), box(vec(0.0), vec(1.0)), maxr(0), maxh(0), G(vec(0.0)), 
		maxNParticle(maxnparticle), verlet_need_update(true),epsFreeze(1.0e-12), outEnergy("log_energy"){
	TRY
	#ifndef VERLET
		pairs=new ParticleContactHolder<CParticle> *[maxNParticle*maxNParticle];
		for(unsigned int i=0; i<maxNParticle*maxNParticle; i++) pairs[i]=NULL;
	#endif
	CATCH
		};
	~CSys();

	void particles_on_grid();
	void initialize(const CConfig &c);
	void solve();
	void forward(double dt);
	void calForces();
	void interactions();
	inline bool interact(unsigned int i, unsigned int j)const; //force from p2 on p1
	inline bool interact(CParticle *p1,CParticle *p2)const;
	inline bool interact(CParticle *p1, GeomObject<tbox> *p2)const;
	vec contactForce(const Contact &c, const vec &dv, CProperty &m, double tmp=1)const;
	vec center_of_mass()const;

	int read_packing(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing2(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing3(string infilename, const vec &shift=vec(), double scale=1);
	void write_packing(string infilename);
	double total_volume();
	//void setup_grid(double d);

	bool add(CParticle *p);
	double min_distance(CParticle *p1, CParticle *p2)const;
	void setup_verlet(CParticle *p);
	void print_verlet(ostream &out);
	void update_verlet();
	inline bool exist(int i);

	double t, tMax, dt, outDt;
	ParticleContainer particles;

	
	#ifndef VERLET
	ParticleContactHolder<CParticle>  **pairs;
	ParticleContactHolder<CParticle>  *  &pair(size_t i, size_t j)const{
		return pairs[j+i*particles.size()];
			}
	#endif

	GeomObject<tbox> box;
	CPlane *sp;
	//CRecGrid *grid;
	double maxr, maxh;
	double verlet_distance, min_verlet_distance, max_verlet_distance, verlet_factor;
	vec G;
	const unsigned maxNParticle;
	bool verlet_need_update;
 	private:
	double epsFreeze;
	ofstream outEnergy;
	};

CSys::~CSys(){
TRY
	#ifndef VERLET
	for(size_t i=0; i<particles.size(); i++)
	for(size_t j=0; j<particles.size(); j++){
		delete pair(i,j);
		}
	delete [] pairs;
	#endif
	ParticleContainer::iterator it;
	for(it=particles.begin(); it!=particles.end(); ++it){
		delete (*it);
		}
CATCH
	}

vec CSys::center_of_mass()const{
TRY
	vec cm(0,0,0);
	double M=0;
	ParticleContainer::const_iterator it;
	for(it=particles.begin(); it!=particles.end(); ++it){
		
		M+=(*it)->get_mass();
		cm+=(*it)->get_mass()*(*it)->x(0);
		}
		cm/=M;
	return cm;
CATCH
	}

void CSys::initialize(const CConfig &c){
TRY
	outDt=c.get_param<double>("outDt");
	maxr=c.get_param<double>("particleSize");

	G=c.get_param<vec>("Gravity");

	cerr<< 2 <<endl;
	particles_on_grid();
	cerr<< "Number of Particles: "<<particles.size() <<endl;

	#ifdef VERLET
	verlet_factor=c.get_param<double>("verletfactor");
	verlet_distance=verlet_factor*maxr;
	min_verlet_distance=verlet_distance/5;
	max_verlet_distance=verlet_distance*5;
	update_verlet();
	#else
	//create for each pair an object for the contact
	for(unsigned int i=0; i<particles.size(); i++){
		for(unsigned int j=0; j<particles.size(); j++){
			pair(i,j)=new ParticleContactHolder<CParticle> (particles.at(i), particles.at(j));
			}
		}
	#endif

	tMax=config.get_param<double>("maxTime");
	dt=config.get_param<double>("timeStep");
CATCH
	}

double CSys::total_volume(){
	vector<CParticle *>::iterator it;
	double v=0;
	for(it=particles.begin(); it!=particles.end(); ++it){
		v+=(*it)->shape->vol();
		}
	return v;
	}

//void CSys::setup_grid(double _d){
		//grid=new CRecGrid(box.corner, box.L, _d*maxr);
		//}

//minimum distances of the surfaces of the particles (putting them in spherical shells)
double CSys::min_distance(CParticle *p1, CParticle *p2)const{
TRY
	return (p1->x(0)-p2->x(0)).abs() - p1->shape->radius -p2->shape->radius;
CATCH
	}

void CSys::print_verlet(ostream &out){
	vector<CParticle *>::iterator it;
	CVerlet<CParticle >::iterator neigh;
	for(it=particles.begin(); it!=particles.end(); ++it){
		out<< (*it)->id <<": ";
		for(neigh=(*it)->vlist.begin(); neigh!=(*it)->vlist.end(); ++neigh ){
			out<< (*neigh).first->id <<" ";
			}
			out<< endl;
		}
	}
//construct the verlet list of all particles
void CSys::update_verlet(){
	verlet_distance*=0.99;
	if(verlet_distance<min_verlet_distance)verlet_distance=min_verlet_distance;

	vector<CParticle *>::iterator it;
	for(it=particles.begin(); it!=particles.end(); ++it){
		if(verlet_need_update==false and ((*it)->x(0)-(*it)->vlist.x).abs() > verlet_distance/2.0-epsilon) verlet_need_update=true;
		}

	if(!verlet_need_update) return;

	for(it=particles.begin(); it!=particles.end(); it++){
		(*it)->vlistold=(*it)->vlist;
		(*it)->vlist.clear();
		setup_verlet(*it);
		}

	verlet_distance*=1.5;
	if(verlet_distance>max_verlet_distance)verlet_distance=max_verlet_distance;

	verlet_need_update=false;
	//cerr<< "Verlet updated at: "<<t <<"\t verlet distance: "<<verlet_distance<<endl;
	}

//construct the verlet list of one particle 
void CSys::setup_verlet(CParticle *p){
TRY
	vector<CParticle *>::iterator it;
	CVerlet<CParticle>::iterator v_it_old;
	for(it=particles.begin(); (*it)!=p; it++){ //checking particles before in the list
	//for(it=particles.begin(); it!=particles.end(); it++){ //checking particles before in the list
		if(min_distance(p, *it) < verlet_distance){
			v_it_old=p->vlistold.find(*it);
			if(v_it_old!= p->vlistold.end())
				p->vlist.insert(*v_it_old);
				else
				p->vlist.add((*it));
				}
		}
		assert((*it)->id == p->id);
	p->vlist.x=p->x(0);//save the position at which the list has been updated
	p->vlist.set=true;
CATCH
	}

bool CSys::add(CParticle *p){
TRY
	ERROR(particles.size()==maxNParticle, "Reached max number of particles.");
	vec force(0.0);
	//if(interact(p, box, force)){
		//ERROR("The particle is intially within the box: "<<p->x(0));
		//return false;
		//}
	if(maxr<p->shape->radius)maxr=p->shape->radius;
	if(maxh<p->x(0)(2))maxh=p->x(0)(2);

	particles.push_back(p);
	particles.back()->id=particles.size()-1;

	//setup_verlet(particles.back());
	//assert(grid);
	//grid->add(p);
	return true;
CATCH
	}

void CSys::calForces(){
TRY
//FROMTIME
	ParticleContainer::iterator it1, it2, ittemp;
	CVerlet<CParticle>::iterator neigh;
	//reset forces
	vec cm=center_of_mass();
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		(*it1)->forces=G*((*it1)->get_mass())-2.0*(*it1)->get_mass()*(*it1)->x(1);//gravity plus damping
		//(*it1)->forces=-G.abs()*((**it1).x(0)-cm)*((*it1)->get_mass())-2.0*(*it1)->get_mass()*(*it1)->x(1);//gravity plus damping
		//(*it1)->forces=-50*G.abs()*((**it1).x(0)-cm)*((*it1)->get_mass())+G*(*it1)->get_mass();
		(*it1)->torques=0.0;
		}

	//interactions
	#ifdef VERLET
	update_verlet();
	#endif
	for(it1=particles.begin(); it1!=particles.end(); ++it1){

		#ifdef VERLET
		ERROR(!(*it1)->vlist.set,"the verlet list was not constructed properly");

		assert(*it1==(*it1)->vlist.self_p);
		for(neigh=(*it1)->vlist.begin(); neigh!=(*it1)->vlist.end(); ++neigh){
			if(interact(*it1,(*neigh).first)){
				if(!neigh->second.in_contact){
					neigh->second.in_contact=true;
					neigh->second.col_time=t;
					}
				}
			else{
				if(neigh->second.in_contact){
					neigh->second.in_contact=false;
					//cerr<< (t-neigh->second.col_time)/dt <<endl;
					}
				}
			}
		#else
		it2=it1; ++it2;
		for(it2=particles.begin(); it2!=it1; ++it2){//without verlet
			if(interact(*it1,*it2)){ }
			}
		#endif

		//the walls
		if(interact(*it1, &box)){  }
		}
//TOTIME
CATCH
};

inline
vec CSys::contactForce(const Contact &c, const vec &dv, CProperty &m, double tmp)const{
TRY
	double proj=(dv*c.n);
	double ksi=c.dx_n;
	

	ksi=(m.stiffness*ksi+tmp*m.damping*proj)*sqrt(ksi); 
	if(ksi<0)ksi=0;//to eliminate artifical attractions
	vec fn=-ksi*c.n;//normal force
	vec ft=-tmp*m.friction*(fn.abs())*((dv - proj*c.n));//dynamic frictions
	//cerr<< ft.abs()/(fn+ft).abs()<<"   "<<ft.abs()<<"  "<<fn.abs()<<endl;

	static int i2=0;
	i2++;
	//if(i2==160)exit(0);
	//cerr<< t<<"  "<<setprecision(14)<<c.x<<"  "<<c.n<<"  "<< c.dx_n << "  dv="<<dv<<"  Fn="<<fn<<endl;
	//cerr<< ft.abs()/fn.abs() <<endl;

	//ERROR((fn+ft).abs()>20, "divergence in force");
	return fn+ft;//visco-elastic Hertz law
CATCH
	}

inline bool CSys::interact(CParticle *p1,CParticle *p2)const{
TRY
	//CParticle *p1=particles.at(i);
	//CParticle *p2=particles.at(j);
#ifdef VERLET
	ShapeContact &overlaps=p1->vlist[p2];
	//ShapeContact overlaps;
#else
	ShapeContact &overlaps=*pair(p1->id,p2->id);
	//ShapeContact overlaps;
#endif
	overlaps.clear();
	CInteraction::overlaps(&overlaps, p1->shape, p2->shape);
	static vec r1, r2, v1, v2, dv, force, torque(0.0);
	//static double proj, ksi;
	if(overlaps.size()==0)return false;
	for(size_t ii=0; ii<overlaps.size(); ii++){
		r1=overlaps(ii).x-p1->x(0);
		r2=overlaps(ii).x-p2->x(0);
		v1=p1->x(1)+cross(p1->w(1), r1);
		v2=p2->x(1)+cross(p2->w(1), r2);
		dv=v1-v2;

		force=contactForce(overlaps(ii), dv, p1->material);
		p1->shape->fixToBody(HomVec(overlaps(ii).x,1));
		//cerr<< (p2->x(0)-p1->x(0)).normalized() <<endl;


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


	if(0)if(t>0.5){
		print_verlet(cout);
		exit(0);
		}

	ParticleContainer::iterator it;

	static double starttime=clock();
	if((10*count)%outPutN==0)
		cout<<(clock()-starttime)/CLOCKS_PER_SEC<< "   "<<t<<endl;

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
	maxh=0;
	for(it=particles.begin(); it!=particles.end(); ++it){
	//	if(!(*it)->frozen) 
		(*it)->calPos(dt);
		if((*it)->x(0)(2)>maxh) maxh=(*it)->x(0)(2);
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

void CSys::solve(){
	try{
	bool stop=false;
	while(true){
		//dt=(double)((int)t+1)*dt0;
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
		double ran=rgen();
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
		double ran=rgen();
		if(ran<0.05)p->shape->identifier=1;
		else if(ran<0.2)p->shape->identifier=2;
		else p->shape->identifier=3;
		p->material.color=stringify(rgen())+stringify(rgen())+stringify(rgen());
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

		force=contactForce(overlaps(i), dv, p1->material, 0);
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

double rand_aspect_ratio(double asphericity, double asphericityWidth){

	double temp=asphericity-2*asphericityWidth;
	int randtry=0;
	while(randtry<10 and (temp<asphericity-asphericityWidth or temp>asphericity+asphericityWidth) ){
		temp=rgen.randNorm(asphericity, asphericityWidth);
		++randtry;
		}

	if(asphericityWidth<1e-3)return exp(asphericity);
	else return exp(rgen.randNorm(asphericity, asphericityWidth) );
	}

void CSys::particles_on_grid(){ 
TRY
	double size=maxr;
	double margin=4.3*size;
	vec x(0.0, 0.0, .0);

	double k=0;

	
	double ee;
	double asphericity=config.get_param<double>("asphericity");
	double asphericityWidth=config.get_param<double>("asphericityWidth");

	while(particles.size()<maxNParticle){
		k=maxh+margin;
		if(k==1+10*margin)break;
		for(double i=margin/2; i<1-margin/2; i+=margin){
			for(double j=1-margin/2; j>margin/2; j-=margin){
				if(particles.size()==maxNParticle)break;
				ee=rand_aspect_ratio(asphericity, asphericityWidth);
				double r=size;//*(1-0.1*rgen());
				double a =r/pow(ee,1./3.);
				double b =a;//*rgen();
				double c =ee*a;//*rgen();

				x(1)=i+size*rgen()/10; 
				x(0)=j+size*rgen()/10;
				x(2)=k+size*rgen()/10; 
				double alpha=rgen()*M_PI;
				Quaternion q=Quaternion(cos(alpha),sin(alpha),0,0)*Quaternion(cos(alpha),0,0,sin(alpha) );
				//CParticle *p = new CParticle(GeomObject<tsphere>(x,size*(1-0.0*rgen())));
				//GeomObject<tellipsoid> E(x, 1-0.0*rgen(), 1-0.0*rgen(),1-0.0*rgen(), size*(1+0.0*rgen()));
				//CParticle *p = new CParticle(E);
				GeomObject<tsphere> E1(x,r);
				//GeomObject<tellipsoid> E2(x, 1, 1, 1, size, q);

				//to implement constant volume (4/3 Pi r^3) while changing the shape
				GeomObject<tellipsoid> E2(x, a,b,c);
				CParticle *p = new CParticle(E2);
				p->w(1)(1)=10*(1-2*rgen());
				p->w(1)(0)=10*(1-2*rgen());
				p->x(1)(1)=2*(1-2*rgen());
				p->x(1)(0)=2*(1-2*rgen());
				add(p);
				
				}
			}
		}
CATCH
}
#endif /* MDSYS_H */
