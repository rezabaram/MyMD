#ifndef MDSYS_H
#define MDSYS_H 
#include "common.h"
#include"celllist.h"
#include"particle.h"
#include"packing.h"
#include"interaction.h"
#include"interaction_force.h"
#include"particlecontact.h"
#include"map_asph_aspect.h"
#include"ibeta_dist.h"

#include<tr1/random>
std::tr1::ranlux64_base_01 eng;

//#define WITH_VERLET

extern MTRand rgen;

typedef CPacking<CParticle> ParticleContainer;

typedef GeomObjectBase * BasePtr;
class CSys{
	CSys();
	public:
	CSys(unsigned long maxnparticle):t(0), outDt(0.01), 
	walls(config.get_param<vec>("boxcorner"), config.get_param<vec>("boxsize"), config.get_param<string>("boundary")),
	maxr(0), maxh(0), maxv(0), G(vec(0.0)),
	maxNParticle(maxnparticle), 
	#ifdef WITH_VERLET
	verlet((&particles)), 
	#else 
	#endif

	epsFreeze(1.0e-12), outEnergy("log_energy"),
	maxRadii(0), top_v(vec(0.0,0.0,0.0)),
	celllist(CCellList<ParticleContainer, CParticle>(&walls))
	{
	TRY
	CATCH
		};
	~CSys();

	void particles_on_grid();
	void add_particle_layer(double z);

	void initialize(const CConfig &c);
	void solve();
	void forward(double &dt);
	void adapt(double &dt);
	void calForces();
	void interactions();
	inline bool interact(unsigned int i, unsigned int j)const; //force from p2 on p1
	inline bool interact(CParticle *p1,CParticle *p2)const;
	inline bool interact(CParticle *p1, BoxContainer *p2);
	vec center_of_mass()const;

	int read_packing(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing2(string infilename, const vec &shift=vec(), double scale=1);
	int read_packing3(string infilename, const vec &shift=vec(), double scale=1);
	void read_radii(vector<vec> &radii );
	void write_packing(string infilename);
	double total_volume();
	//void setup_grid(double d);

	bool add(CParticle *p);
	void remove(ParticleContainer &packing, ParticleContainer::iterator &it);

	inline bool exist(int i);

	void output(string outname);
	void output(ostream &out=std::cout);

	double t, tMax, dt,DT, outDt, outStart, outEnd;

	//contains the pointers to the particles
	ParticleContainer particles;
	

	BoxContainer walls;
	CPlane *sp;
	double minr, maxr, maxh, maxv;
	#ifdef WITH_VERLET
	double verlet_factor;
	#endif
	//bool verlet_need_update;
	vec G;
	const unsigned maxNParticle;
	double fluiddampping;
	#ifdef WITH_VERLET
	CVerletManager<CParticle> verlet;
	#endif
	double Energy, rEnergy, pEnergy, kEnergy;
 	private:
	bool do_read_radii, softwalls, spherize_on;
	string out_name;
	double epsFreeze;
	ofstream outEnergy;
	vector<vec> radii;
	ifstream inputRadii;
	double maxRadii;
	vec top_v;

	CCellList<ParticleContainer, CParticle> celllist;
	
	};

CSys::~CSys(){
TRY
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

double TruncGaussRand(double r, double dr=0.0){
	if(dr<1e-6) return r;
	double x=rgen.randNorm(r, dr) ;
	if(x<r-2*dr or x> r+2*dr)return TruncGaussRand(r, dr);//trying until finding in range (r-dr, r+dr)
	return x;
	}

void CSys::initialize(const CConfig &config){
TRY
	G=config.get_param<vec>("Gravity");
	out_name=config.get_param<string>("output");
	outDt=config.get_param<double>("outDt");
	outStart=config.get_param<double>("outStart");
	outEnd=config.get_param<double>("outEnd");
	fluiddampping=config.get_param<double>("fluiddampping");
	string particleType=config.get_param<string>("particleType");
	string init_method=config.get_param<string>("initialization");
	softwalls=config.get_param<bool>("softwalls");
	spherize_on=config.get_param<bool>("spherize_on");

	double dr=config.get_param<double>("particleSizeWidth");
	DisBetaDistribution ibeta_dist(3,3,dr);



	if(init_method=="restart"){
		particles.parse(config.get_param<string>("input"));
		maxRadii=particles.maxr;
		}
	else if(init_method=="generate"){
		if(particleType=="general" or particleType=="gen1" or particleType=="gen2"
					   or particleType=="gen3" or particleType=="gen4"){
			string fileRadii=config.get_param<string>("radii");
			inputRadii.open(fileRadii.c_str());
			ERROR(!inputRadii.good(), "Unable to open input file: "+fileRadii );
			read_radii(radii);
			}
		else if(particleType=="prolate"){
		//calculating the maximum radius
			double ee=get_aspect_prolate(config.get_param<double>("asphericity"));
			cerr<<config.get_param<double>("asphericity") <<endl;
			double r=config.get_param<double>("particleSize");
			double a =r/pow(ee,1./3.);
			double b =a;
			double c =ee*a;
			maxRadii=max(r,max(a,max(b,c)));
			radii.push_back(vec(a, b, c));

			}
		else if(particleType=="sandstone"){
			
			double rmin=config.get_param<double>("rmin");
			double rmax=config.get_param<double>("rmax");
			std::tr1::uniform_real<double> unif(rmin, rmax);
			
			for(int i=0; i<10000;i++){
				double a = unif(eng);
				double b = unif(eng);
				double c = unif(eng);
				maxRadii=max(a,max(b,c));
				radii.push_back(vec(a, b, c));
				}
			}
		else if(particleType=="oblate"){
			
			double ee=get_aspect_oblate(config.get_param<double>("asphericity"));
			double r=config.get_param<double>("particleSize");
			double a =r/pow(ee,1./3.);
			double b =a;
			double c =ee*a;
			maxRadii=max(r,max(a,max(b,c)));
			radii.push_back(vec(a, b, c));
			}
		else if(particleType=="abc"){
			double eta0=config.get_param<double>("eta");
			double etaW=config.get_param<double>("etaWidth");
			double xi0=config.get_param<double>("xi");
			double xiW=config.get_param<double>("xiWidth");
			double r0=config.get_param<double>("particleSize");
			for(int i=0; i<10000;i++){
				double eta=eta0*TruncGaussRand(1, etaW);
				double xi =xi0*TruncGaussRand(1, xiW);
				double r=r0*ibeta_dist.rnd();
				double a =r*pow(eta,1./3.)/pow(xi,1./3);
				double b =r/(pow(eta,2./3.)*pow(xi, 1/3.));
				double c =r*xi*pow(eta/xi,1./3.);
				maxRadii=max(r,max(a,max(b,c)));
				radii.push_back(vec(a, b, c));
				}
			}
		else
			ERROR(1, "Unknown particle type: "+particleType);
		}
	else{
		ERROR(1, "Unknown initialization method: "+init_method);
		}

	#ifdef WITH_VERLET
	verlet.set_distance(particles.maxr*config.get_param<double>("verletfactor"));
	verlet.update();
	#else
	celllist.setup(2.0*maxRadii);
	celllist.build(particles);
	#endif

	 //add_particle_layer(1.02*maxRadii);

	cerr<< "Number of Particles: "<<particles.size() <<endl;

	tMax=config.get_param<double>("maxTime");
	DT=config.get_param<double>("timeStep");
	dt=DT;

        minr=config.get_param<double>("particleSize");
CATCH
	}

double CSys::total_volume(){
	ParticleContainer::iterator it;
	double v=0;
	for(it=particles.begin(); it!=particles.end(); ++it){
		v+=(*it)->shape->vol();
		}
	return v;
	}


void CSys::remove(ParticleContainer &packing, ParticleContainer::iterator &it){
TRY
        delete (*it);
       packing.erase(it);
CATCH
        }


bool CSys::add(CParticle *p){
TRY
	//ERROR(particles.size()==maxNParticle, "Reached max number of particles.");
	//vec force(0.0);
	//if(interact(p, box, force)){
		//ERROR("The particle is intially within the box: "<<p->x(0));
		//return false;
		//}
	if(maxr<p->shape->radius)maxr=p->shape->radius;
	if(p->top()>maxh) maxh=p->top();

	particles.push_back(p);
	particles.back()->id=particles.TotalParticlesN;
	 ++(particles.TotalParticlesN);

	#ifdef WITH_VERLET
	if(verlet.add_particle(p)){}
	#else
	celllist.add(p);
	#endif

	//setup_verlet(particles.back());
	//assert(grid);
	//grid->add(p);
	return true;
CATCH
	}

void CSys::calForces(){
TRY
//FROMTIME
	//FIXME make it dimensionless

	ParticleContainer::iterator it1, it2, ittemp;
	#ifdef WITH_VERLET
	CVerletList<CParticle>::iterator neigh;
	#endif
	//reset forces
	vec cm=center_of_mass();
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		(*it1)->reset_forces(G*((*it1)->get_mass())-fluiddampping*G.abs()*(*it1)->get_mass()*(*it1)->x(1));//gravity plus damping (coef of dumping is ad hoc)
		(*it1)->reset_torques(vec(0.0));
		}

	//interactions
	#ifdef WITH_VERLET
	verlet.update();
	#endif
	for(it1=particles.begin(); it1!=particles.end(); ++it1){

		#ifdef WITH_VERLET
		ERROR(!(*it1)->vlist.set,"the verlet list was not constructed properly");

		assert(*it1==(*it1)->vlist.self_p);
		for(neigh=(*it1)->vlist.begin(); neigh!=(*it1)->vlist.end(); ++neigh){
			if(Test::interact(*it1,(*neigh).first)){
				//this to calculate contact duration
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
		#endif

		//the walls
		if(interact(*it1, &walls)){  }
		}
	#ifndef WITH_VERLET
	celllist.build(particles);
	celllist.interact();
	#endif
//TOTIME
CATCH
};



void CSys::output(ostream &out){
TRY
	walls.print(out);
	ParticleContainer::iterator it;
	for(it=particles.begin(); it!=particles.end(); ++it){
		out<<**it<<endl;
		}

CATCH
	}

void CSys::output(string outname){
TRY
	ofstream out(outname.c_str());
	output(out);
CATCH
	}

void CSys::forward(double &dt){
TRY
       ParticleContainer::iterator it;
       CParticle *p1;
       for(it=particles.begin(); it!=particles.end(); ++it){
               p1=*it;
               //if((vec2d(p1->x(0)(0),p1->x(0)(1))-vec2d(0.5, 0.5)).abs()>0.7*(1.2-p1->x(0)(2)))p1->expired=true;
               //if(abs(p1->x(0)(0)-0.5)>0.4*(1.2-p1->x(0)(2)))p1->expired=true;
	       //if(p1->x(0)(2)<config.get_param<double>("particleSize")/2.)p1->frozen=true;
               if((*it)->expired){
                       remove(particles, it);//As the side effect, "it" is set to next value
                       if(it==particles.end())break;
                       }
               }
       if(config.get_param<string>("initialization")=="generate" and maxh< 1.+2*maxRadii ) 
               {
               add_particle_layer(maxh+1.02*maxRadii);
               maxh=0;
               }


	static int count=0, outN=0,outPutN=outDt/DT;
	static ofstream out;




	//this is for a messure of performance
	static double starttime=clock();
	if((10*count)%outPutN==0)
		cout<<(clock()-starttime)/CLOCKS_PER_SEC<< "   "<<t<<endl;

	if(count%outPutN==0 and t>=outStart and t<=outEnd){
			stringstream outstream;
			outstream<<out_name<<setw(5)<<setfill('0')<<outN;
			output(outstream.str());
			count=0;
			outN++;
			Energy=rEnergy+kEnergy+pEnergy;
			outEnergy<<setprecision(14)<<t<<"  "<<Energy<<"  "<<kEnergy<<"  "<<pEnergy<<"  "<<rEnergy <<endl;
			rEnergy=0; pEnergy=0; kEnergy=0; Energy=0;
			//for relaxation
			if(t>2 and G.abs()>1){
					G*=0.9;
					dt*=1.08;	
					cerr<<"t: "<<t<<" G: "<<G<<" dt: "<<dt<<endl;
					}
			if(spherize_on)for(it=particles.begin(); it!=particles.end(); ++it){
				(*it)->shape->spherize();
				}
			}
	//bool allforwarded=false;
	maxh=0;
	
	for(it=particles.begin(); it!=particles.end(); ++it){
		if(!(*it)->frozen) 
			(*it)->calPos(dt);

		if((*it)->top()>maxh) {
				maxh=(*it)->top();
				top_v=(*it)->x(1);
				}
		//if(!it->frozen) it->x.gear_predict<4>(dt);
		}


	calForces();
//	if(!allforwarded)foward(dt/2.0, 2);

	Energy=0.0, rEnergy=0, pEnergy=0, kEnergy=0;
	double vtemp;
	maxv=0;
	for(it=particles.begin(); it!=particles.end(); ++it){
	//	if(!(*it)->frozen) 
		(*it)->calVel(dt);
		vtemp=(*it)->x(1).abs()+(*it)->w(1).abs()*(*it)->shape->radius;
		if(vtemp>maxv) 
			maxv=vtemp;

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

void CSys::adapt(double &dt){
	//FIXME Just for trying. 
	//	I dont think adaptive time step can be done without considering
	//	corresponding changes in the intergrator.
	cerr<< dt/DT <<endl;
	if(dt*maxv < minr*1e-3 and dt<DT*100)dt*=1.02;
	if(dt*maxv > minr*1e-3 and dt>DT/10)dt*=0.98;

	}

void CSys::solve(){
	try{
	bool stop=false;
	while(true){
		//dt=(double)((int)t+1)*dt0;
		if(t+dt>tMax ){//to stop exactly at tMax
			dt=tMax-t;
			stop=true;
			}

		//calForces();
		forward(dt);
		t+=dt;
		if(stop or (t>1 and kEnergy<1e-8) ){
			output(out_name+"end");
			cerr<<"Relaxation criterion reached at time="<<t<<": KE= "<<kEnergy<< " < 1e-8"<<endl;
			break;
			}
		//adapt(dt);
		}
	}catch(CException e){
		ERROR(1,"Some error in the solver at t= "+ stringify(t)+"\n\tfrom "+e.where());
		}
	catch(...){
		ERROR(1,"Unknown error at t= "+ stringify(t));
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

inline bool CSys::interact(CParticle *p1, BoxContainer *p2){
TRY
	ShapeContact &overlaps=p1->vlist[p2];
	overlaps.clear();
	CInteraction::overlaps(&overlaps, p1->shape, (GeomObjectBase*)p2);

	static vec dv, r1, force, torque, vt, vn;
	if(overlaps.size()==0)return false;
	for(size_t i=0; i<overlaps.size(); i++){

		r1=overlaps(i).x-p1->x(0);
		dv=p1->x(1)+cross(p1->w(1), r1);

		force=Test::contactForce(overlaps(i), dv, p1->material, 1);
		p1->addforce(force);
		
		if(softwalls)continue;
		torque=cross(r1, force);
		p1->addtorque(torque);

		}

	return true;
CATCH
	}


double rand_aspect_ratio(double asphericity, double asphericityWidth){

	/*
	double temp=asphericity-2*asphericityWidth;
	int randtry=0;
	while(randtry<10 and (temp<asphericity-asphericityWidth or temp>asphericity+asphericityWidth) ){
		temp=rgen.randNorm(asphericity, asphericityWidth);
		++randtry;
		}
	*/

	//FIXME 
	//if large width, some aspect ratios can be too big (then you need to make the times step smaller)
	//or you may introduce a cut off
	
	if(asphericityWidth<1e-5)return exp(asphericity);
	else return exp(rgen.randNorm(asphericity, asphericityWidth) );
	}

void spheroid(double &a, double &b, double &c, double asph, double w){
		double ee;
		ee=rand_aspect_ratio(asph, w);
		a =1/pow(ee,1./3.);
		b =a;//*rgen();
		c =ee*a;//*rgen();
}
void CSys::read_radii(vector<vec> &radii){
	
	double size=config.get_param<double>("particleSize");
	double dr=size*config.get_param<double>("particleSizeWidth");
	DisBetaDistribution ibeta_dist(3,3,dr);



	string line;
	//Parse the line
	double a, b, c;
	while(getline(inputRadii,line)){
	double r=size*ibeta_dist.rnd()*pow(5./3.*M_PI,1./3.);//note 4/3 pi a b c=1 
		stringstream ss(line);
		ss>>a>>b>>c;
		radii.push_back(vec(r*a,r*b,r*c));
		maxRadii=max(maxRadii,r*max(a, max(b,c)));
		}
}


void CSys::add_particle_layer(double z){ 
	double size=config.get_param<double>("particleSize");
	vec x(0.0, 0.0, .0);

	
	//double asphericity=config.get_param<double>("asphericity");
	//double asphericityWidth=config.get_param<double>("asphericityWidth");

	double xtemp=0, ytemp=0;

	static unsigned int nRadii=0;
	for(int i=0;i<celllist.nx;i++){
		for(int j=0;j<celllist.ny;j++){
		if(particles.size()>=maxNParticle)break;
		double a, b, c;
		//spheroid(a, b, c, asphericity, asphericityWidth);
		ERROR(nRadii>=radii.size(), "List of radii doesn't have enough entries");

		int randn=rgen.rand(radii.size());
		if(config.get_param<string> ("particleType") == "gen1") randn=1;
		if(config.get_param<string> ("particleType") == "gen2") randn=2;
		if(config.get_param<string> ("particleType") == "gen3") randn=3;
		if(config.get_param<string> ("particleType") == "gen4") randn=4;
		vec abc=radii.at(randn);

		a=abc(0);b=abc(1);c=abc(2);


		xtemp=(i+0.5)*celllist.dx;
		ytemp= (j+0.5)*celllist.dy;

/*
		i+=2.1*maxRadii;
		if(j<maxRadii)j=1.5*maxRadii;
		if(k<maxRadii)k=1.5*maxRadii;
		if(i>1-1.2*maxRadii){
			i=1.1*maxRadii;
			j+=2.1*maxRadii;
			}
		if(j>1-1.2*maxRadii){
			break;
			}
*/

		x(0)=xtemp+size*rgen()/5;
		x(1)=ytemp+size*rgen()/5; 
		x(2)=z+size*rgen()/5; 
		double alpha=rgen()*M_PI;
                double beta=rgen()*M_PI;
               	double phi=rgen()*M_PI;
               	Quaternion q=Quaternion(cos(alpha),sin(alpha),0,0)*Quaternion(cos(beta),0,0,sin(beta))*Quaternion(cos(phi),0,sin(phi),0 );
		//CParticle *p = new CParticle(CSphere(x,size*(1-0.0*rgen())));
		//CEllipsoid E(x, 1-0.0*rgen(), 1-0.0*rgen(),1-0.0*rgen(), size*(1+0.0*rgen()));
		//CParticle *p = new CParticle(E);
		CSphere E1(x,maxRadii);
		//CEllipsoid E2(x, 1, 1, 1, size, q);

		//to implement constant volume (4/3 Pi r^3) while changing the shape
		CEllipsoid E2(x, a,b,c, q);
		CParticle *p = new CParticle(E2);
		p->w(1)(0)=5.0*(1-2*rgen());
		p->w(1)(1)=5.0*(1-2*rgen());

		p->x(1)(0)=0.3*(1-2*rgen());
		p->x(1)(1)=0.3*(1-2*rgen());
		p->x(1)(2)=top_v(2)+0.3*(1-2*rgen());
		add(p);
		
		}
		}
	}

void CSys::particles_on_grid(){ 
TRY
	double size=config.get_param<double>("particleSize");
	vec x(0.0, 0.0, .0);

	double k=0;
	
	//double asphericity=config.get_param<double>("asphericity");
	//double asphericityWidth=config.get_param<double>("asphericityWidth");

	double i=0, j=0;

	unsigned int nRadii=0;
	read_radii(radii);
	while(particles.size()<maxNParticle){
		if(particles.size()==maxNParticle)break;
		double a, b, c;
		double r=size;
		//spheroid(a, b, c, asphericity, asphericityWidth);
		vec abc=radii.at(nRadii);
		a=r*abc(0);b=r*abc(1);c=r*abc(2);
		++nRadii;

		r=max(r,max(a,max(b,c)));
		i+=2.1*r;
		if(j<r)j=1.5*r;
		if(k<r)k=1.5*r;
		if(i>1-1.2*r){
			i=1.1*r;
			j+=2.1*r;
			}
		if(j>1-1.2*r){
			i=2.1*r;
			j=2.1*r;
			k+=2.1*r;
			}

		x(0)=i+size*rgen()/10;
		x(1)=j+size*rgen()/10; 
		x(2)=k+size*rgen()/10; 
		double alpha=rgen()*M_PI;
		Quaternion q=Quaternion(cos(alpha),sin(alpha),0,0)*Quaternion(cos(alpha),0,0,sin(alpha) );
		//CParticle *p = new CParticle(CSphere(x,size*(1-0.0*rgen())));
		//CEllipsoid E(x, 1-0.0*rgen(), 1-0.0*rgen(),1-0.0*rgen(), size*(1+0.0*rgen()));
		//CParticle *p = new CParticle(E);
		CSphere E1(x,r);
		//CEllipsoid E2(x, 1, 1, 1, size, q);

		//to implement constant volume (4/3 Pi r^3) while changing the shape
		CEllipsoid E2(x, a,b,c);
		CParticle *p = new CParticle(E2);
		p->w(1)(0)=5.0*(1-2*rgen());
		p->w(1)(1)=5.0*(1-2*rgen());

		p->x(1)(0)=0.3*(1-2*rgen());
		p->x(1)(1)=0.3*(1-2*rgen());
		p->x(1)(2)=0.3*(1-2*rgen());
		add(p);
		
		}
CATCH
}

//void CSys::setup_grid(double _d){
		//grid=new CRecGrid(box.corner, box.L, _d*maxr);
		//}
#endif /* MDSYS_H */
