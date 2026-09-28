// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef MDSYS_H
#define MDSYS_H 
#include "common.h"
#include "version.h"
#include"config.h"
#include"celllist.h"
#include"packing.h"
#include"interaction.h"
#include"interaction_force.h"
#include"particlecontact.h"
#include"map_asph_aspect.h"
#include"ibeta_dist.h"
#include"size_dist.h"

extern CConfig config;

#include<random>
std::ranlux48_base eng;      // see the note in include/size_dist.h

extern MTRand rgen;

typedef CPacking<CParticle> ParticleContainer;


typedef GeomObjectBase * BasePtr;
class CSys{
	CSys();
	public:
	// The initialiser list follows the order in which the members are
	// declared.  It did not before, which is harmless here (nothing depends on
	// another member) but is the kind of thing that hides a real bug the day
	// something does.
	CSys(unsigned long maxnparticle=100000000):t(0), outDt(0.01)
	,walls(config.get_param<vec>("boxcorner"), config.get_param<vec>("boxsize"), config.get_param<string>("boundary"))
	,celllist(CCellList<ParticleContainer, CParticle>(&walls))
	,maxr(0), maxh(0), maxv(0), G(vec(0.0))
	,maxNParticle(maxnparticle)
	,epsFreeze(1.0e-12)
	,maxRadii(0), top_v(vec(0.0,0.0,0.0))
	,rainAccum(0), rainBlockedTime(0)
	,relaxation(true), rainFullWarned(false)
	,outEnergy("log_energy")
	{
	TRY
	CATCH
		};
	~CSys();

	void particles_on_grid();
	void add_particle_layer(double z);
	/// 'method rain': try to release one particle at a random point across the
	/// top of the box.  Returns false if that point is already occupied.
	bool rain_one();

	void initialize(const CConfig &c);
	void solve();
	void forward(double &dt);
	void adapt(double &dt);
	void calForces();
	void computeEnergies();
	void interactions();
	inline bool interact(unsigned int i, unsigned int j)const; //force from p2 on p1
	inline bool interact(CParticle *p1,CParticle *p2)const;
	inline bool interact(CParticle *p1, BoxContainer *p2);

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
	BoxContainer walls;
	ParticleContainer particles;
	CCellList<ParticleContainer, CParticle> celllist;
	CSizeDistribution size_dist;
	

	CPlane *sp;
	double minr, maxr, maxh, maxv;
	vec G;
	const unsigned maxNParticle;
	double fluiddampping;
	double Energy, rEnergy, pEnergy, kEnergy;
 	private:
	bool do_read_radii, softwalls, spherize_on;
	double scaling;
	string out_name;
	double epsFreeze;
	vector<vec> radii;
	ifstream inputRadii;
	double maxRadii;
	vec top_v;
	/// 'method rain' state: fractional particles owed by the release rate, and
	/// whether the "box is full" message has been said.
	double rainAccum;
	/// how long the release has been blocked, to tell "nothing fits any more"
	/// from "something happens to be falling past the release height"
	double rainBlockedTime;
	/// the original gravity/dt relaxation ramp; off in config_rain
	bool relaxation;
	bool rainFullWarned;

	ofstream outEnergy;
	};

CSys::~CSys(){
	// Deliberately no TRY/CATCH here.  A destructor is implicitly noexcept
	// since C++11, so the RETHROW in the CATCH block would call std::terminate
	// rather than propagate -- and there is nothing in this body to throw.
	}

double TruncGaussRand(double r, double dr=0.0){
	if(dr<1e-6) return r;
	double x=rgen.randNorm(r, dr) ;
	if(x<r-2*dr or x> r+2*dr)return TruncGaussRand(r, dr);//trying until finding in range (r-dr, r+dr)
	return x;
	}

// Reads the simulated time out of a snapshot header, for 'method restart'.
// Returns 0 if there is no header, which is what snapshots written before the
// header existed look like -- those restarts then run for maxTime as they
// always did.
static double snapshot_time(const string &path){
	ifstream in(path.c_str());
	string line;
	while(getline(in,line)){
		if(line.empty() || line[0]!='#') continue;
		string::size_type p=line.find("t=");
		if(p==string::npos) continue;
		istringstream ss(line.substr(p+2));
		double v;
		if(ss>>v) return v;
		}
	return 0.0;
	}

CMaterial particle_material;
void CSys::initialize(const CConfig &config){
TRY
	G=config.get_param<vec>("Gravity");
	out_name=config.get_param<string>("output");
	outDt=config.get_param<double>("outDt");
	outStart=config.get_param<double>("outStart");
	outEnd=config.get_param<double>("outEnd");
	fluiddampping=config.get_param<double>("fluiddampping");
	string particleType=config.get_param<string>("particleType");
	string simul_method=config.get_param<string>("method");
	softwalls=config.get_param<bool>("softwalls");
	spherize_on=config.get_param<bool>("spherize_on");
	scaling=config.get_param<double>("scaling");
	relaxation=config.get_param<bool>("relaxation");

	double dr=config.get_param<double>("particleSizeWidth");
	DisBetaDistribution ibeta_dist(3,3,dr);

	particle_material.stiffness=paramsDouble("stiffness");
	particle_material.damping=paramsDouble("damping");
	particle_material.friction=paramsDouble("friction");
	particle_material.static_friction=paramsDouble("static_friction");
	particle_material.cohesion=paramsDouble("cohesion");
	particle_material.density=paramsDouble("density");

	// These two are read into the material but no force law ever consults
	// them: Test::contactForce uses only stiffness, damping and friction.
	// Say so rather than letting a config silently ask for physics the solver
	// does not have.  (They are kept rather than deleted because both are
	// meaningful model extensions -- see ROADMAP.md.)
	{
		struct { const char *name; const char *what; } unimpl[] = {
			{"cohesion",        "particle cohesion"},
			{"static_friction", "static friction"},
		};
		for(size_t i=0; i<sizeof(unimpl)/sizeof(unimpl[0]); ++i){
			double given=config.get_param<double>(unimpl[i].name);
			double def=config.get_param<double>(unimpl[i].name, CParamBase::Default);
			if(given!=def)
				WARNING(unimpl[i].what<<" is not implemented: '"<<unimpl[i].name
					<<"' is accepted but has no effect on the simulation");
			}
	}


	size_dist=config.get_param<CSizeDistribution>("SizeDistribution");
	

	if(simul_method=="restart"){
		const string restart_file=config.get_param<string>("input");
		// Resume the clock as well as the state.  Without this a restarted run
		// was given the whole of maxTime again rather than the remainder, which
		// is why a restart from t=0.25 with maxTime=0.3 used to run six times
		// too long.
		t=snapshot_time(restart_file);
		if(t>0)
			cerr<<"Restarting from t="<<t<<endl;
		particles.parse(restart_file);
		maxRadii=particles.maxr;
		ParticleContainer::iterator it1;
		for(it1=particles.begin(); it1!=particles.end(); ++it1){
			(*it1)->set_material(particle_material);
			}
		// Prime the accelerations.  Beeman advances the position from a_n and
		// a_{n-1}, and while the snapshot now carries the velocities it does
		// not carry the acceleration history.  Computing the forces once and
		// using them for both gives a second-order first step instead of a
		// zeroth-order one; without it a restart resumes the trajectory but
		// with a visible one-step kick.
		celllist.setup(2.0*maxRadii);
		calForces();
		for(it1=particles.begin(); it1!=particles.end(); ++it1){
			(*it1)->x(2)=*(*it1)->forces/(*it1)->get_mass();
			(*it1)->x0(2)=(*it1)->x(2);
			// ...and the rotational one.  Priming only the translation left
			// the first step with a zero angular acceleration, which is what
			// the residual difference between a restarted run and a continuous
			// one turned out to be.
			(*it1)->calAngularAccel();
			(*it1)->w0(2)=(*it1)->w(2);
			}

		// NOTE: deliberately no celllist.build() call anywhere here.  The grid
		// is not sized until celllist.setup(), so building earlier ran
		// CCellList::clear() and which() against uninitialised nx/ny/nz and
		// dx/dy/dz -- the first added particle reported
		// "(i,j,k): 2147483647 2147483647 2147483647".  calForces() rebuilds
		// the list from scratch every step anyway.
		}
	else if(simul_method=="Stillinger"){

		double eta=config.get_param<double>("eta");
		double zeta=config.get_param<double>("zeta");
		unsigned int N=config.get_param<unsigned int>("nParticle");
		double dl=1./pow((double)N,1./3.);
		std::uniform_real_distribution<double> unif(0, 1);

		cerr<<"Cell size: "<< dl <<endl;
		celllist.setup(dl);
		double dx=celllist.dx;
		double dy=celllist.dy;
		double dz=celllist.dz;
		double xx,yy,zz;
		xx=0;yy=zz=0.5*dl;
		ofstream testout("orient");
		for(unsigned int i=0; i<N; i++){
			xx+=dl;
			if(xx>walls.L(0)){xx=0.5*dx; yy+=dy;}
			if(yy>walls.L(1)){xx=0.5*dx; yy=0.5*dy; zz+=dz;}
			if(zz>walls.L(2))break;
			vec x=vec(xx-0.5*dl+0.5*unif(eng)*dl,yy-0.5*dl+0.5*unif(eng)*dl, zz-0.5*dl+0.5*unif(eng)*dl);
			double r=0.2*dl*size_dist.get();
			double a =r*pow(zeta,1./3.)/pow(eta,1./3);
			double b =r/(pow(zeta,2./3.)*pow(eta, 1/3.));
			double c =r*eta*pow(zeta/eta,1./3.);

			// (three unused rgen() draws for a hand-built quaternion used to
			// sit here, left over from a commented-out construction.  They
			// were dead but still advanced the random stream.)
			Quaternion q=randomQuaternion();
			//q=q*randomQuaternion();
			testout<<spherical(q.toBody(vec(1,0,0)))<<endl;
			//testout<<spherical(randomDirection())<<endl;
			//testout<<spherical(q.v)<<endl;
			q*=0;
			q.u=1;
			CEllipsoid E2(x, a,b,c, q);
			CParticle *p = new CParticle(E2);
			//p->w(1)(0)=5.0*(1-2*rgen());
			//p->w(1)(1)=5.0*(1-2*rgen());
			//p->w(1)(2)=5.0*(1-2*rgen());

			p->x(1)(0)=0.3*(1-2*rgen());
			p->x(1)(1)=0.3*(1-2*rgen());
			p->x(1)(2)=0.3*(1-2*rgen());
			add(p);
			}
		}
	else if(simul_method=="deposition" or simul_method=="rain"){
		if( particleType=="gen1" or particleType=="gen2"
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
			std::uniform_real_distribution<double> unif(rmin, rmax);
			
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
		else if(particleType=="general"){
			// eta = a/b and zeta = b/c (see the README).  These two shape
			// parameters must be drawn independently: this used to read
			// zeta0/zetaW and then compute zeta from eta0/etaW, so setting
			// zetaWidth did nothing and zeta silently tracked eta.
			double zeta0=config.get_param<double>("zeta");
			double zetaW=config.get_param<double>("zetaWidth");
			double eta0=config.get_param<double>("eta");
			double etaW=config.get_param<double>("etaWidth");
			double r0=config.get_param<double>("particleSize");
			for(int i=0; i<10000;i++){
				double zeta=zeta0*TruncGaussRand(1, zetaW);
				double eta =eta0*TruncGaussRand(1, etaW);
				double r=r0*ibeta_dist.rnd();
				double a =r*pow(zeta,1./3.)/pow(eta,1./3);
				double b =r/(pow(zeta,2./3.)*pow(eta, 1/3.));
				double c =r*eta*pow(zeta/eta,1./3.);
				maxRadii=max(r,max(a,max(b,c)));
				radii.push_back(vec(a, b, c));
				}
			}
		else
			ERROR(1, "Unknown particle type: "+particleType);
		}
	else{
		ERROR(1, "Unknown method: "+simul_method);
		}
	celllist.setup(2.0*maxRadii);

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
	p->set_material(particle_material);
	particles.add(p);
	particles.back()->id=particles.TotalParticlesN;
	 ++(particles.TotalParticlesN);

	celllist.add(p);

	return true;
CATCH
	}

void CSys::computeEnergies(){
	rEnergy=0; pEnergy=0; kEnergy=0;
	ParticleContainer::iterator it;
	for(it=particles.begin(); it!=particles.end(); ++it){
		rEnergy+=(*it)->rEnergy();
		pEnergy+=(*it)->pEnergy(G);
		kEnergy+=(*it)->kEnergy();
		}
	}

void CSys::calForces(){
TRY
//FROMTIME
	//FIXME make it dimensionless

	//reset forces
	ParticleContainer::iterator it1, it2, ittemp;
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		(*it1)->reset_forces(G*((*it1)->get_mass())-fluiddampping*G.abs()*(*it1)->get_mass()*(*it1)->x(1));//gravity plus damping (coef of dumping is ad hoc)
		(*it1)->reset_torques(vec(0.0));
		}

	//interactions
	for(it1=particles.begin(); it1!=particles.end(); ++it1){

		//the walls
		if(interact(*it1, &walls)){  }
		}
	celllist.build(particles);
	celllist.interact();
//TOTIME
CATCH
};



void CSys::output(ostream &out){
TRY
	// A header, so a snapshot says what wrote it and what time it is.  Readers
	// that only know the 6 and 14 records skip it, and it is what lets a
	// restarted run resume the clock instead of starting it at zero again.
	out<<"# ellipmd "<<MYMD_VERSION<<"  t="<<setprecision(14)<<t<<endl;
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
       for(it=particles.begin(); it!=particles.end(); ++it){
               //if((vec2d((*it)->x(0)(0),(*it)->x(0)(1))-vec2d(0.5, 0.5)).abs()>0.7*(1.2-(*it)->x(0)(2)))(*it)->expired=true;
               //if(abs(p1->x(0)(0)-0.5)>0.4*(1.2-p1->x(0)(2)))p1->expired=true;
	       //if(p1->x(0)(2)<config.get_param<double>("particleSize")/2.)p1->frozen=true;
               if((*it)->expired){
                       remove(particles, it);//As the side effect, "it" is set to next value
                       if(it==particles.end())break;
                       }
               }
       if(config.get_param<string>("method")=="rain" and !rainFullWarned)
               {
               // Release on a schedule rather than in layers.  rainAccum
               // carries the fraction of a particle owed from one step to the
               // next, so the rate need not be a whole number per step.
               //
               // rainRate is a ceiling, not a promise.  A particle starts from
               // rest, so one released now has only fallen 1/2 g t^2 by the
               // time the next is due: at 20/s that is 12 mm, far less than a
               // particle, so consecutive releases land inside each other and
               // jam into a column under the lid.  When a spawn point is
               // blocked the release simply waits, and the rate the simulation
               // achieves settles at whatever the fall time allows.
               rainAccum+=dt*config.get_param<double>("rainRate");
               if(rainAccum>4.0) rainAccum=4.0;   // no burst after a stall

               int released=0;
               while(rainAccum>=1.0 and particles.size()<maxNParticle
                                 and released<4){
                       bool placed=false;
                       for(int attempt=0; attempt<64 and !placed; ++attempt)
                               placed=rain_one();
                       if(!placed) break;   // blocked for now; retry next step
                       rainAccum-=1.0;
                       ++released;
                       }

               if(released>0) rainBlockedTime=0.0;
               else if(particles.size()>0){
                       // Nothing placed for a whole simulated second, with
                       // particles already in the box: the pile has reached the
                       // lid.  A fixed number of failed attempts would be
                       // wrong here -- one particle falling past the release
                       // height blocks most of the spawn area for as long as it
                       // takes to clear, which is most of a second.
                       rainBlockedTime+=dt;
                       if(rainBlockedTime>1.0){
                               rainFullWarned=true;
                               WARNING("rain: box is full; released "
                                       <<particles.size()<<" particles");
                               }
                       }
               }
       if(config.get_param<string>("method")=="deposition"
                       and maxh < walls.corner(2)+walls.L(2) ) 
               {
               // Only seed a layer whose particles will land inside the grid:
               // CCellList::which() rejects a particle whose centre is out of
               // bounds, and the old code called this unconditionally with a z
               // that could exceed the lid once the pile reached the top.  That
               // only went unnoticed because nParticle was usually exhausted
               // first; ask for more particles than the box holds and the run
               // died with "Point out of grid".  The jitter added inside
               // add_particle_layer is accounted for here.
               double z= maxh+1.02*maxRadii;
               double jitter=config.get_param<double>("particleSize")/5.0;
               if(z + jitter < walls.corner(2)+walls.L(2)) {
                       add_particle_layer(z);
                       }
               else if(particles.size()<maxNParticle) {
                       // warn once, not once per step, or a long run floods the log
                       static bool warned=false;
                       if(!warned){
                               warned=true;
                               WARNING("deposition: box is full, placed only "
                                       <<particles.size()<<" of "<<maxNParticle
                                       <<" particles requested");
                               }
                       }
               maxh=0;
               }


	// Frame generation runs on a fixed number of *steps*, not on a span of
	// simulated time.  outEvery is outDt/DT, computed once from the nominal
	// time step, so the rate at which frames are written does not change when
	// dt changes -- which it does if the relaxation ramp below is left on.
	// That is what "the image rate is independent of the physical time axis"
	// means here: how often a picture is taken is a property of the run, not of
	// how far the clock has moved.
	//
	// The consequence is worth stating: if dt does change during a run, frames
	// are no longer evenly spaced in simulated time.  Set 'relaxation 0' and dt
	// stays put, and the two are the same thing.
	static long stepCount=0, outN=0;
	static const long outEvery=(long)(outDt/DT);
	static ofstream out;

	//this is for a messure of performance
	static double starttime=clock();
	if(outEvery>0 and stepCount%(outEvery/10+1)==0)
		cout<<(clock()-starttime)/CLOCKS_PER_SEC<< "   "<<t<<endl;

	if(outEvery>0 and stepCount%outEvery==0 and t>=outStart and t<=outEnd){
			stringstream outstream;
			outstream<<out_name<<setw(5)<<setfill('0')<<outN;
			output(outstream.str());
			// The energies are accumulated at the *end* of forward(), so on the
			// very first snapshot they have never been computed and the old
			// code wrote uninitialised doubles (~1e-314) as the t=0 line of
			// log_energy.  Recompute them from the current state instead.
			if(outN==0) computeEnergies();
			outN++;
			Energy=rEnergy+kEnergy+pEnergy;
			outEnergy<<setprecision(14)<<t<<"  "<<Energy<<"  "<<kEnergy<<"  "<<pEnergy<<"  "<<rEnergy <<endl;
			rEnergy=0; pEnergy=0; kEnergy=0; Energy=0;
			// The original relaxation protocol: once t>2, weaken gravity by
			// 10% and lengthen the step by 8%, once per output, until gravity is
			// down to 1.  It is a compression/relaxation device, and because it
			// steps once per frame it makes the physics depend on outDt.  Set
			// 'relaxation 0' for a run with constant gravity.
			if(relaxation and t>2 and G.abs()>1){
					G*=0.9;
					dt*=1.08;	
					cerr<<"t: "<<t<<" G: "<<G<<" dt: "<<dt<<endl;
					}
			if(spherize_on)for(it=particles.begin(); it!=particles.end(); ++it){
				(*it)->shape->spherize();
				}
			static bool done=false;
			if(scaling>1.0000001){
				if(!done && particles.totalVolume()>0.5 ){done=true;scaling=1.01;}
				for(it=particles.begin(); it!=particles.end(); ++it){
				(*it)->scale(scaling);
				//(*it)->x(1)=0;
				//(*it)->w(1)=0;
				}
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
	++stepCount;
	double vtemp;
	maxv=0;
	for(it=particles.begin(); it!=particles.end(); ++it){
	//	if(!(*it)->frozen) 
		(*it)->calVel(dt);
		vtemp=(*it)->x(1).abs()+(*it)->w(1).abs()*(*it)->shape->radius;
		if(vtemp>maxv) 
			maxv=vtemp;

		rEnergy+=(*it)->rEnergy();
		pEnergy+=(*it)->pEnergy(G);
		kEnergy+=(*it)->kEnergy();
		//cout<< it->x <<"  "<<it->size<< " cir"<<endl;
		}
	
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
			if(stop)
				cerr<<"Reached maxTime at t="<<t<<endl;
			else
				cerr<<"Relaxation criterion reached at time="<<t<<": KE= "<<kEnergy<< " < 1e-8"<<endl;
			break;
			}
		//adapt(dt);
		}
	}catch(CException &e){
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
	// see the note in interaction_force.h: this used to index CParticle::vlist
	// by the walls pointer
	static ShapeContact overlaps;
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


bool CSys::rain_one(){
TRY
	// Same shape sampling as deposition, but one particle at a time and at a
	// random point rather than on a grid.
	const vec abc=radii.at(rgen.rand(radii.size()));
	const double rr=tmax(abc(0), tmax(abc(1), abc(2)));

	// The centre has to stay inside the grid for CCellList::which(), and the
	// particle should appear right at the top, so the highest valid centre is
	// one radius below the lid.
	const double lo0=walls.corner(0)+rr, hi0=walls.corner(0)+walls.L(0)-rr;
	const double lo1=walls.corner(1)+rr, hi1=walls.corner(1)+walls.L(1)-rr;
	const double z  =walls.corner(2)+walls.L(2)-rr;
	ERROR(hi0<lo0 or hi1<lo1 or z<walls.corner(2)+rr,
	      "rain: a particle of radius "+stringify(rr)
	      +" does not fit in a box of "+stringify(walls.L(0))+" x "
	      +stringify(walls.L(1))+" x "+stringify(walls.L(2)));

	vec x(lo0+(hi0-lo0)*rgen(), lo1+(hi1-lo1)*rgen(), z);

	Quaternion q=randomQuaternion();

	// Is the point free?  Broad phase on the circumscribed spheres, then the
	// exact ellipsoid test.
	//
	// The broad phase alone is not enough here, and that is the whole reason
	// for the second stage: the circumscribed radius of these spheroids is
	// about 1.4x their equivalent-sphere radius, so a single particle falling
	// past the top of a 1x1 box blocks most of it and the rain stalls after a
	// handful of particles.
	CEllipsoid E;
	bool have_E=false;
	ShapeContact ovs;
	ParticleContainer::iterator it1;
	for(it1=particles.begin(); it1!=particles.end(); ++it1){
		const double need=rr+(*it1)->shape->radius;
		if(((*it1)->x(0)-x).abs2()>=need*need) continue;
		if(!have_E){ E=CEllipsoid(x, abc(0), abc(1), abc(2), q); have_E=true; }
		if(doOverlap(ovs, E, *static_cast<CEllipsoid*>((*it1)->shape)))
			return false;
		}
	if(!have_E) E=CEllipsoid(x, abc(0), abc(1), abc(2), q);

	CParticle *p=new CParticle(E);
	// Released, not thrown: CFreedom::init() already zeroed the linear
	// velocity, so the particle starts falling from rest.  The spin is what
	// stops the packing being a stack of identically aligned particles.
	const double spin=config.get_param<double>("rainSpin");
	p->w(1)(0)=spin*(1-2*rgen());
	p->w(1)(1)=spin*(1-2*rgen());
	p->w(1)(2)=spin*(1-2*rgen());
	add(p);
	return true;
CATCH
	return false;
	}

void CSys::add_particle_layer(double z){ 
	double size=config.get_param<double>("particleSize");
	vec x(0.0, 0.0, .0);

	
	//double asphericity=config.get_param<double>("asphericity");
	//double asphericityWidth=config.get_param<double>("asphericityWidth");

	double xtemp=0, ytemp=0;

	for(int i=0;i<celllist.nx;i++){
		for(int j=0;j<celllist.ny;j++){
		if(particles.size()>=maxNParticle)break;
		double a, b, c;
		//spheroid(a, b, c, asphericity, asphericityWidth);

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
		// Safety net: the seeded centre must be inside the grid or
		// CCellList::which() aborts the run.  See the gate in CSys::forward.
		for(int k=0; k<3; ++k){
			double lo=walls.corner(k)+1e-9;
			double hi=walls.corner(k)+walls.L(k)-1e-9;
			if(x(k)<lo)x(k)=lo;
			if(x(k)>hi)x(k)=hi;
			}
		double alpha=rgen()*M_PI;
                double beta=rgen()*M_PI;
               	double phi=rgen()*M_PI;
               	Quaternion q=Quaternion(cos(alpha),sin(alpha),0,0)*Quaternion(cos(beta),0,0,sin(beta))*Quaternion(cos(phi),0,sin(phi),0 );
		//CParticle *p = new CParticle(CSphere(x,size*(1-0.0*rgen())));
		//CEllipsoid E(x, 1-0.0*rgen(), 1-0.0*rgen(),1-0.0*rgen(), size*(1+0.0*rgen()));
		//CParticle *p = new CParticle(E);
		//CSphere E1(x,maxRadii);
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
