#ifndef DEFINE_PARAMS_H
#define DEFINE_PARAMS_H 
#include"baseconfig.h"
#include"vec.h"
#include"size_dist.h"
#include<string>
using namespace std;

class CConfig : public CBaseConfig{
        public:
        CConfig(string filename){
                define_parameters();
                parse(filename);
                }

        CConfig(){
                define_parameters();
                }

	void define_parameters()
	{
	       add_param<CSizeDistribution>("SizeDistribution", CMonoDist(1.0));

	       add_param<vec>("Gravity", vec(0.0, 0.0, -10.0));
	       add_param<vec>("boxcorner", vec(0.0, 0.0, 0.0));
	       add_param<vec>("boxsize", vec(1.0, 1.0, 2.0));

		//controlling the output
	       add_param<double>("outStart", 0.00);
	       add_param<double>("outEnd", 1000.00);
	       add_param<double>("outDt", 0.02);
	       add_param<string>("output", "out");
	       add_param<string>("outDensity", "density");

	       add_param<double>("stiffness", 5.0e+02); 
	       add_param<double>("damping", 5); 
	       add_param<double>("fluiddampping", 0.05); 
	       add_param<double>("friction", 0.2); 
	       add_param<double>("static_friction", 0); 
	       add_param<double>("friction_threshold", 0); 
	       add_param<double>("cohesion", 0); 
	       add_param<double>("density", 1.0); 
	       add_param<double>("particleSize", 0.05); 
	       add_param<double>("rmin", 0.05); 
	       add_param<double>("rmax", 0.05); 

	       add_param<double>("particleSizeWidth", 0); 
	       add_param<double>("timeStep", 0.00001); 
	       add_param<double>("maxTime", 4.0); 
	       add_param<double>("verletfactor", 0.1); 
	       add_param<unsigned int>("nParticle", 5); 
	       add_param<string>("particleType", "general"); 
	       add_param<string>("method", "generate"); 
	       add_param<double>("e", 0.5); 
	       add_param<double>("eta", 1.0); 
	       add_param<double>("etaWidth", 0.0); 
	       add_param<double>("xi", 1.0); 
	       add_param<double>("xiWidth", 0); 
	       add_param<double>("asphericity", -0.5); 
	       add_param<double>("asphericityWidth", 0.1); 
	       add_param<double>("scaling", 1.0); 

	       add_param<string>("boundary", "wall"); 

	       add_param<bool>("read_radii", false); 
	       add_param<bool>("softwalls", false); 
	       add_param<bool>("spherize_on", false); 
	       add_param<string>("input", "input.dat"); 
	       add_param<string>("radii", "radii.dat"); 
	}
};

#endif /* DEFINE_PARAMS_H */
