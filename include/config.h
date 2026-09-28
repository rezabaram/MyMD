// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


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

	/// Set one parameter from its textual form, as it would appear in a config
	/// file.  Used by the command line's --set KEY=VALUE, which is the point:
	/// a parameter sweep should not have to write a config file per point.
	void set(string name, string value){
		validate_param(name);
		istringstream ss(value);
		params.find(name)->second->parse(ss);
		}

	void parse(string fname);
	// CBaseConfig declares print(ostream&, out_type)const.  A differently
	// shaped print() here would hide it rather than overload it, so
	// config.print(cout) would not compile on a CConfig.  Re-expose it.
	using CBaseConfig::print;
	void print(string fname)const;

	void parse(istream &inputFile) {
		string line;
		string vname;

		//Parse the line
		while(getline(inputFile,line))
		{

		line = line.substr( 0, line.find(comm) );

		//Insert the line string into a stream
		stringstream ss(line);

		//Read up to first whitespace
		ss >> vname;

		
		//Skip to next line if this one starts with a # (i.e. a comment)
		if(vname.find("#",0)==0) continue;

		// A blank line, or one holding nothing but a comment, leaves vname
		// empty; without this it was reported as an unknown parameter.
		if(vname.empty()) continue;

		if(!isValidParam(vname)){
			cerr<< "Warning: "<<vname<<" is not a valid parameter or keyword" <<endl;
			continue;
			}

		params[vname]->parse(ss);
		}
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

	       add_param<double>("stiffness", 5.0e+02); 
	       add_param<double>("damping", 5); 
	       // 'method rain': particles are released one at a time from the top
	       // of the box at zero linear velocity, with a small random spin.
	       // rainRate is a ceiling in particles per unit of simulated time; the
	       // rate actually achieved is limited by how fast a particle released
	       // from rest clears the release height, about sqrt(4*particleSize/g).
	       add_param<double>("rainRate", 4.0);
	       add_param<double>("rainSpin", 2.0);
	       add_param<double>("fluiddampping", 0.05); 
	       add_param<double>("friction", 0.2); 
	       add_param<double>("static_friction", 0); 
	       add_param<double>("cohesion", 0); 
	       add_param<double>("density", 1.0); 
	       add_param<double>("particleSize", 1); 
	       add_param<double>("rmin", 0.05); 
	       add_param<double>("rmax", 0.05); 

	       add_param<double>("particleSizeWidth", 0); 
	       add_param<double>("timeStep", 0.00001); 
	       add_param<double>("maxTime", 10.0); 
	       add_param<unsigned int>("nParticle", 5); 
	       add_param<string>("particleType", "general"); 
	       add_param<string>("method", "deposition"); 
	       add_param<double>("zeta", 1.0); 
	       add_param<double>("zetaWidth", 0.0); 
	       add_param<double>("eta", 1.0); 
	       add_param<double>("etaWidth", 0); 
	       add_param<double>("asphericity", -0.5); 
	       add_param<double>("asphericityWidth", 0.1); 
	       add_param<double>("scaling", 1.0); 

	       add_param<string>("boundary", "solid"); 

	       // The original relaxation protocol: after t>2, gravity is weakened by
	       // 10% and the time step lengthened by 8% once per output until
	       // gravity reaches 1.  Kept as the default so existing configurations
	       // behave as they always did; turn it off for a plain constant-gravity
	       // run, which is what the rain animation uses.
	       add_param<bool>("relaxation", true); 
	       add_param<bool>("softwalls", false); 
	       add_param<bool>("spherize_on", false); 
	       add_param<string>("input", "input.dat"); 
	       add_param<string>("radii", "radii.dat"); 
	}
};

void CConfig::print(string outname)const{
	ofstream outputFile(outname.c_str());
	if(!outputFile.good() ){
		cerr << "WARNING: Unable to open input file: " << outname << endl;
		return;
		}
	CBaseConfig::print(outputFile);
}
void CConfig::parse(string infilename) {

	ifstream inputFile(infilename.c_str());

	if(!inputFile.good())
	{
	// This used to warn and carry on with the compiled-in defaults, which are
	// not a working configuration -- particleSize=1 inside a 1x1x2 box puts
	// every particle outside the grid, and the run dies later with a confusing
	// "Point out of grid".  Better to say so here.
	ERROR(1, "cannot open config file '"+infilename+"'\n"
	         "\tPass one as the second argument:  ./ellipmd <seed> <config-file>\n"
	         "\tOr use one of the examples:        make run CONFIG=config_quick");
	}
	parse(inputFile);
	inputFile.close();
}
#endif /* DEFINE_PARAMS_H */
