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

	void parse(string fname);
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
	       add_param<double>("maxTime", 10.0); 
	       add_param<double>("verletfactor", 0.1); 
	       add_param<unsigned int>("nParticle", 5); 
	       add_param<string>("particleType", "general"); 
	       add_param<string>("method", "deposition"); 
	       add_param<double>("e", 0.5); 
	       add_param<double>("eta", 1.0); 
	       add_param<double>("etaWidth", 0.0); 
	       add_param<double>("xi", 1.0); 
	       add_param<double>("xiWidth", 0); 
	       add_param<double>("asphericity", -0.5); 
	       add_param<double>("asphericityWidth", 0.1); 
	       add_param<double>("scaling", 1.0); 

	       add_param<string>("boundary", "solid"); 

	       add_param<bool>("read_radii", false); 
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
	cerr << "WARNING: Unable to open input file: " << infilename << endl;
	return;
	}
	parse(inputFile);
	inputFile.close();
}
#endif /* DEFINE_PARAMS_H */
