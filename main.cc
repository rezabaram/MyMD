#include"include/main.h"
using namespace std;


long RNGSeed;
extern MTRand rgen;

string config_file="config";

void Initialize(){
	rgen.seed(RNGSeed);
	define_parameters();
	config.parse(config_file);
	//config.parse(cin);
	}

void Run(){

	CSys sys(config.get_param<unsigned int>("nParticle"));
	sys.initialize(config);
	sys.solve();
	}

int main(int pi, char **params){
	if(pi==1)
		RNGSeed=0;
	if(pi>=2)
		RNGSeed=313*atoi(params[1])+1;
	if(pi>=3)
		config_file=(string)params[2];

	cerr<<"RNG Seed: "<<RNGSeed<<endl;
	cerr<<"Config file: "<<config_file<<endl;
	
	try {
	Initialize();
	Run();
	//Shutdown();
	return 0;
	} catch(CException e)
	{
	e.Report();
	return 1;
	}
return 0;
}
