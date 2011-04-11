#include<iostream>
#include<stdlib.h>
#include"include/main.h"
#include"include/mdsys.h"
using namespace std;


long RNGSeed;
extern MTRand rgen;

void Initialize(){
	rgen.seed(RNGSeed);
	define_parameters();
	//config.parse("config");
	config.parse(cin);
	}

void Run(){

	CSys sys(config.get_param<unsigned int>("nParticle"));
	sys.initialize(config);
	sys.solve();
	}

int main(int pi, char **params){
	if(pi==1)
		RNGSeed=0;
	else
		RNGSeed=313*atoi(params[1])+1;

	cerr<<"RNG Seed: "<<RNGSeed<<endl;
	
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
