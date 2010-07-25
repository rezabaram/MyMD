#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;


long RNGSeed;
MTRand rgen(RNGSeed);

void Initialize(){
	rgen.seed(RNGSeed);
	define_parameters();
	config.parse("config");
	}

void Run(){

	CSys sys(config.get_param<size_t>("nParticle"));
	sys.initialize(config);
	sys.solve();
}

int main(int pi, char **params){
	if(pi==1)
		RNGSeed=0;
	else
		RNGSeed=313*atoi(params[1])+1;
	
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
