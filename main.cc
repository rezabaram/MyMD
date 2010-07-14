#include<iostream>
#include<stdlib.h>
#include"mdsys.h"
using namespace std;


long RNGSeed;
void Initialize(){
	RNGSeed=0;
	define_parameters();
	config.parse("config");
	}

void Run(){

	CSys sys(config.get_param<size_t>("nParticle"));
	sys.initialize(config);
	sys.solve();
}

int main(){
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
