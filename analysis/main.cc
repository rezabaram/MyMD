#include<iostream>
#include"packing.h"
using namespace std;

long RNGSeed;
extern MTRand rgen;

void Initialize(){
	rgen.seed(RNGSeed);
	cerr<< "RNG Seed: "<< RNGSeed <<endl;
	//define_parameters();
	//config.parse("config");
	}

void Run(){

	CPacking packing;
	packing.parse("out00005");
	cerr<< packing.packFraction(vec(0.3,0.3,0.3),vec(.7,.7,.7), 10000 ) <<endl;
	cerr<< packing.totalVolume()<<endl;
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
