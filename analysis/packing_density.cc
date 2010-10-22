#include<iostream>
#include"packing.h"
using namespace std;

long RNGSeed;
extern MTRand rgen;
CPacking<CEllipsoid> packing;

void Initialize(int n_params, char **params){
	rgen.seed(RNGSeed);
	cerr<< "RNG Seed: "<< RNGSeed <<endl;

	ERROR(n_params!=2, "Usage: convert2raster input-file");
	
	ifstream inputFile(params[1]);
	ERROR(!inputFile.good(), "Unable to open input file");
	
	packing.parse(inputFile);
	}

void Run(){

	cerr<< packing.size() <<endl;
	cerr<< packing.packFraction(vec(0.1,0.1,0.1),vec(.9,.9,.9), 10000 ) <<endl;
	cerr<< packing.totalVolume()<<endl;
}

int main(int n_params, char **params){
	if(n_params==1)
		RNGSeed=0;
	else
		RNGSeed=313*atoi(params[1])+1;
	
	try {
	Initialize(n_params, params);
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
