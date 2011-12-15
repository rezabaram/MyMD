#include<iostream>
#include"../include/define_params.h"
#include"../include/particle.h"
#include"../include/packing.h"
using namespace std;

long RNGSeed;
extern MTRand rgen;
CPacking<CParticle> packing;


void Initialize(int n_params, char **params){
	rgen.seed(RNGSeed);
//	cerr<< "RNG Seed: "<< RNGSeed <<endl;

	ERROR(n_params!=2, "Usage: convert2raster input-file");
	
	ifstream inputFile(params[1]);
	ERROR(!inputFile.good(), "Unable to open input file"+(string)(params[1]));
	
	define_parameters();
	packing.parse(inputFile, false);
	}

class SMap
	{
	public:
	SMap(unsigned long int n){
		}
	
	//void step
 	private:
	Vec<3,unsigned long int> index;
	};

void Run(){
	
	packing.BuildGrid();
	vec x(10,13,17);
	packing.test(x);

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
