#include<iostream>
#include"../include/particle.h"
#include"../include/packing.h"
#include"../include/slice.h"
using namespace std;

long RNGSeed;
extern MTRand rgen;
CPacking<CParticle> packing;

string	inputFileName;

void Initialize(int n_params, char **params){
	rgen.seed(RNGSeed);
//	cerr<< "RNG Seed: "<< RNGSeed <<endl;

	ERROR(n_params!=2, "Usage: convert2raster input-file");
	
	inputFileName=(params[1]);
	ifstream inputFile(params[1]);
	ERROR(!inputFile.good(), "Unable to open input file"+(string)(params[1]));
	
	packing.parse(inputFile, true);
	}

void Run(){
	
	vec x1 (0, 0, 0.10000);
	vec x2 (1, 1, 0.9);
	double overalscale=1;
	double scale=1;
	while(1){
		cerr<< "Give scale: ";
		cin>>scale;
		if(scale<=0)break;
		overalscale*=scale;
		cout<< packing.packFraction(x1,x2, 100000, scale ) <<endl;
		}
	cerr<<"Final scale factor: "<<overalscale  <<endl;
	ofstream output((inputFileName+"_scaled").c_str());
	packing.print(output);
	CSlice slice(vec2d(0.0),vec2d(1.0),1./256);
	double density;
	if(0)for(double z=0.0; z<1.01; z+=0.02){
		density=packing.get_slice(slice, z, 1000000);
		slice.reset();
		cout<< z<<"  "<<density <<endl;
	}
	density=packing.get_slice_uniform(slice, 0.5, 256);
	slice.export_ppm("out.ppm");
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
