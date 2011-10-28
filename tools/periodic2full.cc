#include<iostream>
#include<fstream>
#include<string>
#include"../include/define_params.h"
#include"../include/particle.h"
#include"../include/packing.h"
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#define MAXSTRING 160
#define RANK  3
#define DIM1 128 
#define DIM2  128
#define DIM3  128


void Initialize(int n_params, char **params){
	//lattice.out_hdf("packing.h5");
//	cerr<< "RNG Seed: "<< RNGSeed <<endl;

	ERROR(n_params!=2, "Usage: periodic2full <input-file> ");
	
	string input=(string)(params[1]);

	CPacking<CParticle> packing;
	packing.parse(input, true);
	string output=input+"full";
	ofstream out(output.c_str());
	packing.print(out);
	}

void Run(){
	//cout<< packing.packFraction(vec(0.1,0.1,0.1),vec(.9,.9,.9), 100000 ) <<endl;
	//cerr<< packing.totalVolume()<<endl;
}

int main(int n_params, char **params){
	
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
