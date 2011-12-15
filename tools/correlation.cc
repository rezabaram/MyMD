

#include<iostream>
#include<fstream>
#include<string>
#include"../include/define_params.h"
#include"../include/particle.h"
#include"../include/packing.h"
#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
using namespace std;

#include <hdf5.h>
#define MAXSTRING 160
#define RANK  3
#define DIM1 128 
#define DIM2  128
#define DIM3  128

typedef float dtype;

class CLattice{
	public:
	CLattice(long int _nx, long int _ny, long int _nz, string input): nx(_nx), ny(_ny), nz(_nz){
		//allocate();
		ifstream inputFile(input.c_str());
		ERROR(!inputFile.good(), "Unable to open input file: "+input);
		define_parameters();
		packing.parse(inputFile, true);
		cerr<< packing.size()  <<endl;
		packing.BuildGrid();
		}
	~CLattice(){
		//deallocate();
		}
	
	void deallocate(){
		
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++)
				delete [] data[i][j];
			}

		for(int i=0; i<nx; i++){
			delete [] data[i];
			}
		delete [] data;
		}
	void allocate(){
		data=new dtype **[nx];
		for(int i=0; i<nx; i++){
			data[i]=new dtype *[ny];
			for(int j=0; j<ny; j++)
				data[i][j]=new dtype [nz];
			}

		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++)
				data[i][j]=new dtype [nz];
			}

		}
	void correlation(double r){
		static double phi=1-packing.packFraction(vec(0,0,0),vec(1,1,1), 1000000 );
		double dx=1.0/nx, dy=1.0/ny, dz=1.0/nz;
		double n=0;
		int  i, j,k;
		vec x1, x2;
		double kapa1, kapa2, nom=0, denom=0;
		 for ( i = 0; i < nx; i++ ) 
		 for ( j = 0; j < ny; j++ )
		 for ( k = 0; k < nz; k++ ){
			//if(k*dz<0.2 or k*dz >0.9)continue;
			n++;
		        //data[i][j][k] = packing.is_in_void(vec((i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz)); 
		        x1=vec((i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz); 
			double temp=(j+0.5)*dx+r;
			if(temp>1)temp--;//periodic condition
		        //x2=vec(temp, (j+0.5)*dy, (k+0.5)*dz); 
		        x2=vec((i+0.5)*dx, temp, (k+0.5)*dz); 
			kapa1=kapa2=0;
		        if( packing.is_in_void(x1)) kapa1=1;
		        if( packing.is_in_void(x2)) kapa2=1;
			nom+=(kapa1-phi)*(kapa2-phi);
			denom+=(kapa1-phi)*(kapa1-phi);
			}
		cout<< r*2250. <<"   "<< nom/denom <<endl;
	}
	

CPacking<CParticle> packing;
private:
long int nx, ny, nz;
dtype ***data;
};




void Initialize(int n_params, char **params){
	//lattice.out_hdf("packing.h5");
//	cerr<< "RNG Seed: "<< RNGSeed <<endl;

	ERROR(n_params!=3, "Usage: coord2xdr <input-file> <exponent>");
	
	string input=(string)(params[1]);

	int e=atoi(params[2]);
	cerr<< "Resolution: "<< e <<"X"<<e<<"X"<<e<<endl;
	CLattice lattice(e, e, e,input);
	lattice.allocate();

	for(double r=0; r<0.2; r+=0.005){
		lattice.correlation(r);
		}
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
