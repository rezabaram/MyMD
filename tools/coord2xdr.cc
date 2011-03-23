

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
class CXdr{
	public:
	CXdr(string filename){
		file.open(filename.c_str());
		char   command[MAXSTRING];
		sprintf(command,"touch %s ", filename.c_str());
		system(command);
		FILE   *f2=fopen(filename.c_str(), "w");
		if(NULL==f2){
			fprintf(stderr,"Unable to open output file, aborting... ");
			exit(1);
			}
	  	xdrstdio_create(&xdrs,f2,XDR_ENCODE);
		};
		template<typename T>
		void add(T &x){
			xdr_float(&xdrs,&x);//data[k][j][i]);
			}
	private:
	XDR xdrs;
	ofstream file;
	};

class CLattice{
	public:
	CLattice(int _nx, int _ny, int _nz, string input):m_xdr(new CXdr(input+".xdr")), nx(_nx), ny(_ny), nz(_nz){
		allocate();
		ifstream inputFile(input.c_str());
		ERROR(!inputFile.good(), "Unable to open input file: "+input);
		define_parameters();
		packing.parse(inputFile);
		packing.BuildGrid();
		}
	~CLattice(){
		deallocate();
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

		int  i, j,k;
		 for ( i = 0; i < nx; i++ ) 
		 for ( j = 0; j < ny; j++ )
		 for ( k = 0; k <nz; k++ )
		     data[i][j][k] = i+j+k;
		}

	void out_xdr(){
	 double dx=1.0/nx, dy=1.0/ny, dz=1.0/nz;
	  /* Create xdrfile, stop if it already exists/overwrite? */
	    for(int i=0;i<nx;++i) {
	    for(int j=0;j<ny;++j) {
	    for(int k=0;k<nz;++k) {
		/* Check output of xdr file*/
		float x=1;
		if(packing.is_in_void(vec(i*dx, j*dy, k*dz)))x=0;
		m_xdr->add(x);
	      }
	    }
	  }
	}
	void out_hdf(char *FILENAME)
	{
	     hid_t   file, dataset;
	     hid_t  status, fid, cparms;
	     hsize_t fdim[] = {DIM1, DIM2, DIM3};
	     hsize_t chunk_dims[] = {32,32,32};
	     cparms=H5Pcreate(H5P_DATASET_CREATE);
	     status=H5Pset_chunk(cparms, RANK, chunk_dims);

	     herr_t ret;

		int  i, j,k;
		 for ( i = 0; i < DIM1; i++ ) 
		 for ( j = 0; j < DIM2; j++ )
		 for ( k = 0; k <DIM3; k++ ){
		     data[i][j][k] = 0;
			}


	     file = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	     fid = H5Screate_simple (RANK, fdim, NULL);

	     dataset = H5Dcreate (file, "packing", H5T_NATIVE_INT, fid, cparms);

	     ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	     ret = H5Dclose (dataset);
	     ret = H5Sclose (fid);
	}

CXdr *m_xdr;
CPacking<CParticle> packing;
private:
int nx, ny, nz;
dtype ***data;
};




void Initialize(int n_params, char **params){
	//lattice.out_hdf("packing.h5");
//	cerr<< "RNG Seed: "<< RNGSeed <<endl;

	ERROR(n_params!=3, "Usage: coord2xdr <input-file> <exponent>");
	
	string input=(string)(params[1]);

	int n=atoi(params[2]);
	int e=1;
	e=e<<n;
	CLattice lattice(e, e, e,input);
	lattice.out_xdr();
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