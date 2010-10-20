#include<iostream>
using namespace std;
#include"../exception.h"
#include"packing.h"

 
int main(int n_params, char **params){
	try {
	ERROR(n_params!=2, "Usage: convert2raster input-file");
	
	ifstream inputFile(params[1]);
	ERROR(!inputFile.good(), "Unable to open input file");
	
	CPacking<GeomObjectBase> packing;
	packing.parse(inputFile);
	packing.printRaster3D(cout);

	return 0;
	} catch(CException e)
	{
	e.Report();
	return 1;
	}
}
