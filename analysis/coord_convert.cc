#include<iostream>
using namespace std;
#include"../exception.h"
#include"packing.h"

 
int main(int n_params, char **params){
	try {
	ERROR(n_params!=3, "Usage: convert [-r/-p] input-file");
	ERROR(params[1][0]!='-', "Usage: convert [-r/-p] input-file");
	
	ifstream inputFile(params[2]);
	ERROR(!inputFile.good(), "Unable to open input file");
	
	CPacking<CEllipsoid> packing;
	packing.parse(inputFile);
	if(params[1][1]=='r')
		packing.printRaster3D(cout);
	else if(params[1][1]=='p')
		packing.printEuler(cout);
	else
		ERROR(params[1][0]!='-' or (params[1][1]!='r' and params[1][1]!='p'), "Usage: convert [-r/-p] input-file");

	return 0;
	} catch(CException e)
	{
	e.Report();
	return 1;
	}
}
