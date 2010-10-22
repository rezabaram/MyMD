#include<iostream>
#include"packing.h"
using namespace std;
#include"../ellips_contact.h"

CPacking<CEllipsoid> packing;

void Initialize(int n_params, char **params)
{
	if(n_params!=2){
		cerr<< "Usage: "+(string)params[0]+"  input-file"<<endl;
		exit(1);
		}
	
	ifstream inputFile(params[1]);
	ERROR(!inputFile.good(), "Unable to open input file");
	
	packing.parse(inputFile);
}

void Run()
{
	packing.BuildContactNetwork();
	//packing.output_contact_network(cout);
}

int main(int n_params, char **params){
	
	try 
	{
		Initialize(n_params, params);
		Run();
		//Shutdown();
		return 0;
	} 
	catch(CException e)
	{
	e.Report();
	return 1;
	}
return 0;
}
