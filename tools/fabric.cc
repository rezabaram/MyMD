#include<iostream>
#include"../include/packing.h"
#include"../include/particle.h"
#include"../include/ellips_contact.h"
#include"../include/define_params.h"
using namespace std;

CPacking<CParticle> packing;

void Initialize(int n_params, char **params)
{
	define_parameters();
	config.parse("config");

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
	ofstream out("network");
	packing.output_contact_network(out);
	cerr<< packing.contacts.cal_fabric_tensor() <<endl;
	
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