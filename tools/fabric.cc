#include<iostream>
#include"../include/packing.h"
#include"../include/particle.h"
#include"../include/ellips_contact.h"
using namespace std;

CPacking<CParticle> packing;

void Initialize(int n_params, char **params)
{

	if(n_params!=2){
		cerr<< "Usage: "+(string)params[0]+"  input-file"<<endl;
		exit(1);
		}
	
	ifstream inputFile(params[1]);
	ERROR(!inputFile.good(), "Unable to open input file");
	
	packing.parse(inputFile, true);

}

void Run()
{
	packing.BuildContactNetwork();
	packing.contacts.print_eigen(cout, vec(0.0,0.0,0.20), vec(1.0,1.0,.95));
	//packing.contacts.print_branch_vectors(cout);
	//packing.contacts.print_branch_vectors(cout, vec(0.0,0.0,0.15), vec(1,1,.95));
	//packing.print_particle_axes(cout, vec(0.0,0.0,0.15), vec(1,1,.95));;
	//cout<<packing.packFraction(vec(0.0,0.0,0.15), vec(1,1,0.95))<<"\t";
	//cout<<packing.avg_contact_number()<<endl;
	
	//ofstream out("network");
	//packing.output_contact_network(out);
	//cerr<< <<endl;
	
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
