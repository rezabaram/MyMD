#ifndef IBETA_DIST_H
#define IBETA_DIST_H 
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<math.h>
#include"MersenneTwister.h"
#include"exception.h"
//this class is to generate a distribution for particles radii
//whose volumes
using namespace std;
class DisBetaDistribution
	{
	template<typename T>
	double d_ibeta(T a, T b, double x, double min=0, double max=1)
		{
		//return (x<=min || x>=max)?0:pow(x-min,a-1)*pow(max-x,b-1);
		return (x<=min || x>=max)?0:pow(x-min,a-1)*pow(max-x,b-1)/(x*x*x);
		//return (x<=min || x>=max)?0:ibeta_derivative(a,b,(x-min)/(max-min))/(x*x*x);
		}
	public:
	template<typename T>
	DisBetaDistribution(T a, T b, double width=0.5, int _nbins=100, int npool=10000){
		ERROR(width>1, "The width of the distribution is larger than its median.");
		if(width<=1e-10){
			dx=0;
			min=1;
			max=1;
			return;
			}
		min=1-width;
		max=1+width;
		nbins=_nbins;
		bins=new double[nbins];
		dx=(max-min)/double(nbins);
		double sum=0;
		for(int i=0; i<=nbins; i++){
			bins[i]=d_ibeta(a,b,min+dx*(double)i, min, max);
			sum+=bins[i]*dx;
			}
		//normalize
		for(int i=0; i<=nbins; i++)
			bins[i]/=sum;

		for(int i=0; i<=nbins; i++){
			int n=npool*bins[i];
			for(int j=0; j<n; j++){
				pool.push_back(i);
				}
			}
		
		}

		double rnd(){
			if(dx<1e-10)return min;
			return min+dx*pool.at((int)(rgen.rand()*pool.size()));
			}
			
		void print(){
			for(int i=0; i<=nbins; i++){
				cout<< (min+dx*i) <<"\t"<<bins[i] <<endl;
				}
			}
	
 	private:
	double *bins;
	vector<int> pool;
	double dx;
	double min, max;
	int nbins;
	};
#endif /* IBETA_DIST_H */
