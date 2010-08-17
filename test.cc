#include<iostream>
using namespace std;
#include"/home/reza/workstation/mysrc/CStat.h"
#include"polynom.h"
//#include"eigen.h"
#include"vec.h"
#include"MersenneTwister.h"
#include"common.h"
void polytest(){

vector<double> coefs;
coefs.push_back(1);
coefs.push_back(-2);
coefs.push_back(3);
coefs.push_back(4);
coefs.push_back(1);

//CPolynom<4,double> poly(coefs);
CCubic cube(2*3, -2*2, 2*5, 2*1);
//cube.plot(cout, -3, 3, 0.01);
//cube.print_roots(cerr);

bool b=true;
if(0)for(int i=0; i<100000; ++i){
	//CQuartic poly(1, 10*drand48(), drand48(), drand48(), 10*drand48());
	CCubic poly(1, 10*drand48(), drand48(), drand48());
	poly.solve();
	poly.print(cout);
	b= b and poly.check_roots();
	if(!b)break;
	}
cerr<< b <<endl;

CCubic poly(1 , 0.485717 , 0.0786405 , 0.690939);
poly.solve();
//poly.check_roots();
poly.print_roots(cerr);
poly.plot(cout, -3, 3, 0.01);
}

void gaussElimTest(){
	drand48();
	const size_t N=15, M=15;
	Matrix A(N,M);
	Matrix v(M,1);
	
	int k=0;
	for(int i=0; i<N; ++i){
		v(i,0)=+(int)(5*(1-2*drand48()));
	for(int j=0; j<M; ++j){
		A(i,j)=+(int)(50*(drand48()));
		if(j==3)A(i,j)=0;
		++k;
		//A(i,j)=k;
		}
		}

	Matrix B=A;
	cerr<< A <<endl;
	cerr<< "--------" <<endl;
	//cerr<< !A*v <<endl;
	//A=GaussEliminate2(A,v);
 	to_reduced_row_echelon_form(A,v);
	//A.swapRow(2,3);
	//v.swapRow(2,3);
	cerr<< "--------" <<endl;
	cerr<< A <<endl;
	cerr<< "--------" <<endl;
	cerr<< v <<endl;
	
	}


void mysort(int val[], int ind[] , int n){

//const int N=100;
//int ind[N];//for holding indeces 
//int val[N];//for holding indeces 
//for(int i=0; i<N; ++i){
	//ind[i]=i;
	//val[i]=100*drand48();
	//}

	bool swapped;

	do{
		swapped = false;
		for(int i=0; i<n-1; ++i){
		      if (val[ind[i]] > val[ind[i+1]]){
			swap( ind[i], ind[i+1] );
			swapped = true;
			}
			}
		--n;
	}while (swapped);


//	for(int i=0; i<N; ++i){
//	cerr<< val[ind[i]] <<endl;
//	}

}


int main(int argc, char **argv){
long seed=0;
MTRand rgen(seed);
CStat stat;
stat.set_hist(-3, 8, 500);

for(int i=0; i<10000000; i++){
	stat<<rgen.randNorm(2.5, 0.5);
	}
stat.print_hist("hist.dat");

return 0;
}

