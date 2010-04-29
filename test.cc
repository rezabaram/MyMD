#include<iostream>
using namespace std;
#include"polynom.h"
#include"group.h"
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

void listtest(){
	CList<int> list;
	for(int i=0; i<100; ++i){
	list.push_back(i);
	}

	}
void gaussElimTest(){
	Matrix M(3,3);
	
	M(0,0)=0;
	M(0,1)=1;
	M(0,2)=-1;
	M(1,0)=-3;
	M(1,1)=-1;
	M(1,2)=2;
	M(2,0)=-2;
	M(2,1)=1;
	M(2,2)=2;
	vec v(8, -11, -3);

	cerr<< M <<endl;
	M=GaussEliminate(M, v);
	cerr<< M <<endl;
	}

int main(int argc, char **argv){
gaussElimTest();
return 0;
}
