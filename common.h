#ifndef COMMON_H
#define COMMON_H 


//#define ERROR(x)  std::cerr<<"Error: In file " __FILE__<<" line "<<__LINE__<<":  "<<x<<std::endl;
//#include"matrix.h"

#include<assert.h>
#include<time.h>
#include<memory>
#include<limits>
#include<iomanip>
#include<sstream>
#include<ostream>
#include<vector>
#include<list>
#include"exception.h"
#include"config.h"

//#include"vec3d.h"
typedef size_t indexType;
#include"vec.h"

//#include"vec3d_policy.h"
//#include"vec3d.h"
//typedef vec3d<double> vec;

#include"quaternion.h"
#include"MersenneTwister.h"
#include"matrix.h"
using namespace math;
typedef matrix<double> Matrix;
#include"polynom.h"
#include"shapes.h"

#define FROMTIME double CLOCKSTART=clock(); cerr<<"Time in function "<<__FUNCTION__<<": ";
#define TOTIME cerr<<(clock()-CLOCKSTART)/CLOCKS_PER_SEC<<endl;

ofstream *gout=NULL;

extern CConfig &config; // don't forget "&" or you get a vicious bug, which took me one day to find

//typedef vec3d<double> vec;

double aG[]={0.0, 0, -2.0};
vec G(aG);

void define_parameters()
{
	config.add_param<vec>("Gravity", G);
	config.add_param<double>("outDt", 0.02);
	config.add_param<double>("stiffness", 5.0e+05); 
	config.add_param<double>("damping", 1); 
	config.add_param<double>("friction", 0); 
	config.add_param<double>("cohesion", 0); 
	config.add_param<double>("density", 10000.0); 
	config.add_param<double>("particleSize", 0.05); 
	config.add_param<double>("timeStep", 0.00001); 
	config.add_param<double>("maxTime", 4.0); 
	config.add_param<double>("verletfactor", 0.5); 
	config.add_param<size_t>("nParticle", 5); 
	config.add_param<string>("particleShape", "sphere"); 
	config.add_param<double>("e", 0.5); 
	config.add_param<double>("asphericity", -0.5); 
	config.add_param<double>("asphericityWidth", 0.1); 
}

template <class T>
string stringify(T x, int width=15, const char ch=' ')
 {
   std::ostringstream o;
   if (!(o << setw(width)<<setfill(ch)<<x))
     cerr<<"Bad coversion to string"<<endl;
   return o.str();
 }

using namespace math;


vec randomVec(const vec &x1, const vec &x2){
	return vec( x1(0)+(x2(0)-x1(0))*rgen(), x1(1)+(x2(1)-x1(1))*rgen(), x1(0)+(x2(2)-x1(2))*rgen());
	}

//obtaining the polynomial | B + lambda A | = 0
CQuartic characteristicPolynom(const Matrix &AB){

    	return CQuartic(
		1,

		 -AB.Tr(),

-(AB(0,1)*AB(1,0)) + AB(0,0)*AB(1,1) - AB(0,2)*AB(2,0) - AB(1,2)*AB(2,1) + AB(0,0)*AB(2,2) + AB(1,1)*AB(2,2) - AB(0,3)*AB(3,0) - AB(1,3)*AB(3,1) - AB(2,3)*AB(3,2) + AB(0,0)*AB(3,3) + AB(1,1)*AB(3,3) + AB(2,2)*AB(3,3),

AB(0,2)*AB(1,1)*AB(2,0) - AB(0,1)*AB(1,2)*AB(2,0) - AB(0,2)*AB(1,0)*AB(2,1) + AB(0,0)*AB(1,2)*AB(2,1) + AB(0,1)*AB(1,0)*AB(2,2) - AB(0,0)*AB(1,1)*AB(2,2) + AB(0,3)*AB(1,1)*AB(3,0) - AB(0,1)*AB(1,3)*AB(3,0) + AB(0,3)*AB(2,2)*AB(3,0) - 
   AB(0,2)*AB(2,3)*AB(3,0) - AB(0,3)*AB(1,0)*AB(3,1) + AB(0,0)*AB(1,3)*AB(3,1) + AB(1,3)*AB(2,2)*AB(3,1) - AB(1,2)*AB(2,3)*AB(3,1) - AB(0,3)*AB(2,0)*AB(3,2) - AB(1,3)*AB(2,1)*AB(3,2) + AB(0,0)*AB(2,3)*AB(3,2) + AB(1,1)*AB(2,3)*AB(3,2) + AB(0,1)*AB(1,0)*AB(3,3) - AB(0,0)*AB(1,1)*AB(3,3) + AB(0,2)*AB(2,0)*AB(3,3) + AB(1,2)*AB(2,1)*AB(3,3) - AB(0,0)*AB(2,2)*AB(3,3) - AB(1,1)*AB(2,2)*AB(3,3),

		AB.Det() ); 
}


Matrix  GaussEliminate (Matrix M, Matrix &v){
size_t i = 0, j = 0, maxi;
size_t m = M.RowNo(), n = M.ColNo();

while (i < m and j < n) {
	//Find pivot in column j, starting in row i:
	maxi = i;
	for (indexType k = i+1; k<m; ++k){
		if (fabs(M(k,j)) > fabs(M(maxi,j))){
			maxi = k;
			}
		}

	if( fabs(M(maxi,j)) > epsilon){
		M.swapRow(i, maxi);
		v.swapRow(i, maxi);
		//divide each entry in row i by M[i,j]
		double Mij=M(i,j);
		for(indexType jj=0; jj<n; ++jj){
			M(i,jj)/=Mij;
			v(i,0)/=Mij;
			}

		for(indexType u = i+1; u< m; ++u){
			// subtract A[u,j] * row i from row u
			double Muj=M(u,j);
			for(indexType jj=0; jj<n; ++jj){
				M(u,jj)-=Muj*M(i,jj);
				}
			v(u,0)-=Muj*v(i,0);
			}
		++i;
		}
	++j;

if(i>1)return M;
}

/*
for(indexType k=m-1; k<m; ++k){
	if(M(k,k)
      for(int jj=0; jj<n; ++jj){
		M(u,jj)=M(u,jj)-Muj*M(i,jj);
		}
	v(u,0)=v(u,0)-Muj*v(i,0);
	}
*/

return M;

}

//function ToReducedRowEchelonForm(Matrix M) is

void to_reduced_row_echelon_form(Matrix& A, Matrix &v) {
 
  typedef size_t index_type;
size_t max_row = A.RowNo()-1, max_column = A.ColNo()-1;

  index_type lead = 0;
 
  for (index_type row = 0; row <= max_row; ++row)
  {
    if (lead > max_column)
      return;
    index_type i = row;
    while (fabs(A (i, lead)) <epsilon)
    {
      ++i;
      if (i > max_row)
      {
        i = row;
        ++lead;
        if (lead >max_column)
          return;
      }
    }

    A.swapRow(i, row);
    v.swapRow(i, row);
    double Mrl=A(row, lead);
    A.multiply_row( row, 1/Mrl);
    v.multiply_row( row, 1/Mrl);
    for (i = 0; i <= max_row; ++i) {
      if (i != row){
	double Mil=-A(i, lead);
	for (size_t col = 0; col <=max_column; ++col) A (i, col) += Mil * A(row, col);
	//A(i,lead)=0;
	v (i, 0) += Mil * v(row, 0);
	}
        //add_multiple_row(A, i, row, -mt.element(A, i, lead));
    }
  }
}
 
#endif /* COMMON_H */
