#ifndef EIGSYS_H
#define EIGSYS_H 
#include"exception.h"
#include"matrix.h"
#include"polynom.h"
#include"eigen.h"

class CEigSys{
	public:
	CEigSys(const Matrix &_M):M(_M), charPolynom(1,0,0,0),solved(false){
		//charPolynom=getCharPolynomial();
		}

	CCubic getCharPolynomial()const;
	void solve(){
		//charPolynom.solve();
		eigens(M,eigvals, eigvecs);
		solved=true;
		}
	void print_eigens(ostream &out){
		if(!solved)solve();
		//out<< real(charPolynom.root(0)) <<"\t"<< real(charPolynom.root(1)) <<"\t"<< real(charPolynom.root(2)) <<endl;
		for(unsigned int i=0; i<eigvals.size(); i++){
			out<<eigvals.at(i)<<"   "<<spherical(eigvecs.at(i))<<"   ";
			}
		}

	const Matrix &M;
	vector<double> eigvals;
	vector<vec3d> eigvecs;
	CCubic charPolynom;
	bool solved;
	};

CCubic CEigSys::getCharPolynomial()const{
	return CCubic( -1,
			M(0,0) + M(1,1) + M(2,2), 
			M(0,1)*M(1,0) - M(0,0)*M(1,1) + M(0,2)*M(2,0) + M(1,2)*M(2,1) - M(0,0)*M(2,2) - M(1,1)*M(2,2), 
			-M(0,2)*M(1,1)*M(2,0) + M(0,1)*M(1,2)*M(2,0) + M(0,2)*M(1,0)*M(2,1) - M(0,0)*M(1,2)*M(2,1) - M(0,1)*M(1,0)*M(2,2) + M(0,0)*M(1,1)*M(2,2)
			);
	}
#endif /* EIGSYS_H */
