#ifndef EIGSYS_H
#define EIGSYS_H 
#include"exception.h"
#include"matrix.h"
#include"polynom.h"

class CEigSys{
	public:
	CEigSys(const Matrix &_M):M(_M), charPolynom(1,0,0,0){
		charPolynom=getCharPolynomial();
		}

	CCubic getCharPolynomial()const;
	void solve(){
		charPolynom.solve();
		}
	void print_eigen_vals(ostream &out){
		out<< real(charPolynom.root(0)) <<"\t"<< real(charPolynom.root(1)) <<"\t"<< real(charPolynom.root(2)) <<endl;
		}

	const Matrix &M;
	CCubic charPolynom;
	};

CCubic CEigSys::getCharPolynomial()const{
	return CCubic( -1,
			M(0,0) + M(1,1) + M(2,2), 
			M(0,1)*M(1,0) - M(0,0)*M(1,1) + M(0,2)*M(2,0) + M(1,2)*M(2,1) - M(0,0)*M(2,2) - M(1,1)*M(2,2), 
			-M(0,2)*M(1,1)*M(2,0) + M(0,1)*M(1,2)*M(2,0) + M(0,2)*M(1,0)*M(2,1) - M(0,0)*M(1,2)*M(2,1) - M(0,1)*M(1,0)*M(2,2) + M(0,0)*M(1,1)*M(2,2)
			);
	}
#endif /* EIGSYS_H */
