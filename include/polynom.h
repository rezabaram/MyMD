#ifndef QUARTIC_H
#define QUARTIC_H 
#include<fstream>
#include<iomanip>
#include<math.h>
#include<vector>
#include<math.h>
#include<complex>
#include"exception.h"

using namespace std;


inline 
double myabs2(const complex<double> &c){
	return sqrt(c.imag()*c.imag()+c.real()*c.real());
	}

inline 
double myabs(const complex<double> &c){
	return sqrt(myabs2(c));
	}


//this is just to impose an interface and to reuse common features. 
//but not to make all polynomials accessible via a common pointer (which is usually the goal of polymorphism)

const double epsilon=1e-8;

template<int order, typename T=double>
class CPolynom{
	public:
  	CPolynom(const T a): solved(false){//just to enforce the first coefs be non-zero
		ERROR(fabs(a)<epsilon, "The leading coefficient cannot be zero");
		coefs.push_back(a);
		for(size_t i=0; i<order; ++i){
			coefs.push_back(T(0.0));
			}//total no. of coefs = order+1
		}
		
  	CPolynom(const vector<T> _coefs): solved(false){
		ERROR(_coefs.size()!=order+1, "Polynomial not created. The number of coefficents does not correspond the order.")
		ERROR(fabs(_coefs.at(0))<epsilon, "Polynomial not created. The coefficent of highest term should not be zero. ");

		coefs=_coefs;
		};

	virtual ~CPolynom(){}
	bool is_root(complex<double> x, double e=epsilon){
		return (myabs2((*this)(x))<e);
		}

	bool check_roots(){
		if(!solved)return false;
		bool b=true;
		for(int i=0; i<roots.size(); ++i){
			b= (b and is_root(roots.at(i)));
			}
		return b;
		}

	
	void print(ostream &out){
		
		out<<"[";
		for(size_t i=0; i<coefs.size()-1; ++i){
			out<<coefs.at(i)<<" , ";
			}
		out<<coefs.back()<<"]\n";
		}

	void print_roots(ostream &out=std::cerr){
		if(!solved and !solve()){
			WARNING("No root to print out; Polynomial has not been solved.")
			return;
			}
		
		out<<"Roots: [";
		for(size_t i=0; i<roots.size()-1; ++i){
			out<<roots.at(i)<<" , ";
			}
		out<<roots.back()<<"]\n";
		}

	T operator()(T x)const;
	complex<double> operator()(complex<double> x)const;
	complex<double> root(size_t i);
	virtual bool solve(){
		ERROR(true, "Solve method not implemented for this polynomial.");
		return -1;
		};

	void plot(ostream &out, T min, T max, T dx)const;
	

	protected:
	vector<T> coefs;
	vector<complex<double> > roots;
	bool solved;
	};

template<int order, typename T>
inline T CPolynom<order, T>::operator() (T x)const{
	static double xx, term, result;
	result=0;
	xx=1.0;
	
	for(size_t i=0; i<order+1; ++i){
		result+=xx*coefs.at(order-i);
		xx=x*xx;
		}
	return result;
	}

template<int order, typename T>
inline complex<double> CPolynom<order, T>::operator()(complex<double> x)const{
	static complex<double> xx, term, result;
	result=0;
	xx=1.0;
	
	for(size_t i=0; i<order+1; ++i){
		result+=xx*coefs.at(order-i);
		xx=x*xx;
		}
	return result;
	}

template<int order, typename T>
inline complex<double> CPolynom<order, T>::root(size_t i){
	ERROR(!solved and !solve(), "The polynomial has no roots;");//FIXME write the correct error message
	ERROR(i>=roots.size(), "Index out of range.");

	return roots.at(i);
	}

template<int order, typename T>
void CPolynom<order, T>::plot(ostream &out, T min, T max, T dx)const{

	for(T x=min; x<=max; x+=dx){
		out<<setprecision(12)<<x<<"  "<<(*this)(x)<<endl;
		}
	}


//-------------------- Cubic polynomial--------

class CQuadratic: public CPolynom<2, double> {
	public:
  	CQuadratic(const vector<double> _coefs)
		:CPolynom<2,double>(_coefs)
		{}
	CQuadratic(double _a, double _b, double _c);

	double max_root();
	bool solve();
	private:
	};

CQuadratic::CQuadratic(double _a, double _b, double _c)
	:CPolynom<2,double>(_a){
		coefs.at(0)=_a;
		coefs.at(1)=_b;
		coefs.at(2)=_c;
		}

inline
bool CQuadratic::solve(){//returning the number of real roots;
	if(solved)return true;
	static double delta;
	static double a, b, c;
	a=coefs.at(0);
	b=coefs.at(1);
	c=coefs.at(2);

	delta = b*b - 4.0*a*c;

        if(delta>=0){// two real roots.
		roots.push_back(complex<double> ( (-b-sqrt(delta))/2/a, 0));
		roots.push_back(complex<double> ( (-b+sqrt(delta))/2/a, 0));
	  	}

	else{//two complex roots
		roots.push_back(complex<double> (-b/2/a,-sqrt(-delta)/2/a));
		roots.push_back(complex<double> (-b/2/a,sqrt(-delta)/2/a));
		}
	
	solved=true;
	return true;
	}
inline
double CQuadratic::max_root(){
	if(!solved)solve();

	if( fabs(roots.back().imag() > epsilon)){
		WARNING("Returning only real part of a complex root.");
		}
	return roots.back().real();
	}

//-------------------- Cubic polynomial--------

class CCubic: public CPolynom<3, double> {
	public:
  	CCubic(const vector<double> _coefs)
		:CPolynom<3,double>(_coefs)
		{}
	CCubic(double _a, double _b, double _c, double _d);

	bool solve();
	private:
	};


CCubic::CCubic(double _a, double _b, double _c, double _d)
	:CPolynom<3,double>(_a){
		coefs.at(0)=_a;
		coefs.at(1)=_b;
		coefs.at(2)=_c;
		coefs.at(3)=_d;
		}

inline
bool CCubic::solve(){//returning the number of real roots;
	if(solved)return true;
	static double p, q, D;
	static double a, b, c, d;
	a=coefs.at(0);
	b=coefs.at(1);
	c=coefs.at(2);
	d=coefs.at(3);

	double ap=1.0/a;
	double ap2=ap*ap;
	double ap3=ap2*ap;
	double b2=b*b;
	double b3=b*b2;
	
	p =  +c*ap - b2*ap2 / 3.0;
        q = ( 2.0*b3*ap3 - 9.0*b*c*ap2 + 27*d*ap ) / 27.0;
        D = p*p*p/27 + q*q/4;
        if(D>=0){// three real roots.
		double temp1=-q/2 + sqrt(D);
		double temp2=-q/2 - sqrt(D);
		double u = pow(fabs(temp1), 1.0/3.0);
		double v = pow(fabs(temp2), 1.0/3.0);
		if(temp1<0)u=-u;
		if(temp2<0)v=-v;
		//here are the roots
		// y1 = u + v
		// y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
		// y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
		roots.push_back( complex<double>( -b/a/3+(u+v), .0) );
		roots.push_back( complex<double>( -b/a/3-(u+v)/2.0 , (u-v)*sqrt(3.)/2.0 ) );
		roots.push_back( complex<double>( -b/a/3-(u+v)/2.0 ,-(u-v)*sqrt(3.)/2.0 ) );
	  	}

	else{//three distinct real roots
		double phi = acos(-q/2.0/sqrt(fabs(p*p*p)/27));
		roots.push_back( complex<double>(-b/a/3 + 2 * sqrt(fabs(p)/3) * cos(phi/3), 0) );
		roots.push_back( complex<double>(-b/a/3 - 2 * sqrt(fabs(p)/3) * cos((phi+M_PI)/3), 0) );
		roots.push_back( complex<double>(-b/a/3 - 2 * sqrt(fabs(p)/3) * cos((phi-M_PI)/3), 0) );
		}
	
	solved=true;
	return true;
	}

//-------------------- Quartic polynomial--------

class CQuartic: public CPolynom<4, double> {
	public:
  	CQuartic (const vector<double> _coefs)
		:CPolynom<4,double>(_coefs)
		{}
	CQuartic (double _a, double _b, double _c, double _d, double _e);

	bool solve();
	private:
	};


CQuartic::CQuartic (double _a, double _b, double _c, double _d, double _e)
	:CPolynom<4,double>(_a){
		coefs.at(0)=_a;
		coefs.at(1)=_b;
		coefs.at(2)=_c;
		coefs.at(3)=_d;
		coefs.at(4)=_e;
		}


inline
bool CQuartic::solve(){//returning the number of real roots;
	if(solved)return true;

	//solution from wiki, Ferrari's method: http://en.wikipedia.org/wiki/Quartic_function#Solving_a_quartic_equation
	double ap=1.0/coefs.at(0);
	//double ap2=ap*ap;
	//double ap3=ap2*ap;
	//double ap4=ap2*ap2;
	static double a, b, c, d, e;
	a=coefs.at(0);
	b=coefs.at(1)/a;
	c=coefs.at(2)/a;
	d=coefs.at(3)/a;
	e=coefs.at(4)/a;
	a=1;

	double b2=b*b;
	double b3=b*b2;
	double b4=b2*b2;

	double f = c - (3*b2/8); 
	double g = d + (b3 / 8) - (b*c/2); //  d + (b3 / 8) - (b*c/2)
	double h = e - (3*b4/256.0) + (b2 * c/16.0) - (b*d/4);

	if(fabs(f) < epsilon and fabs(g) < epsilon and fabs(h)<epsilon){//this happens for (x-a)^4=0
		roots.push_back(-b/4);
		roots.push_back(-b/4);
		roots.push_back(-b/4);
		roots.push_back(-b/4);
		return true;
		}

	// Y3 + (f/2)*Y2 + ((f2 -4*h)/16)*Y -g2/64 = 0
	CCubic cube(1,(f/2),(f*f -4*h)/16, -g*g/64);
	
	complex<double> r1=cube.root(0);
	complex<double> r2=cube.root(1);
	complex<double> r3=cube.root(2);

	//choosing two non-zero roots
	int i1=-1;
	complex<double> p, q;
	if(myabs(cube.root(0))>0){
		    p=sqrt(cube.root(0));
		i1=0;
		}
	else if(myabs(cube.root(1))>0){
		p=sqrt(cube.root(1));
		i1=1;
		}
	else if(myabs(cube.root(2))>0){
		p=sqrt(cube.root(2));
		i1=2;
		}
	else ERROR(true, "Impossible happened");
		

	if(myabs(cube.root(0))>0 and i1!=0){
		q=sqrt(cube.root(0));
		}
	else if(myabs(cube.root(1))>0 and i1!=1){
		q=sqrt(cube.root(1));
		}
	else if(myabs(cube.root(2))>0 and i1!=2){
		q=sqrt(cube.root(2));
		}
	else  ERROR(true, "Impossible happened");

	complex<double> r= -g/(8.0*p*q);
	complex<double> s= ap*coefs.at(1)/4;

	roots.push_back( p + q + r -s);
	roots.push_back( p - q - r -s);
	roots.push_back(-p + q - r -s);
	roots.push_back(-p - q + r -s);

	solved=true;
	return true;
/*
	double alfa, beta, gamma;
	alfa=-3.0/8.0*b2*ap2 +coefs.at(2)*ap;
	beta=1.0/8.0*b3*ap3 - 0.5*coefs.at(1)*coefs.at(2)*ap2 + coefs.at(3)*ap;
	gamma=-3.0/256.0*b4*ap4 + 1.0/16.0*coefs.at(2)*b2*ap3 - 1.0/4* coefs.at(1)*coefs.at(3)*ap2 + coefs.at(4)*ap;

	if(fabs(beta)<epsilon){// the case beta=0. 
		double delta1=alfa*alfa-4*gamma;

		if(delta1<0)return 0;//no roots

		delta1=sqrt(delta1);
		double delta2=-alfa+delta1;
		double delta3=alfa*alfa-4*gamma;

		roots.push_back(0);
		}
*/
	}

#endif
