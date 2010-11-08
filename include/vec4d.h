//=======================================================================
//======  Class vec4d Template  =====================================
//=======================================================================
#ifndef CVEC4D_H
#define CVEC4D_H
#include<ostream>
#include<math.h>
#include<complex>
using namespace std;

template<class T>
class vecElemPolicy{
	public:
	 static T mul (const T &a, const T &b){ return a*b;}
};
template<>
class vecElemPolicy<complex<double> >{
	public:
	 static complex<double> mul (const complex<double> &a, const complex<double> &b){ 
			return a*conj(b);}
};


template<class T, typename policy = vecElemPolicy<T> >
class vec4d{
 public:
  vec4d(const T vec[]){
	for(indexType i=0; i<4; i++){
		x[i]=vec[i];
		}
	};
   vec4d(const T &x0, const T &x1, const T &x2, const T &x3=0){
	x[0]=x0; 
	x[1]=x1; 
	x[2]=x2;
	x[3]=x3;
	}


explicit vec4d(T a){
	x[0]=(T)a;
	x[1]=(T)a;
	x[2]=(T)a;
	x[3]=(T)a;
	};

  vec4d(){
	x[0]=(T)0;
	x[1]=(T)0;
	x[2]=(T)0;
	x[3]=(T)0;
	};

  ~vec4d(){
	}

  T &operator[](indexType i){
	return x[i];
	};

  const T &operator[](indexType i)const{
	return x[i];
	};

  T &operator()(indexType i){
	return x[i];
	};
  const T &operator()(indexType i)const {
	return x[i];
	};

  T operator*(const vec4d<T, policy> &p)const;
  vec4d<T, policy> operator^(const vec4d<T, policy> &p)const;
  //const vec4d<T, policy> operator*(T a)const ;
  vec4d<T, policy> &operator=(T a);
  vec4d<T, policy> &operator*=(T a);
  vec4d<T, policy> &operator/=(T a);
  vec4d<T, policy> &operator+=(T a);
  vec4d<T, policy> &operator-=(T a);
  vec4d<T, policy> &operator+=(const vec4d<T, policy> &p);
  vec4d<T, policy> &operator-=(const vec4d<T, policy> &p);
  vec4d<T, policy> operator*(T a)const ;
  vec4d<T, policy> operator/(T a)const ;
  vec4d<T, policy> operator+(const vec4d<T, policy> &p)const ;
  vec4d<T, policy> operator-(const vec4d<T, policy> &p)const ;
  T abs2()const;
  T abs()const;
  T distance(const vec4d<T, policy> p)const;
  vec4d<T, policy> & normalize();
  vec4d<T, policy> & normalized()const;

  template<class U, class p >
  friend std::ostream & operator<< (std::ostream &out, const vec4d<U, p> &v);
  template<class U, class p>
  friend std::istream & operator>>(std::istream &in, vec4d<U, p> &v);
private:
  T x[4];
  static const indexType dim;
};
template<class T, class policy>
const indexType vec4d<T, policy>::dim=4;



template<class T, class policy>
inline
  std::ostream & operator<< (std::ostream &out, const vec4d<T, policy> &v){
	out<<v.x[0]<<"  "<<v.x[1]<< "  "<<v.x[2]<< "  "<<v.x[3];
	return out;
	}

template<class T, class policy>
  std::istream & operator>>(std::istream &in, vec4d<T, policy> &v){
		in>>v.x[0]>>v.x[1] >>v.x[2] >>v.x[3];
	return in;
	}

template<class T, class policy>
vec4d<T, policy> vec4d<T, policy>::operator^(const vec4d<T, policy> &p)const{
	return vec4d<T, policy>(x[0]*p.x[0], x[1]*p.x[1], x[2]*p.x[2], x[3]*p.x[3]);
	}

template<class T, class policy>
inline
  T vec4d<T, policy>::operator*(const vec4d<T, policy> &p)const{
	T prod;
	//prod=x[0]*p.x[0]+ x[1]*p.x[1]+x[2]*p.x[2];
	prod=policy::mul(x[0],p.x[0])+ policy::mul(x[1],p.x[1]) +policy::mul(x[2],p.x[2]) +policy::mul(x[3],p.x[3]);
	return prod;
	}

template<class T, class policy>
inline
  vec4d<T, policy> &vec4d<T, policy>::operator=(T a){
	x[0]=a; x[1]=a; x[2]=a; x[3]=a;
	return *this;
	}

template<class T, class policy>
inline
  vec4d<T, policy> &vec4d<T, policy>::operator*=(T a){
	x[0]*=a; x[1]*=a; x[2]*=a; x[3]*=a;
	return *this;
	}
template<class T, class policy>
inline
  vec4d<T, policy> &vec4d<T, policy>::operator/=(T a){
	x[0]/=a; x[1]/=a; x[2]/=a; x[3]/=a;
	return *this;
	}
template<class T, class policy>
inline
  vec4d<T, policy> &vec4d<T, policy>::operator+=(T a){
	x[0]+=a; x[1]+=a; x[2]+=a; x[3]+=a;
	return *this;
	}
template<class T, class policy>
inline
  vec4d<T, policy> &vec4d<T, policy>::operator-=(T a){
	x[0]-=a; x[1]-=a; x[2]-=a; x[3]-=a;
	return *this;
	}

template<class T, class policy>
inline
  vec4d<T, policy> &vec4d<T, policy>::operator+=(const vec4d<T, policy> &p){
	x[0]+=p.x[0];
	x[1]+=p.x[1];
	x[2]+=p.x[2];
	x[3]+=p.x[3];
	return *this;
	}

template<class T, class policy>
inline
  vec4d<T, policy> &vec4d<T, policy>::operator-=(const vec4d<T, policy> &p){
	x[0]-=p.x[0];
	x[1]-=p.x[1];
	x[2]-=p.x[2];
	x[3]-=p.x[3];
	return *this;
	}

template<class T, class policy>
inline
  vec4d<T, policy> vec4d<T, policy>::operator*(T a)const {
	vec4d<T, policy> pp(*this);
		pp*=(a);
	return pp;
	}


template<class T, class policy>
inline
  vec4d<T, policy> vec4d<T, policy>::operator/(T a)const {
	vec4d<T, policy> pp(*this);
		pp/=(a);
	return pp;
	}

template<class T, class policy>
inline
  vec4d<T, policy> vec4d<T, policy>::operator+(const vec4d<T, policy> &p)const {
	vec4d<T, policy> pp(*this);
	pp+=p;
	return pp;
	}


template<class T, class policy>
inline
  vec4d<T, policy> vec4d<T, policy>::operator-(const vec4d<T, policy> &p)const {
	vec4d<T, policy> pp(*this);
	pp-=p;
	return pp;
	}

template<class T, class policy>
inline
  T vec4d<T, policy>::abs2()const{
  	return (*this)*(*this);
	}

template<class T, class policy>
inline
  T vec4d<T, policy>::abs()const{
  	return sqrt(abs2());
	}

template<class T, class policy>
inline
  vec4d<T, policy> & vec4d<T, policy>::normalize(){
	T d=1./abs();
	x[0]*=d; x[1]*=d; x[2]*=d; x[3]*=d;
	return *this;
	}
template<class T, class policy>
inline
  vec4d<T, policy> & vec4d<T, policy>::normalized()const{
	return vec4d<T, policy> (*this).normalize();
	}

template<class T, class policy>
inline
  const vec4d<T, policy> operator *(double a, vec4d<T, policy> v){
	v*=(a);
	return v;
	};

template<class T, class policy>
inline
  const vec4d<T, policy> operator -(vec4d<T, policy> v){
	v*=(T)(-1.0);
	return v;
	};

//template<>
//vec4d::vec4d(double , double, double){}

template<class T, class policy>
inline
vec4d<T, policy> cross(const vec4d<T, policy> u, const vec4d<T, policy> v){
	//ignoring the last component
	return vec4d<T, policy> (policy::mul(u(1),v(2))-policy::mul(u(2),v(1)),policy::mul(u(2),v(0))-policy::mul(u(0),v(2)),policy::mul(u(0),v(1))-policy::mul(u(1),v(0)), 0);
	}

#endif

