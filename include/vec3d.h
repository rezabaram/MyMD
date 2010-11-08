//=======================================================================
//======  Class vec3d Template  =====================================
//=======================================================================
#ifndef CVEC3D_H
#define CVEC3D_H
#include<ostream>
#include<math.h>


template<class T>
class vec3d{
 public:
  vec3d(const T vec[]){
	for(indexType i=0; i<3; i++){
		x[i]=vec[i];
		}
	};
   vec3d(const T &x0, const T &x1, const T &x2){
	x[0]=x0; x[1]=x1; x[2]=x2;
	}


explicit vec3d(T a){
	x[0]=(T)a;
	x[1]=(T)a;
	x[2]=(T)a;
	};

  vec3d(){
	x[0]=(T)0;
	x[1]=(T)0;
	x[2]=(T)0;
	};

  ~vec3d(){
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

  T operator*(const vec3d<T> &p)const;
  vec3d<T> operator^(const vec3d<T> &p)const;
  //const vec3d<T> operator*(T a)const ;
  vec3d<T> &operator=(T a);
  vec3d<T> &operator*=(T a);
  vec3d<T> &operator/=(T a);
  vec3d<T> &operator+=(T a);
  vec3d<T> &operator-=(T a);
  vec3d<T> &operator+=(const vec3d<T> &p);
  vec3d<T> &operator-=(const vec3d<T> &p);
  vec3d<T> operator*(T a)const ;
  vec3d<T> operator/(T a)const ;
  vec3d<T> operator+(const vec3d<T> &p)const ;
  vec3d<T> operator-(const vec3d<T> &p)const ;
  T abs2()const;
  T abs()const;
  T distance(const vec3d<T> p)const;
  vec3d<T> & normalize();
  vec3d<T>  normalized()const;

  template<class U>
  friend std::ostream & operator<< (std::ostream &out, const vec3d<U> &v);
  template<class U>
  friend std::istream & operator>>(std::istream &in, vec3d<U> &v);
private:
  T x[3];
};



template<class T>
inline
  std::ostream & operator<< (std::ostream &out, const vec3d<T> &v){
	out<<v.x[0]<<"  "<<v.x[1]<<"  "<<v.x[2];
	return out;
	}

template<class T>
  std::istream & operator>>(std::istream &in, vec3d<T> &v){
		in>>v.x[0]>>v.x[1]>>v.x[2];
	return in;
	}

template<class T>
vec3d<T> vec3d<T>::operator^(const vec3d<T> &p)const{
	return vec3d<T>(x[0]*p[0], x[1]*p[1], x[2]*p[2]);
	}

template<class T>
inline
  T vec3d<T>::operator*(const vec3d<T> &p)const{
	T prod;
	prod=x[0]*p[0]+ x[1]*p[1]+x[2]*p[2];
	return prod;
	}

template<class T>
inline
  vec3d<T> &vec3d<T>::operator=(T a){
	x[0]=a; x[1]=a; x[2]=a;
	return *this;
	}

template<class T>
inline
  vec3d<T> &vec3d<T>::operator*=(T a){
	x[0]*=a; x[1]*=a; x[2]*=a;
	return *this;
	}
template<class T>
inline
  vec3d<T> &vec3d<T>::operator/=(T a){
	x[0]/=a; x[1]/=a; x[2]/=a;
	return *this;
	}
template<class T>
inline
  vec3d<T> &vec3d<T>::operator+=(T a){
	x[0]+=a; x[1]+=a; x[2]+=a;
	return *this;
	}
template<class T>
inline
  vec3d<T> &vec3d<T>::operator-=(T a){
	x[0]-=a; x[1]-=a; x[2]-=a;
	return *this;
	}

template<class T>
inline
  vec3d<T> &vec3d<T>::operator+=(const vec3d<T> &p){
	x[0]+=p[0];
	x[1]+=p[1];
	x[2]+=p[2];
	return *this;
	}

template<class T>
inline
  vec3d<T> &vec3d<T>::operator-=(const vec3d<T> &p){
	x[0]-=p[0];
	x[1]-=p[1];
	x[2]-=p[2];
	return *this;
	}

template<class T>
inline
  vec3d<T> vec3d<T>::operator*(T a)const {
	vec3d<T> pp(*this);
		pp*=(a);
	return pp;
	}


template<class T>
inline
  vec3d<T> vec3d<T>::operator/(T a)const {
	vec3d<T> pp(*this);
		pp/=(a);
	return pp;
	}

template<class T>
inline
  vec3d<T> vec3d<T>::operator+(const vec3d<T> &p)const {
	vec3d<T> pp(*this);
	pp+=p;
	return pp;
	}


template<class T>
inline
  vec3d<T> vec3d<T>::operator-(const vec3d<T> &p)const {
	vec3d<T> pp(*this);
	pp-=p;
	return pp;
	}

template<class T>
inline
  T vec3d<T>::abs2()const{
  	return (*this)*(*this);
	}

template<class T>
inline
  T vec3d<T>::abs()const{
  	return sqrt(abs2());
	}

template<class T>
inline
  vec3d<T> & vec3d<T>::normalize(){
	T d=1./abs();
	x[0]*=d; x[1]*=d; x[2]*=d;
	return *this;
	}
template<class T>
inline
  vec3d<T>  vec3d<T>::normalized()const{
	return vec3d<T> (*this).normalize();
	}

template<class T>
inline
  const vec3d<T> operator *(double a, vec3d<T> v){
	v*=(a);
	return v;
	};

template<class T>
inline
  const vec3d<T> operator -(vec3d<T> v){
	v*=(T)(-1.0);
	return v;
	};

//template<>
//vec3d::vec3d(double , double, double){}

template<class T>
inline
vec3d<T> cross(const vec3d<T> u, const vec3d<T> v){
	return vec3d<T> (u(1)*v(2)-u(2)*v(1),u(2)*v(0)-u(0)*v(2),u(0)*v(1)-u(1)*v(0));
	}

#endif

