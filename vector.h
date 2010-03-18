//=======================================================================
//======  Class CVector Template  =====================================
//=======================================================================
#ifndef CVECTOR_H
#define CVECTOR_H
#include<ostream>
#include<math.h>


 typedef unsigned short indexType;
template<indexType Dim, class T>
class CVector{
 public:
  CVector(const T vec[]){
	for(indexType i=0; i<Dim; i++){
		x[i]=vec[i];
		}
	};

/*FIXME  although in this case (using array, and not allocating memory dynamically) one doesn't need to
	 define the copy and assignment constructors, I could not define one which works propery.
	 How should the copy constructor of a template class be??!! 

  CVector(const CVector & v){ 
	for(indexType i=0; i<Dim; i++){
		x[i]=v(i);
		}
	}

  CVector(CVector<Dim, T> const & v){
	for(indexType i=0; i<Dim; i++){
		x[i]=v(i);
		}
	};

  void operator=(CVector const & v){
	for(indexType i=0; i<Dim; i++){
		x[i]=v.x[i];
		}
	};
*/

explicit CVector(T a){
	for(indexType i=0; i<Dim; i++){
		x[i]=(T)a;
		}
	};

  CVector(indexType j, T a){//set jth component, and the rest zero
	for(indexType i=0; i<Dim; i++){
		x[i]=(T)0;
		}
		x[j]=(T)a;
	};

  CVector(){
	for(indexType i=0; i<Dim; i++){
		x[i]=(T)0;
		}
	};
  ~CVector(){
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

  T operator*(const CVector<Dim, T> &p)const;
  //const CVector<Dim,T> operator*(T a)const ;
  CVector<Dim,T> &operator=(T a);
  CVector<Dim,T> &operator*=(T a);
  CVector<Dim,T> &operator/=(T a);
  CVector<Dim,T> &operator+=(T a);
  CVector<Dim,T> &operator-=(T a);
  CVector<Dim,T> &operator+=(const CVector<Dim, T> &p);
  CVector<Dim,T> &operator-=(const CVector<Dim, T> &p);
  CVector<Dim,T> operator*(T a)const ;
  CVector<Dim,T> operator/(T a)const ;
  CVector<Dim,T> operator+(const CVector<Dim, T> &p)const ;
  CVector<Dim,T> operator-(const CVector<Dim, T> &p)const ;
  T abs2()const;
  T abs()const;
  T distance(const CVector<Dim,T> p)const;
  CVector<Dim,T> & normalize();
  CVector<Dim,T> & normalized()const;

  template<indexType _Dim, class U>
  friend std::ostream & operator<< (std::ostream &out, const CVector<_Dim, U> &v);
  template<indexType _Dim, class U>
  friend std::istream & operator>>(std::istream &in, CVector<_Dim, U> &v);
private:
static const indexType dim=Dim;
  T x[Dim];
};



template<indexType Dim, class T>
  std::ostream & operator<< (std::ostream &out, const CVector<Dim,T> &v){
	out<<v.x[0];
	for(indexType i=1; i<Dim; i++){
		out<<"  "<<v.x[i];
		}
	return out;
	}

template<indexType Dim, class T>
  std::istream & operator>>(std::istream &in, CVector<Dim,T> &v){
	for(indexType i=0; i<Dim; i++){
		in>>v.x[i];
		}
	return in;
	}

template<indexType Dim, class T>
  T CVector<Dim,T>::operator*(const CVector<Dim, T> &p)const{
	T prod;
	prod=0;
	for(indexType i=0; i<Dim; i++){
		prod+=x[i]*p.x[i];
		}
	return prod;
	}

template<indexType Dim, class T>
  CVector<Dim,T> &CVector<Dim, T>::operator=(T a){
	for(indexType i=0; i<Dim; i++){
		x[i]=a;
		}
	return *this;
	}

template<indexType Dim, class T>
  CVector<Dim,T> &CVector<Dim, T>::operator*=(T a){
	for(indexType i=0; i<Dim; i++){
		x[i]*=a;
		}
	return *this;
	}
template<indexType Dim, class T>
  CVector<Dim,T> &CVector<Dim, T>::operator/=(T a){
	for(indexType i=0; i<Dim; i++){
		x[i]/=a;
		}
	return *this;
	}
template<indexType Dim, class T>
  CVector<Dim,T> &CVector<Dim, T>::operator+=(T a){
	for(indexType i=0; i<Dim; i++){
		x[i]+=a;
		}
	return *this;
	}
template<indexType Dim, class T>
  CVector<Dim,T> &CVector<Dim, T>::operator-=(T a){
	for(indexType i=0; i<Dim; i++){
		x[i]-=a;
		}
	return *this;
	}

template<indexType Dim, class T>
  CVector<Dim,T> &CVector<Dim, T>::operator+=(const CVector<Dim, T> &p){
	for(indexType i=0; i<Dim; i++){
		x[i]+=p.x[i];
		}
	return *this;
	}

template<indexType Dim, class T>
  CVector<Dim,T> &CVector<Dim, T>::operator-=(const CVector<Dim, T> &p){
	for(indexType i=0; i<Dim; i++){
		x[i]-=p.x[i];
		}
	return *this;
	}

template<indexType Dim, class T>
  CVector<Dim,T> CVector<Dim, T>::operator*(T a)const {
	CVector<Dim,T> pp(*this);
		pp*=(a);
	return pp;
	}


template<indexType Dim, class T>
  CVector<Dim,T> CVector<Dim, T>::operator/(T a)const {
	CVector<Dim,T> pp(*this);
		pp/=(a);
	return pp;
	}

template<indexType Dim, class T>
  CVector<Dim,T> CVector<Dim, T>::operator+(const CVector<Dim, T> &p)const {
	CVector<Dim,T> pp(*this);
	pp+=p;
	return pp;
	}


template<indexType Dim, class T>
  CVector<Dim,T> CVector<Dim, T>::operator-(const CVector<Dim, T> &p)const {
	CVector<Dim,T> pp(*this);
	pp-=p;
	return pp;
	}

template<indexType Dim, class T>
  T CVector<Dim,T>::abs2()const{
  	return (*this)*(*this);
	}

template<indexType Dim, class T>
  T CVector<Dim,T>::abs()const{
  	return sqrt(abs2());
	}

template<indexType Dim, class T>
  CVector<Dim,T> & CVector<Dim,T>::normalize(){
	T d=1./abs();
	for(indexType i=0; i<Dim; i++){
		x[i]*=d;
		}
	return *this;
	}
template<indexType Dim, class T>
  CVector<Dim,T> & CVector<Dim,T>::normalized()const{
	return CVector<Dim,T> (*this).normalize();
	}

template<indexType Dim, class T>
  const CVector<Dim,T> operator *(double a, CVector<Dim,T> v){
	v*=(a);
	return v;
	};

template<indexType Dim, class T>
  const CVector<Dim,T> operator -(CVector<Dim,T> v){
	v*=(T)(-1.0);
	return v;
	};

typedef CVector<4,double> vec4d;
typedef CVector<3,double> vec3d;
typedef CVector<2,double> vec2d;

//template<>
//vec3d::CVector(double , double, double){}

vec3d cross(const vec3d u, const vec3d v){
	double crossprod[]={u(1)*v(2)-u(2)*v(1),u(2)*v(0)-u(0)*v(2),u(0)*v(1)-u(1)*v(0)};
	return vec3d(crossprod);
	}

#endif

