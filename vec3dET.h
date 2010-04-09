//=======================================================================
//======  Class sVec3d Template  =====================================
//=======================================================================
#ifndef Vec3d_ET_H
#define Vec3d_ET_H
#include<ostream>
#include<math.h>


template<typename, typename>
class vec3d;

typedef unsigned short indexType;
template<class T>
class sVec3d{
 friend class vec3d<T, sVec3d>;
 public:
  sVec3d(const T vec[]){
	x[0]=vec[0]; x[1]=vec[1]; x[2]=vec[2];
	};
   sVec3d(T x0, T x1, T x2){
	x[0]=x0; x[1]=x1; x[2]=x2;
	}


sVec3d(T a){
	x[0]=(T)a;
	x[1]=(T)a;
	x[2]=(T)a;
	};


  sVec3d(){
	x[0]=(T)0;
	x[1]=(T)0;
	x[2]=(T)0;
	};
  ~sVec3d(){
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

  T operator*(const sVec3d<T> &p)const;
  //const sVec3d<T> operator*(T a)const ;
  sVec3d<T> &operator=(T a);
  sVec3d<T> &operator*=(T a);
  sVec3d<T> &operator/=(T a);
  sVec3d<T> &operator+=(T a);
  sVec3d<T> &operator-=(T a);
  sVec3d<T> &operator+=(const sVec3d<T> &p);
  sVec3d<T> &operator-=(const sVec3d<T> &p);
  sVec3d<T> operator*(T a)const ;
  sVec3d<T> operator/(T a)const ;
  sVec3d<T> operator+(const sVec3d<T> &p)const ;
  sVec3d<T> operator-(const sVec3d<T> &p)const ;
  T abs2()const;
  T abs()const;
  T distance(const sVec3d<T> p)const;
  sVec3d<T> & normalize();
  sVec3d<T> & normalized()const;

  template<class U>
  friend std::ostream & operator<< (std::ostream &out, const sVec3d<U> &v);
  template<class U>
  friend std::istream & operator>>(std::istream &in, sVec3d<U> &v);
private:
  T x[3];
};



template<class T>
  std::ostream & operator<< (std::ostream &out, const sVec3d<T> &v){
	out<<v.x[0]<<"  "<<v.x[1]<<"  "<<v.x[2];
	return out;
	}

template<class T>
  std::istream & operator>>(std::istream &in, sVec3d<T> &v){
		in>>v.x[0]>>v.x[1]>>v.x[2];
	return in;
	}

template<class T>
  T sVec3d<T>::operator*(const sVec3d<T> &p)const{
	T prod;
	prod=x[0]*p.x[0]+ x[1]*p.x[1]+x[2]*p.x[2];
	return prod;
	}

template<class T>
  sVec3d<T> &sVec3d<T>::operator=(T a){
	x[0]=a; x[1]=a; x[2]=a;
	return *this;
	}

template<class T>
  sVec3d<T> &sVec3d<T>::operator*=(T a){
	x[0]*=a; x[1]*=a; x[2]*=a;
	return *this;
	}
template<class T>
  sVec3d<T> &sVec3d<T>::operator/=(T a){
	x[0]/=a; x[1]/=a; x[2]/=a;
	return *this;
	}
template<class T>
  sVec3d<T> &sVec3d<T>::operator+=(T a){
	x[0]+=a; x[1]+=a; x[2]+=a;
	return *this;
	}
template<class T>
  sVec3d<T> &sVec3d<T>::operator-=(T a){
	x[0]-=a; x[1]-=a; x[2]-=a;
	return *this;
	}

template<class T>
  sVec3d<T> &sVec3d<T>::operator+=(const sVec3d<T> &p){
	x[0]+=p.x[0];
	x[1]+=p.x[1];
	x[2]+=p.x[2];
	return *this;
	}

template<class T>
  sVec3d<T> &sVec3d<T>::operator-=(const sVec3d<T> &p){
	x[0]-=p.x[0];
	x[1]-=p.x[1];
	x[2]-=p.x[2];
	return *this;
	}

template<class T>
  sVec3d<T> sVec3d<T>::operator*(T a)const {
	sVec3d<T> pp(*this);
		pp*=(a);
	return pp;
	}


template<class T>
  sVec3d<T> sVec3d<T>::operator/(T a)const {
	sVec3d<T> pp(*this);
		pp/=(a);
	return pp;
	}

template<class T>
  sVec3d<T> sVec3d<T>::operator+(const sVec3d<T> &p)const {
	sVec3d<T> pp(*this);
	pp+=p;
	return pp;
	}


template<class T>
  sVec3d<T> sVec3d<T>::operator-(const sVec3d<T> &p)const {
	sVec3d<T> pp(*this);
	pp-=p;
	return pp;
	}

template<class T>
  T sVec3d<T>::abs2()const{
  	return (*this)*(*this);
	}

template<class T>
  T sVec3d<T>::abs()const{
  	return sqrt(abs2());
	}

template<class T>
  sVec3d<T> & sVec3d<T>::normalize(){
	T d=1./abs();
	x[0]*=d; x[1]*=d; x[2]*=d;
	return *this;
	}

template<class T>
  sVec3d<T> & sVec3d<T>::normalized()const{
	return sVec3d<T> (*this).normalize();
	}

template<class T>
  const sVec3d<T> operator *(double a, sVec3d<T> v){
	v*=(a);
	return v;
	};

template<class T>
  const sVec3d<T> operator -(sVec3d<T> v){
	v*=(T)(-1.0);
	return v;
	};

//template<>
//sVec3d::sVec3d(double , double, double){}

template<class T>
sVec3d<T> cross(const sVec3d<T> u, const sVec3d<T> v){
	return sVec3d<T> (u(1)*v(2)-u(2)*v(1),u(2)*v(0)-u(0)*v(2),u(0)*v(1)-u(1)*v(0));
	}


using namespace std;
template<typename T, typename Rep = sVec3d<T> >
class vec3d{
	private:
	Rep expr_rep;
	public:
	vec3d():expr_rep(Rep()){
		}
	vec3d(Rep const &v):expr_rep(v){
		}
	vec3d& operator= (vec3d const&v){
		expr_rep[0]=v[0]; expr_rep[1]=v[1]; expr_rep[2]=v[2];
		return *this;
		}
	template<typename T2, typename Rep2>
	vec3d& operator= (vec3d<T2, Rep2> const&v){
		expr_rep[0]=v[0]; expr_rep[1]=v[1]; expr_rep[2]=v[2];
		return *this;
		}
	T& operator [](size_t i){
		return expr_rep[i];
		}

	const T &operator [](size_t i)const{
		return expr_rep[i];
		}

	T& operator ()(size_t i){
		return expr_rep[i];
		}

	const T &operator ()(size_t i)const{
		return expr_rep[i];
		}

	Rep &rep(){
		return expr_rep; 
		}
	Rep const & rep()const {
		return expr_rep; 
		}

template<typename U, typename Rep2>
  	friend std::ostream & operator<< (std::ostream &out, const vec3d<U,Rep2> &v);
template<typename U, typename Rep2>
  	friend std::istream & operator>>(std::istream &in, vec3d<U,Rep2> &v);
	};
template<typename T, typename Rep>
 std::ostream & operator<< (std::ostream &out, const vec3d<T,Rep> &v){
	out<<v.expr_rep;
	return out;
	}

template<typename T, typename OPT1, typename OPT2>
class Add{
	OPT1 op1;
	OPT2 op2;
	public:
	Add(OPT1 const &a, OPT2 const &b):op1(a), op2(b) {}

	T operator [](size_t i)const{
		return op1[i]+op2[i];
		}
	};

template<typename T, typename OPT1, typename OPT2>
inline
vec3d<T, Add<T, OPT1, OPT2> > operator + (vec3d<T, OPT1> const &a, vec3d<T, OPT2> const &b ){
	return vec3d<T, Add<T,OPT1,OPT2> >(Add<T, OPT1, OPT2> (a.rep(), b.rep()) ) ;
	}

template<typename T, typename Rep>
 std::istream & operator>>(std::istream &in, vec3d<T,Rep> &v){
	in>>v.expr_rep;
	return in;
	}

#endif

