//=======================================================================
//======  Class Vec Template  =====================================
//=======================================================================
#ifndef CVec_H
#define CVec_H
#include<ostream>
#include<math.h>
#include<assert.h>
#include<complex>

#define check_index

#ifdef check_index
#define CHECK_INDEX(i, D) assert(0<=i and i< D)
#endif 

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


template<indexType _dim=3, class T=double, typename policy = vecElemPolicy<T> >
class Vec{
 public:
  Vec(const T v[]){
	for(indexType i=0; i<dim; i++){
		x[i]=v[i];
		}
	};
	
        explicit Vec(const T &x0, const T &x1=(T)0, const T &x2=(T)0, const T &x3=(T)0){

		if(dim>0){
			x[0]=x0; 
			}
		else{
			return;
			}	
		if(dim>1){
			x[1]=x1; 
			}
		else{
			return;
			}	
		if(dim>2){
			x[2]=x2; 
			}
		else{
			return;
			}	
		if(dim>3){
			x[3]=x3; 
			}
		else{
			return;
			}	

		for(indexType i=4; i<dim; i++){
			x[i]=(T)0;
			}
		}

	  Vec(){
		for(indexType i=0; i<dim; i++){
			x[i]=(T)0;
			}
		};

	  ~Vec(){
		}

	  T &operator[](indexType i){
		CHECK_INDEX(i, dim);
		return x[i];
		};

	  const T &operator[](indexType i)const{
		CHECK_INDEX(i, dim);
		return x[i];
		};

	  T &operator()(indexType i){
		CHECK_INDEX(i, dim);
		return x[i];
		};
	  const T &operator()(indexType i)const {
		CHECK_INDEX(i, dim);
		return x[i];
		};

	  T operator*(const Vec<_dim, T, policy> &p)const;
	  Vec<_dim, T, policy> operator^(const Vec<_dim, T, policy> &p)const;
	  //const Vec<_dim, T, policy> operator*(T a)const ;
	  Vec<_dim, T, policy> &operator=(T a);
	  Vec<_dim, T, policy> &operator*=(T a);
	  Vec<_dim, T, policy> &operator/=(T a);
	  Vec<_dim, T, policy> &operator+=(T a);
	  Vec<_dim, T, policy> &operator-=(T a);
	  Vec<_dim, T, policy> &operator+=(const Vec<_dim, T, policy> &p);
	  Vec<_dim, T, policy> &operator-=(const Vec<_dim, T, policy> &p);
	  Vec<_dim, T, policy> operator*(T a)const ;
	  Vec<_dim, T, policy> operator/(T a)const ;
	  Vec<_dim, T, policy> operator+(const Vec<_dim, T, policy> &p)const ;
	  Vec<_dim, T, policy> operator-(const Vec<_dim, T, policy> &p)const ;
	  T abs2()const;
	  T abs()const;
	  T distance(const Vec<_dim, T, policy> p)const;
	  Vec<_dim, T, policy> & normalize();
	  Vec<_dim, T, policy> & normalized()const;

	  template<indexType _dim2, class U, class p >
	  friend std::ostream & operator<< (std::ostream &out, const Vec<_dim2, U, p> &v);
	  template<indexType _dim2, class U, class p>
	  friend std::istream & operator>>(std::istream &in, Vec<_dim2, U, p> &v);
	static const unsigned int dim;
	private:
	  T x[_dim];
	};
	
	template<indexType _dim, class T, class policy>
	const unsigned int Vec<_dim, T, policy>::dim=_dim;



	template<indexType _dim, class T, class policy>
	inline
	  std::ostream & operator<< (std::ostream &out, const Vec<_dim, T, policy> &v){
		out<<v.x[0];
		for(int i=1; i<v.dim; i++){
			out<<"  "<<v.x[i];
			}
		return out;
		}

	template<indexType _dim, class T, class policy>
	  std::istream & operator>>(std::istream &in, Vec<_dim, T, policy> &v){
		for(int i=0; i<v.dim; i++){
			in>>v.x[i];
			}
		return in;
		}

	template<indexType _dim, class T, class policy>
	Vec<_dim, T, policy> Vec<_dim, T, policy>::operator^(const Vec<_dim, T, policy> &p)const{
	//	return Vec<_dim, T, policy>(p.x);
		}

	template<indexType _dim, class T, class policy>
	inline
	  T Vec<_dim, T, policy>::operator*(const Vec<_dim, T, policy> &p)const{
		T prod=(T)0;
		for(int i=0; i<dim; i++){
			prod+=policy::mul(x[i],p.x[i]);
			}
		return prod;
		}

	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> &Vec<_dim, T, policy>::operator=(T a){
		for(int i=0; i<dim; i++){
			x[i]=a; 
			}
		return *this;
		}

	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> &Vec<_dim, T, policy>::operator*=(T a){
		for(int i=0; i<dim; i++){
			x[i]*=a; 
			}
		return *this;
		}
	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> &Vec<_dim, T, policy>::operator/=(T a){
		assert(a!=0);
		for(int i=0; i<dim; i++){
			x[i]/=a; 
			}
		return *this;
		}
	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> &Vec<_dim, T, policy>::operator+=(T a){
		for(int i=0; i<dim; i++){
			x[i]+=a; 
			}
		return *this;
		}
	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> &Vec<_dim, T, policy>::operator-=(T a){
		for(int i=0; i<dim; i++){
			x[i]-=a; 
			}
		return *this;
		}

	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> &Vec<_dim, T, policy>::operator+=(const Vec<_dim, T, policy> &p){
		for(int i=0; i<dim; i++){
			x[i]+=p.x[i];
			}
		return *this;
		}

	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> &Vec<_dim, T, policy>::operator-=(const Vec<_dim, T, policy> &p){
		for(int i=0; i<dim; i++){
			x[i]-=p.x[i];
			}
		return *this;
		}

	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> Vec<_dim, T, policy>::operator*(T a)const {
		Vec<_dim, T, policy> pp(*this);
			pp*=(a);
		return pp;
		}


	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> Vec<_dim, T, policy>::operator/(T a)const {
		Vec<_dim, T, policy> pp(*this);
			pp/=(a);
		return pp;
		}

	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> Vec<_dim, T, policy>::operator+(const Vec<_dim, T, policy> &p)const {
		Vec<_dim, T, policy> pp(*this);
		pp+=p;
		return pp;
		}


	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> Vec<_dim, T, policy>::operator-(const Vec<_dim, T, policy> &p)const {
		Vec<_dim, T, policy> pp(*this);
		pp-=p;
		return pp;
		}

	template<indexType _dim, class T, class policy>
	inline
	  T Vec<_dim, T, policy>::abs2()const{
		return (*this)*(*this);
		}

	template<indexType _dim, class T, class policy>
	inline
	  T Vec<_dim, T, policy>::abs()const{
		return sqrt(abs2());
		}

	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> & Vec<_dim, T, policy>::normalize(){
		T d=abs();
		assert(d!=0);
		d=1/d;
		
		(*this)*=d;
		return *this;
		}
	template<indexType _dim, class T, class policy>
	inline
	  Vec<_dim, T, policy> & Vec<_dim, T, policy>::normalized()const{
		return Vec<_dim, T, policy> (*this).normalize();
		}

	template<indexType _dim, class T, class policy>
	inline
	  const Vec<_dim, T, policy> operator *(double a, Vec<_dim, T, policy> v){
		v*=(a);
		return v;
		};

	template<indexType _dim, class T, class policy>
	inline
	  const Vec<_dim, T, policy> operator -(Vec<_dim, T, policy> v){
		v*=(T)(-1.0);
		return v;
		};

	//template<>
	//Vec::Vec(double , double, double){}

	template<indexType _dim, class T, class policy>
	inline
	Vec<_dim, T, policy> cross(const Vec<_dim, T, policy> u, const Vec<_dim, T, policy> v){
		//only first 3 components
		assert(u.dim>=3 and v.dim>=3);
		return Vec<_dim, T, policy> (policy::mul(u(1),v(2))-policy::mul(u(2),v(1)),policy::mul(u(2),v(0))-policy::mul(u(0),v(2)),policy::mul(u(0),v(1))-policy::mul(u(1),v(0)));
	}

#endif

