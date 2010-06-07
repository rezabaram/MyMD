#ifndef ELLIPSOID_H
#define ELLIPSOID_H 

#include"matrix.h"
#include"geombase.h"

////FIXME added interaction between classes. this is not so good. make a better struture
#include"eigen.h"
#include"polynom.h"
#include"ray.h"
typedef GeomObject<tellipsoid> CEllipsoid;
CQuadratic intersect (const CRay<HomVec> &ray, const CEllipsoid &E);
////

using namespace math;
typedef matrix<double> Matrix;

template<typename T>
inline
T tmax(const T &a, const T &b){
	if(a>b)return a;
	return b;
	}

typedef matrix<double> Matrix;
void quaternionToMatrix(const Quaternion &q, Matrix &M){//got from wikipedia
TRY
	double Nq = q.abs2();
	static double s;
	if( Nq > 0.0) s = 2.0/Nq; else s = 0.0;
	double X = q.v(0)*s,   Y = q.v(1)*s,  Z = q.v(2)*s;

	double wX = q.u*X,    wY = q.u*Y,    wZ = q.u*Z;
	double xX = q.v(0)*X, xY = q.v(0)*Y, xZ = q.v(0)*Z;
	double yY = q.v(1)*Y, yZ = q.v(1)*Z, zZ = q.v(2)*Z;

	M(0,0)=1.0-(yY+zZ); M(1,0)=xY-wZ ;       M(2,0)= xZ+wY;
	M(0,1)= xY+wZ;      M(1,1)=1.0-(xX+zZ);  M(2,1)=yZ-wX;
	M(0,2)=xZ-wY;       M(1,2)= yZ+wX;       M(2,2)=1.0-(xX+yY);
	return;
CATCH
}


template <size_t dim>
inline
Vec<dim> operator *(const Matrix &M, const Vec<dim> &v){
TRY
	assert(M.RowNo()>=v.dim);
	assert(M.IsSquare());
	Vec<dim> u;
	
	for(size_t i=0; i<v.dim; i++){
		u(i)=0;
		for(size_t j=0; j<v.dim; j++)
			u(i)+=v(j)*M(i,j);
		}
	return u;
CATCH
	}

template <size_t dim>
inline
Vec<dim> operator *(const Vec<dim> &v, const Matrix &M){
TRY
	assert(M.RowNo()>=v.dim);
	assert(M.IsSquare());
	Vec<dim> u;
	
	for(size_t i=0; i<v.dim; i++){
		u(i)=0;
		for(size_t j=0; j<v.dim; j++)
			u(i)+=v(j)*M(j,i);
		}
	return u;
CATCH
	}

template <class T>
Matrix & operator +=(Matrix &M, const double d){
TRY
	assert(M.RowNo() == M.ColNo());
	for(indexType i=0; i<M.RowNo(); ++i){
		M(i,i)+=d;
		}
	return M;
CATCH
	}


template<>
class GeomObject<tellipsoid>: public GeomObjectBase{
	public:	
		
	GeomObject(const vec &v,double _a, double _b, double _c, double _r=1, const Quaternion &_q=Quaternion(1,0,0,0)):GeomObjectBase(v,tellipsoid, _q), a(_a), b(_b), c(_c), R(_a,_b,_c), P(HomVec(0,0,0,1)) {
		identifier=14;
		radius=tmax(a, tmax(b,c));
		a*=_r;
		b*=_r;
		c*=_r;
		radius*=_r;
		setup();

		//put a point on the surface of the ellipse	
		
		CRay<HomVec> ray(HomVec(0,0,0, 1) , HomVec(0,0,1, 1) );
		CQuadratic quartic(intersect(ray, *this));
		//the roots are sorted ascending
		//P0= ray(quartic.root(0).real()); //on the surface of E1
		//quartic.print_roots(cerr);
		//ray2.print(cerr);
		//cerr<< P0 <<endl;
		rotateTo(q);
		update_tranlation_mat();
		}

	~GeomObject(){}
	HomVec toWorld(const HomVec &point)const{
		return (!(rotat_mat*trans_mat))*point;
		}
	HomVec toBody(const HomVec &point)const{
		return (rotat_mat*trans_mat)*point;
		}
	void fixToBody(const HomVec &point){
		P=point;
		P0=(rotat_mat*trans_mat)*point;
		}

	void setup(){
		//Elements of the rotational matrix

		double temp=0;
		for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j){
		  if(i==j)temp=1.0;
		  else temp=0;
		  
		  rotat_mat(i,j)= temp;
		  scale_mat(i,j)=temp;
		  inv_scale_mat(i,j)=temp;
		  inert_mat(i,j)= temp;
		  trans_mat(i,j)= temp;
		}

		  //Elements of the scaling matrix
		  scale_mat(0,0)=1.0/(a*a);
		  scale_mat(1,1)=1.0/(b*b);
		  scale_mat(2,2)=1.0/(c*c);
		  scale_mat(3,3)=-1;//in homogeneous formulation
		  ellip_mat=scale_mat;

		  inv_scale_mat(0,0)=(a*a);
		  inv_scale_mat(1,1)=(b*b);
		  inv_scale_mat(2,2)=(c*c);

		  inv_scale_vec(0)=(a*a);
		  inv_scale_vec(1)=(b*b);
		  inv_scale_vec(2)=(c*c);

		  //angular moment of inertial
		  inert_mat(0,0)=0.2*(b*b+c*c);
		  inert_mat(1,1)=0.2*(a*a+c*c);
		  inert_mat(2,2)=0.2*(a*a+b*b);



		}
	vec gradient (const vec &X)const {
		return 2.0*ellip_mat*(X-Xc);
		}
	double operator() (const HomVec &X)const {
		return X*ellip_mat*X;
		}

	Matrix inv()const{
		return (~rotat_mat*inv_scale_mat*rotat_mat); 
		}

	void update_tranlation_mat(){
	TRY	
		trans_mat(0,3)=-Xc(0);
		trans_mat(1,3)=-Xc(1);
		trans_mat(2,3)=-Xc(2);
		trans_mat(3,3)=1;
		ellip_mat=~trans_mat*~rotat_mat*scale_mat*rotat_mat*trans_mat;
		//P=HomVec(0.1,0.1,0.1,1);
		P=(!(rotat_mat*trans_mat))*P0;
	CATCH
		}

	void moveto(const vec &v){
	TRY
		Xc=v;
		update_tranlation_mat();
	CATCH
		}

	void rotateTo(const Quaternion &q){
	TRY
		quaternionToMatrix(q, rotat_mat);
		//Matrix temp=(rotat_mat*(~rotat_mat));
		//cerr<< temp <<endl;
		//assert(fabs((rotat_mat*(~rotat_mat)).Det-1) < 0.0001);
		//ellip_mat=~rotat_mat*scale_mat*rotat_mat;
		update_tranlation_mat();
		//P=ellip_mat*P0;
	CATCH
		}

	double I(vec n){//FIXME only in special coordinate system
	TRY
		n.normalize();
		return n*inert_mat*n;
	CATCH
		}

	double vol(){
		return 4.0/3.0*M_PI*a*b*c;
		}

	void rotate(const vec& n , double alpha){//FIXME
		//q.setRotation(n, alpha);
		//orientation=q.rotate(orientation);
		}

	void shift(const vec& v){
		Xc+=v;
		update_tranlation_mat();
		}

	void scale(double scale){
		a*=scale;
		b*=scale;
		c*=scale;
		radius*=scale;
		setup();
		}

	vec point_to_plane(const CPlane &plane)const{//FIXME need to be obtimized
	TRY
		double alpha;
		alpha=(plane.n*(this->inv())*plane.n);
	
		ERROR(alpha<0, "Impossible happened");

		alpha=1/sqrt(alpha);
		static vec m(0.0);
		m=alpha*(~rotat_mat*(inv_scale_vec^(rotat_mat*plane.n))); //this more efficient form of m=(alpha*(!ellip_mat)*plane.n);
		double d1=plane.normal_from_point(Xc+m).abs();
		double d2=plane.normal_from_point(Xc-m).abs();

		if(d1<d2)return Xc+m;
		else return Xc-m;
	CATCH
		}

	void print(std::ostream &out)const{
		//FIXME temporary 
		CSphere S(P.project(), 0.01);
		S.print(out); 
	//	return;
		out<<endl;
		///

		out<< identifier<< "   ";
		out<< Xc<< "  "<<radius+0.001<<"  ";
		out<< ellip_mat(0,0) << "  " <<ellip_mat(1,1)<< "  "<<ellip_mat(2,2)<< "  ";
		out<< ellip_mat(1,0) << "  " <<ellip_mat(1,2)<< "  "<<ellip_mat(0,2)<< "  ";
		out<< 0 << "  " << 0 <<  "  " <<0<< "  ";
		out<<-1;
		return ;
		}
	
	void parse(std::istream &in){//FIXME
			ERROR(true,"not implemented");
			//in>>identifier;
			}

	Matrix rotat_mat;
	Matrix trans_mat;
	Matrix scale_mat;
	Matrix inv_scale_mat;
	vec    inv_scale_vec;

	Matrix ellip_mat;
	Matrix inert_mat;

	double a,b,c;
	vec R;
	HomVec P, P0;
	private:
	GeomObject<tellipsoid> (const GeomObject<tcomposite> & p);//not allow copies
	GeomObject<tellipsoid> ();
	};

////FIXME added interaction between classes. this is not so good. make a better struture
CQuadratic intersect (const CRay<HomVec> &ray, const CEllipsoid &E){
	double a=ray.n*E.ellip_mat*ray.n;
	double b=ray(0)*E.ellip_mat*ray.n+ray.n*E.ellip_mat*ray(0);
	double c=ray(0)*E.ellip_mat*ray(0);
	return CQuadratic(a, b, c);
}
//----------------
#endif /* ELLIPSOID_H */
