#ifndef ELLIPSOID_H
#define ELLIPSOID_H 

#include"matrix.h"
#include"geombase.h"

////FIXME added interaction between classes. this is not so good. make a better struture
#include"polynom.h"
#include"ray.h"
#include"plane.h"

class CEllipsoid;
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
		//for(size_t j=v.dim; j<M.RowNo(); j++)
			//u(i)+=M(i,j);
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
		//for(size_t j=v.dim; j<M.RowNo(); j++)
			//u(i)+=M(j,i);
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


class CEllipsoid: public GeomObjectBase
	{
	public:	
		
	CEllipsoid():GeomObjectBase(vec(0,0,0),tellipsoid, Quaternion(1,0,0,0)), a(1), b(1), c(1){
		identifier=14;
		}

	CEllipsoid(const vec &v,double _a, double _b, double _c, const Quaternion &_q=Quaternion(1,0,0,0)):GeomObjectBase(v,tellipsoid, _q), a(_a), b(_b), c(_c)
		{
		identifier=14;
		radius=tmax(a, tmax(b,c));
		setup();

		rotateTo(q);
		update_tranlation_mat();
		}

	~CEllipsoid(){}
	HomVec toWorld(const HomVec &point)const
		{
		return (!(rotat_mat*trans_mat))*point;
		}

	HomVec toBody(const HomVec &point)const
		{
		return (rotat_mat*trans_mat)*point;
		}
	void fixToBody(const HomVec &point)
		{
		P=point;
		P0=(rotat_mat*trans_mat)*point;
		}

	virtual GeomObjectBase *clone(){
		return new CEllipsoid(*this);
		}

	void mat_init()
		{
		double temp=0;
		for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			{
			if(i==j)temp=1.0;
			else temp=0;
			  
			rotat_mat(i,j)= temp;
			scale_mat(i,j)=temp;
			inv_scale_mat(i,j)=temp;
			inert_mat(i,j)= temp;
			trans_mat(i,j)= temp;
			}
		}

	void setup(){

		mat_init();

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

	vec gradient (const vec &X)const 
		{
		return 2.0*ellip_mat*(X-Xc);
		}

	HomVec gradient (const HomVec &X)const 
		{
		return 2.0*ellip_mat*X;
		}

	double operator() (const vec &X)const {
		return (X-Xc)*ellip_mat*(X-Xc)-1;
		}
	double operator() (const HomVec &X)const {
		return X*ellip_mat*X;
		}

	Matrix inv()const{
		return (~rotat_mat*inv_scale_mat*rotat_mat); 
		}

	void update_tranlation_mat() {
	TRY	
		trans_mat(0,3)=-Xc(0);
		trans_mat(1,3)=-Xc(1);
		trans_mat(2,3)=-Xc(2);
		trans_mat(3,3)=1;
		static Matrix tempmat;
		tempmat=rotat_mat*trans_mat;
		ellip_mat=(~tempmat)*scale_mat*tempmat;
		//P=HomVec(0.1,0.1,0.1,1);
		//P=(!(rotat_mat*trans_mat))*P0;
	CATCH
		}

	void moveto(const vec &v){
	TRY
		Xc=v;
		update_tranlation_mat();
	CATCH
		}

	void rotateTo(const Quaternion &q) {
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

	double vol()const{
		return 4.0/3.0*M_PI*a*b*c;
		}

	void rotate(const vec& n, double alpha){//FIXME
		//q.setRotation(n, alpha);
		//orientation=q.rotate(orientation);
		}

	void shift(const vec& v){
		Xc+=v;
		update_tranlation_mat();
		}

	void scale(double scale)
		{
		a*=scale;
		b*=scale;
		c*=scale;
		radius*=scale;
		setup();
		}

	bool doesHit(const CPlane &plane)const {//FIXME needs to be obtimized
	TRY
		if(fabs(plane(this->Xc)) > this->radius) return false;
		double alpha;
		alpha=(plane.n*(this->inv())*plane.n);
	
		ERROR(alpha<=0, "Impossible happened");

		alpha=1/sqrt(alpha);
		static vec m(0.0);
		m=alpha*(~rotat_mat*(inv_scale_vec^(rotat_mat*plane.n))); //this is more efficient form of m=(alpha*(!ellip_mat)*plane.n);

		if(((Xc+m-plane.Xc)*plane.n)*((Xc-m-plane.Xc)*plane.n) <= 0)return true;//extreme points lie on both side of the plane
		else return false;
	CATCH
		}

	vec point_to_plane(const CPlane &plane)const{//FIXME needs to be obtimized
	TRY
		double alpha;
		alpha=(plane.n*(this->inv())*plane.n);
	
		ERROR(alpha<0, "Impossible happened");

		alpha=1/sqrt(alpha);
		static vec m(0.0);
		m=alpha*(~rotat_mat*(inv_scale_vec^(rotat_mat*plane.n))); //this is more efficient form of m=(alpha*(!ellip_mat)*plane.n);
		double d1=plane.normal_from_point(Xc+m).abs();
		double d2=plane.normal_from_point(Xc-m).abs();

		if(d1<d2)return Xc+m;
		else return Xc-m;
	CATCH
		}


	virtual const void print_coord_sys(ostream &out){
		vec3d n1, n2, n3;
		if(a>=b and b>=c){
			n1=q.toWorld(vec3d(1,0,0));
			n2=q.toWorld(vec3d(0,1,0));
			n3=q.toWorld(vec3d(0,0,1));
			}
		else if(a>=c and c>=b){
			n1=q.toWorld(vec3d(1,0,0));
			n2=q.toWorld(vec3d(0,0,1));
			n3=q.toWorld(vec3d(0,1,0));
			}
		else if(b>=c and c>=a){
			n1=q.toWorld(vec3d(0,1,0));
			n2=q.toWorld(vec3d(0,0,1));
			n3=q.toWorld(vec3d(1,0,0));
			}
		else if(b>=a and a>=c){
			n1=q.toWorld(vec3d(0,1,0));
			n2=q.toWorld(vec3d(1,0,0));
			n3=q.toWorld(vec3d(0,0,1));
			}
		else if(c>=a and a>=b){
			n1=q.toWorld(vec3d(0,0,1));
			n2=q.toWorld(vec3d(1,0,0));
			n3=q.toWorld(vec3d(0,1,0));
			}
		else if(c>=b and b>=a){
			n1=q.toWorld(vec3d(0,0,1));
			n2=q.toWorld(vec3d(0,1,0));
			n3=q.toWorld(vec3d(1,0,0));
			}
		else{
			ERROR(1, "Non-handled case.");
			}

		//n1(0)=1-2.0*drand48(); n1(1)=1-2.0*drand48(); n1(2)=1-2.0*drand48();
		//n2(0)=1-2.0*drand48(); n2(1)=1-2.0*drand48(); n2(2)=1-2.0*drand48();
		//n3(0)=1-2.0*drand48(); n3(1)=1-2.0*drand48(); n3(2)=1-2.0*drand48();
/*
		double theta=M_PI*drand48();	
		double phi=2*M_PI*drand48();	
		n1=vec3d(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)); 
		theta=M_PI*drand48();	
		phi=2*M_PI*drand48();	
		n2=vec3d(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)); 
		theta=M_PI*drand48();	
		phi=M_PI*(1-2*drand48());	
		n3=vec3d(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)); 
		//n1.normalize();
		//n2.normalize();
		//n3.normalize();
		//out<< n1 <<"\t"<< n2 <<"\t"<<n3 <<endl;
*/
		out<< spherical(n1) <<"\t"<< spherical(n2) <<"\t"<<spherical(n3) <<endl;
		return;
		}
	void print(std::ostream &out)const{
		out<<setprecision(12)<< identifier<< "   ";
		out<< Xc<< "  "<< a<<"  "<<b<<"  "<< c <<"  ";
		out<< q;
		}

	void print_in_euler(std::ostream &out)const{
		out<<setprecision(12)<< identifier<< "   ";
		out<< Xc<< "  "<< a<<"  "<<b<<"  "<< c <<"  ";
		out<< this->euler();
		}

	void printRaster3D(std::ostream &out)const{
		//FIXME temporary 
//		CSphere S(P.project(), 0.01);
//		S.print(out); 
	//	return;
//		out<<endl;
		///
		out<<setprecision(12)<< identifier<< "   ";
		out<< Xc<< "  "<<radius+0.00001<<"  ";
		out<< ellip_mat(0,0) << "  " <<ellip_mat(1,1)<< "  "<<ellip_mat(2,2)<< "  ";
		out<< ellip_mat(1,0) << "  " <<ellip_mat(1,2)<< "  "<<ellip_mat(0,2)<< "  ";
		out<< 0 << "  " <<  0 <<  "  " << 0 << "  ";
		out<<-1;
		//out<< ellip_mat(0,3) << "  " << ellip_mat(1,3) <<  "  " <<ellip_mat(2,3)<< "  ";
		//out<<ellip_mat(3,3);
		//cerr<< trans_mat<<endl;
		return ;
		}
	
	void parse(std::istream &in){//FIXME
		in>>Xc;
		in>> a >> b >> c;
		q.parse(in);

		radius=tmax(a, tmax(b,c));
		setup();

		rotateTo(q);
		update_tranlation_mat();

		}

	friend std::ostream & operator<< (std::ostream &out, const CEllipsoid &E);
	Matrix rotat_mat;
	Matrix trans_mat;
	Matrix scale_mat;
	Matrix inv_scale_mat;
	vec    inv_scale_vec;

	Matrix ellip_mat;
	Matrix inert_mat;

	double a,b,c;
	
	private:
	//CEllipsoid (const CEllipsoid & p);//not allow copies

	};
  
std::ostream & operator<< (std::ostream &out, const CEllipsoid &E){
		E.print(out);
		return out;
		};

////FIXME added interaction between classes. this is not so good. make a better struture
CQuadratic intersect (const CRay<HomVec> &ray, const CEllipsoid &E){
	double a=ray.n*E.ellip_mat*ray.n;
	double b=ray(0)*E.ellip_mat*ray.n+ray.n*E.ellip_mat*ray(0);
	double c=ray(0)*E.ellip_mat*ray(0);
	return CQuadratic(a, b, c);
}
bool intersect (size_t i, HomVec &X, const CRay<HomVec> &ray, const CEllipsoid &E){
	CQuadratic q=intersect(ray, E);
	if(q.root(0).imag()>0.0000001) return false;
	X=ray(q.root(i).real());
	return true;
}
//----------------
#endif /* ELLIPSOID_H */
