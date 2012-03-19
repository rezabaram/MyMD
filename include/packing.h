// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef PACKING_H
#define PACKING_H 
#include"exception.h"
#include"common.h"
#include"shapes.h"
#include"grid.h"
#include"slice.h"
#include"ellips_contact.h"
#include"MersenneTwister.h"
#include"interaction.h"
#include"contact_network.h"
#include"material.h"


//member classes for contact network

template < class T>
class CPacking : public list<T *>
	{
	//void push_back ( const T* x );
	public:
	CPacking():contacts(ContactNetwork<T>(this)), TotalParticlesN(0)
		{
		maxr=0.0;
		grid_built=false;
		}
	~CPacking(){
		typename list<T *>::iterator it;
		if(grid_built)delete grid;
		for(it=this->begin(); it!=this->end(); it++)
			delete (*it);
		}

	void add( T *x){
		push_back (x);
		};

	void print(std::ostream& out, bool raster=false)const;
	void print_particle_axes(ostream &out,const vec3d &c1=vec3d(0,0,0),const vec3d &c2=vec3d(1,1,1) );
	void printRaster3D(std::ostream& out)const{print(out, true);}
	void printEuler(std::ostream& out)const;
	void parse(string infilename, bool periodic=false);
	void parse(istream &inputFile, bool periodic=false);
	void output_contact_network(ostream &out);
	void BuildContactNetwork();
	void BuildGrid(){
		if(grid_built){
			//WARNING("The grid is already built.");
			return;
			}
		grid=(new CRecGrid<T>(this, vec(0.0, 0.0, 0.0), vec(1.5, 1.5, 1.5), 1.1*maxr));
		grid->build();
		grid_built=true;
		};
	double packFraction(const vec &x1, const vec &x2, unsigned long N=10000, double scale=1);
	double get_slice_uniform(CSlice &slice,  double z=0.5, unsigned long N=100);
	double totalVolume();
	double get_slice(CSlice &slice,  double z=0.5, unsigned long N=10000);
	bool is_in_void(const vec &x);

	double avg_contact_number(){
		return contacts.avg_contact_number();
		}
	ContactNetwork<T> contacts;
	
	CRecGrid<T> * grid;
	long TotalParticlesN;
	double maxr;
 	private:
	bool grid_built;
	};

template < class T>
double CPacking<T>::totalVolume(){
	double v=0;
	typename CPacking<T>::const_iterator it;
	for(it=this->begin(); it!=this->end(); it++){
		if(! (*it)->is_shadow) v+=(**it).vol();
		}
	return v;
	}

//calculating packing fraction inside cube with x1, x2 az opposite corners
template < class T>
double CPacking<T>::packFraction(const vec &x1, const vec &x2, unsigned long N, double scale){
TRY
	BuildGrid();
	double n=0, nt=0;
	vec x;
	//vec l=x2-x1;
	//double boxvolume=fabs(l(0)*l(1)*l(2));
	double d=1;
	typename CNode3D<T>::const_iterator it;

	typename CPacking<T>::iterator it2;
	if(fabs(scale-1)>1e-10) 
		for(it2=this->begin(); it2!=this->end(); it2++){
			(*it2)->shape->scale(scale);
		}

	for(unsigned long i=0; i<N; i++){
		++nt;
		x=randomVec(x1,x2);
		CNode3D<T> *node=grid->which(x);
		assert(node);
		for(it=node->begin(); it!=node->end(); it++){
			d=(*(**it).shape)(x);
			if(d<0) {
				++n;
				break;
				}
			}
		}
	return n/nt;
CATCH
	}

template < class T>
double CPacking<T>::get_slice_uniform(CSlice &slice,  double z, unsigned long N){
	vec2d x1=slice.corner;
	vec2d x2=slice.corner+slice.L;
	
	BuildGrid();
	double d;
	double n=0, nt=0;
	double dx=1.0/(double)N;
	vec2d x=x1;
	for(;;x+=vec2d(dx,0)){
		if(x(0)>x2(0)){
			x(0)=x1(0);
			x(1)+=dx;
			if(x(1)>x2(1))break;
		}
	
		++nt;
		bool in_void=true;
		vec x3d=vec(x(0),x(1),z);
		CNode3D<T> *node=grid->which(x3d);
		assert(node);
		typename CNode3D<T>::const_iterator it;
		for(it=node->begin(); it!=node->end(); it++){
			d=(*(**it).shape)(x3d);
			if(d<0) {
				++n;
				in_void=false;
				break;
				}
			}
		if(in_void)slice.add(x, 0);
		else slice.add(x, 1);
		}

	return n/nt;
}
template < class T>
double CPacking<T>::get_slice(CSlice &slice,  double z, unsigned long N){
	vec2d x1=slice.corner;
	vec2d x2=slice.corner+slice.L;
	
	BuildGrid();
	double d;
	double n=0, nt=0;
	for(unsigned long i=0; i<N; i++){
		++nt;
		bool in_void=true;
		vec2d x=randomVec2d(x1,x2);
		vec x3d=vec(x(0),x(1),z);
		CNode3D<T> *node=grid->which(x3d);
		assert(node);
		typename CNode3D<T>::const_iterator it;
		for(it=node->begin(); it!=node->end(); it++){
			d=(*(**it).shape)(x3d);
			if(d<0) {
				++n;
				in_void=false;
				break;
				}
			}
		if(in_void)slice.add(x, 0);
		else slice.add(x, 1);
		}

	return n/nt;
}

template < class T>
bool CPacking<T>::is_in_void(const vec &x){
TRY
	double d=1;
	typename CNode3D<T>::const_iterator it;

		CNode3D<T> *node=grid->which(x);
		assert(node);
		for(it=node->begin(); it!=node->end(); it++){
			d=(*(**it).shape)(x);
			if(d<=0) {
				return false;
				}
			}
	return true;
CATCH
	}

template < class T>
void CPacking<T>::print(std::ostream& out, bool raster)const{
		typename CPacking::const_iterator it;
		for(it=this->begin(); it!=this->end(); it++){
			if(raster) (*it)->shape->printRaster3D(out);
			else (*it)->shape->print(out);
			out<<endl;
			}
		}

template < class T>
void CPacking<T>::printEuler(std::ostream& out)const{
		typename CPacking::const_iterator it;
		for(it=this->begin(); it!=this->end(); it++){
			(*it)->shape->print_in_euler(out);
			out<<endl;
			}
		}

template < class T>
void CPacking<T>::parse(string infilename, bool periodicity) {
	ifstream inputFile(infilename.c_str());
	ERROR(!inputFile.good(), "Unable to open input file: ")

	parse(inputFile, periodicity);
	inputFile.close();
	}

template < class T>
void CPacking<T>::parse(istream &inputFile, bool periodic) {

	string line;

	//Parse the line
	while(getline(inputFile,line))
	{

	//Insert the line string into a stream
	stringstream ss(line);

	int id;
	ss>>id;
	if(id==14){
		CEllipsoid shape;
		shape.parse(ss);
		T *p=new T(shape);
		//p->shape->parse(ss);
		this->push_back(p);
		if(p->shape->radius>maxr)maxr=p->shape->radius;

		//when periodic
		if(periodic){
		if(shape.Xc(0)-shape.radius<0){
			p=new T(shape, true);
			p->shift(vec(1,0,0));
			this->push_back(p);
			}
		else if(shape.Xc(0)+shape.radius>1){
			p=new T(shape, true);
			p->shift(vec(-1,0,0));
			this->push_back(p);
			}
		if(shape.Xc(1)-shape.radius<0){
			p=new T(shape, true);
			p->shift(vec(0,1,0));
			this->push_back(p);
			}
		else if(shape.Xc(1)+shape.radius>1){
			p=new T(shape, true);
			p->shift(vec(0,-1,0));
			this->push_back(p);
			}
		//corner
		if(shape.Xc(0)-shape.radius<0 and shape.Xc(1)-shape.radius<0){
			p=new T(shape, true);
			p->shift(vec(1,1,0));
			this->push_back(p);
			}
		else if(shape.Xc(0)-shape.radius<0 and shape.Xc(1)+shape.radius>1){
			p=new T(shape, true);
			p->shift(vec(1,-1,0));
			this->push_back(p);
			}
		else if(shape.Xc(0)+shape.radius>1 and shape.Xc(1)-shape.radius<0){
			p=new T(shape, true);
			p->shift(vec(-1,1,0));
			this->push_back(p);
			}
		else if(shape.Xc(0)+shape.radius>1 and shape.Xc(1)+shape.radius>1){
			p=new T(shape, true);
			p->shift(vec(-1,-1,0));
			this->push_back(p);
			}
		}
		}

	else{
		//WARNING("Skip reading shape with id: "<< id);
		continue;
		}
		
	}
}
template < class T>
void CPacking<T>::BuildContactNetwork() {
	contacts.build();
	}

template<class T>
void CPacking<T>::output_contact_network(std::ostream &out)
{
TRY
	contacts.print(out);
CATCH
}


template<class T>
void CPacking<T>::print_particle_axes(ostream &out,const vec3d &c1,const vec3d &c2){
	typename CPacking<T>::const_iterator it1;
	for(it1=this->begin();  it1!=this->end(); it1++){
			vec3d x=(*it1)->shape->Xc;
			if( x(0)<c1(0) or 
			x(0)>c2(0) or 
			x(1)<c1(1) or 
			x(1)>c2(1) or 
			x(2)<c1(2) or 
			x(2)>c2(2) ) continue;
		     (*it1)->shape->print_coord_sys(out);
		}
	return;
	}

#endif /* PACKING_h */
