#ifndef CELLLIST_H
#define CELLLIST_H 
#include<list>
#include"exception.h"
#include"interaction_force.h"
#include"box.h"
using namespace std;


template<typename TParticle>
class CCell : public list<TParticle *>{


	public:
	CCell(){
		//adding itself to the list 
		//of neighbors
		neighs[0]=this;
		shifts[0]=0.0;
		}
	void interact(){
	
		typename CCell<TParticle>::iterator it1, it2;
		for(it1=this->begin(); it1!=this->end(); it1++){
			for(it2=this->begin(); it2!=it1; it2++){
				//for the same cell
				Test::interact(*it1, *it2);
				}

			for(int k=1; k<14; k++){
				//for neighbouring cells
				if(neighs[k]==NULL)continue;
				if(shifted[k])(**it1).shift(shifts[k]);
				for(it2=neighs[k]->begin(); it2!=neighs[k]->end(); it2++){
					Test::interact(*it1, *it2);
					}
				if(shifted[k])(**it1).shift(-shifts[k]);
				}
			}
		}

	void add(TParticle *p){
		push_back(p);
		}
	template<typename T, typename U>
	friend 
	class CCellList;

	private:
	CCell<TParticle> *neighs[14];
	vec shifts[14];
	bool shifted[14];
	};


template<typename TParticleContainer, typename TParticle>
class CCellList
	{
	public:
	CCellList(CBox *box):nodes(NULL),periodic_x(false),periodic_y(false), periodic_z(false), is_set_up(false) {
		c=box->corner;
		diag=box->L;
		if(box->btype=="periodic_x")periodic_x=true;
		else if(box->btype=="periodic_xy"){
			periodic_x=true;
			periodic_y=true;
			}
		else if(box->btype=="periodic_xyz"){
			periodic_x=true;
			periodic_y=true;
			periodic_z=true;
			}
		else{
			ERROR(1,"Boundary type not defined")
			}

		}

	~CCellList(){
		delete [] nodes;
		}
	void setup(double d){
		if(is_set_up)return;
		
		ERROR(d<1e-10, "Invalid grid size: "+stringify(d));
		cerr<< "Constructing the grid ... ";
		nx=max(1,(int)floor(diag(0)/d));
		ny=max(1,(int)floor(diag(1)/d));
		nz=max(1,(int)floor(diag(2)/d));
		dx=diag(0)/(double)nx;
		dy=diag(1)/(double)ny;
		dz=diag(2)/(double)nz;
		if(nodes!=NULL)delete [] nodes;
		cerr<< nx*ny*nz <<endl;
		nodes=new CCell<TParticle>[nx*ny*nz];
		build_neighbors();
		cerr<< "done: "<<nx<<" X "<<ny <<" X "<<nz<<endl;
		is_set_up=true;
		}
	void clear(){
		for(int i=0; i<nx*ny*nz; i++){
		nodes[i].clear();
		}
		
		}

	void build(TParticleContainer &p){
		clear();
		typename TParticleContainer::iterator it;
		for(it=p.begin(); it!=p.end(); it++){
			add((*it));
			}
		}
	void update(TParticleContainer &p){
		typename TParticleContainer::iterator it;
		for(it=p.begin(); it!=p.end(); it++){
			update((*it));
			}
		}
	CCell<TParticle> *which(TParticle &p){
		vec &x=p.shape->Xc;
		vec xp=x-c;
		int i=(int)(floor(xp(0)/dx));
		int j=(int)(floor(xp(1)/dy));
		int k=(int)(floor(xp(2)/dz));
		if(!periodic_x) ERROR(i<0 or i>=nx, "Point out of grid: "+stringify(x)+stringify(i)+" " +stringify(j)+" "+stringify(k));
		if(!periodic_y) ERROR(j<0 or j>=ny, "Point out of grid: "+stringify(x)+stringify(i)+" " +stringify(j)+" "+stringify(k));
		if(!periodic_z) ERROR(k<0 or k>=nz, "Point out of grid: "+stringify(x)+stringify(i)+" " +stringify(j)+" "+stringify(k));
		vec shift(0,0,0);
		bool is_shifted=false;
		CCell<TParticle> *cell=boundary_mask(i,j, k, shift, is_shifted);
		if(is_shifted)p.shift(shift);
		return cell;
		}

	void add(TParticle *p){
		CCell<TParticle> *cell=which(*p);
		cell->add(p);
		p->cell=cell;
		}

	void update(TParticle *p){
		CCell<TParticle> *cell=which(*p);
		if(cell==p->cell)return;
		ERROR(!p->cell, "Particle's cell pointer not specified");
		(static_cast<CCell<CParticle>*>(p->cell))->remove(p);
		cell->add(p);
		p->cell=cell;
		}

	void interact(){
		for(int i=0; i<nx*ny*nz; i++){
			nodes[i].interact();
			}
		}

	void print(ostream &out=cout)const{
		out<< "lines"<<endl;
		double x=c(0);
		while(x<=c(0)+diag(0)+1e-10){
			out<< x<<"  "<< c(1)<<"\t"<< x <<"  "<<c(1)+diag(1)<<endl;
			x+=dx;
			}
		double y=c(1);
		while(y<=c(1)+diag(1)+1e-10){
			out<< c(0)<<"  "<< y<<"\t"<< c(0)+diag(0)<<"  "<<y<<endl;
			y+=dy;
			}
		}

	CCell<TParticle> *boundary_mask(int i,int j, int k,  vec &shift, bool &is_shifted);
	CCell<TParticle> *node(int i, int j, int k)const{
		return &nodes[nz*(ny*i+j)+k];
		}
 	//private:
	void build_neighbors();
	int nx, ny, nz;
	vec c, diag;
	double dx, dy, dz;
	CCell<TParticle> *nodes;
 	private:
	bool periodic_x, periodic_y, periodic_z, is_set_up;
	TParticleContainer particles;
	};


template<typename TParticleContainer, typename TParticle>
CCell<TParticle> *CCellList<TParticleContainer, TParticle>::boundary_mask(int i,int j, int k, vec &shift, bool &is_shifted){
	
		cerr<< i<<"  "<<j<<"  "<<k <<endl;
		shift*=0;
		is_shifted=false;
		if(i >= nx){
			if(periodic_x){
				i-=nx;
				shift(0)-=diag(0);
				is_shifted=true;
				}
			else return NULL;
			}
		else if(i < 0){
			if(periodic_x){
				i+=nx;
				shift(0)+=diag(0);
				is_shifted=true;
				}
			else return NULL;
			}
		if(j >= ny){
			if(periodic_y){
				j-=ny;
				shift(1)-=diag(1);
				is_shifted=true;
				}
			else return NULL;
				}
		else if(j < 0){
			if(periodic_y){
				j+=ny;
				shift(1)+=diag(1);
				is_shifted=true;
				}
			else return NULL;
				}

		if(k >= nz){
			if(periodic_z){
				k-=nz;
				shift(2)-=diag(2);
				is_shifted=true;
				}
			else return NULL;
				}

		else if(k < 0){
			if(periodic_z){
				k+=nz;
				shift(2)+=diag(2);
				is_shifted=true;
				}
			else return NULL;
				}
	
		return node(i,j, k);
		}

//only half of the neighbours 
//needs to be saved (should be carefully chosen such that it is the complement of the other half)
template<typename TParticleContainer, typename TParticle>
void CCellList<TParticleContainer, TParticle>::build_neighbors(){
	for(int i=0; i<nx; i++){
	for(int j=0; j<ny; j++){
	for(int k=0; k<nz; k++){
		vec shift(0,0,0);	
		bool is_shifted=false;
		CCell<TParticle> *p, *pij;
		pij=node(i,j,k);
		
		p=boundary_mask(i+1,j-1,k, shift, is_shifted);
		pij->neighs[1]=p;
		pij->shifts[1]=shift;
		pij->shifted[1]=is_shifted;

		p=boundary_mask(i+1,j,k, shift, is_shifted);
		pij->neighs[2]=p;
		pij->shifts[2]=shift;
		pij->shifted[2]=is_shifted;
		
		p=boundary_mask(i,j+1,k, shift, is_shifted);
		pij->neighs[3]=p;
		pij->shifts[3]=shift;
		pij->shifted[3]=is_shifted;

		p=boundary_mask(i+1,j+1,k, shift, is_shifted);
		pij->neighs[4]=p;
		pij->shifts[4]=shift;
		pij->shifted[4]=is_shifted;

		// ----------------

		p=boundary_mask(i-1,j-1,k+1, shift, is_shifted);
		pij->neighs[5]=p;
		pij->shifts[5]=shift;
		pij->shifted[5]=is_shifted;

		p=boundary_mask(i-1,j,  k+1, shift, is_shifted);
		pij->neighs[6]=p;
		pij->shifts[6]=shift;
		pij->shifted[6]=is_shifted;
		
		p=boundary_mask(i-1,j+1,k+1, shift, is_shifted);
		pij->neighs[7]=p;
		pij->shifts[7]=shift;
		pij->shifted[7]=is_shifted;

		p=boundary_mask(i,j-1,	k+1, shift, is_shifted);
		pij->neighs[8]=p;
		pij->shifts[8]=shift;
		pij->shifted[8]=is_shifted;

		p=boundary_mask(i,j,  	k+1, shift, is_shifted);
		pij->neighs[9]=p;
		pij->shifts[9]=shift;
		pij->shifted[9]=is_shifted;
		
		p=boundary_mask(i,j+1,	k+1, shift, is_shifted);
		pij->neighs[10]=p;
		pij->shifts[10]=shift;
		pij->shifted[10]=is_shifted;

		p=boundary_mask(i+1,j-1,k+1, shift, is_shifted);
		pij->neighs[11]=p;
		pij->shifts[11]=shift;
		pij->shifted[11]=is_shifted;

		p=boundary_mask(i+1,j,  k+1, shift, is_shifted);
		pij->neighs[12]=p;
		pij->shifts[12]=shift;
		pij->shifted[12]=is_shifted;
		
		p=boundary_mask(i+1,j+1,k+1, shift, is_shifted);
		pij->neighs[13]=p;
		pij->shifts[13]=shift;
		pij->shifted[13]=is_shifted;

		}
		}
		}
	}

#endif /* CELLLIST_H */
