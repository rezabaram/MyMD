#ifndef CONTACT_NETWORK_H
#define CONTACT_NETWORK_H 
#include<vector>
#include<list>
#include"matrix.h"
#include"eigsys.h"

using namespace std;

template<class T>
class TContact : public Contact{ 
	public: 
	//TContact(T *_p1, T *_p2, vec _x, vec _n):Contact(_x, _n), p1(_p1), p2(_p2),l(_x-_p1->shape->Xc){} ; 
	TContact(T *_p1, T *_p2, const Contact &bc):Contact(bc), p1(_p1), p2(_p2),l(bc.x-_p1->shape->Xc), n(bc.n),fn(bc.n*pow(bc.dx_n,1.5) ){} ; 
	T * p1, *p2; 
	vec l, n, fn;
	}; 

template<class T>
class TNode : public vector<TContact<T> >{
	public:
	TNode(T *_p):p(_p),branch_fabric_M(Matrix(3,3)), normal_fabric_M(Matrix(3,3)){};
	Matrix &branch_fabric_tensor()
	{
		bool calculated=false;
		if(calculated) return branch_fabric_M;
		typename TNode<T>::iterator it;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				branch_fabric_M(i,j)=0;

		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				for(it=this->begin(); it!=this->end(); it++)
					branch_fabric_M(i,j)+=(*it).l(i)*(*it).l(j);
		calculated=true;
		return branch_fabric_M;
	}

	Matrix &normal_fabric_tensor()
	{
		bool calculated=false;
		if(calculated) return branch_fabric_M;
		typename TNode<T>::iterator it;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				normal_fabric_M(i,j)=0;

		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				for(it=this->begin(); it!=this->end(); it++)
					normal_fabric_M(i,j)+=(*it).n(i)*(*it).n(j);
		calculated=true;
		return normal_fabric_M;
	}

	void print_branch_vectors(ostream &out,const vec3d &c1=vec3d(0,0,0),const vec3d &c2=vec3d(1,1,1)){
		typename TNode<T>::iterator it;
		for(it=this->begin(); it!=this->end(); it++){
			vec3d x=it->x;
			if( x(0)<c1(0) or 
			x(0)>c2(0) or 
			x(1)<c1(1) or 
			x(1)>c2(1) or 
			x(2)<c1(2) or 
			x(2)>c2(2) ) continue;
			out<<spherical((*it).n)<<endl;
			}
			//if((*it).n*vec3d(0,0,1)>0)out<<spherical((*it).n)<<endl;
		}

	
	T *p;
	private:
	Matrix branch_fabric_M;
	Matrix normal_fabric_M;
	};

template<class T>
class ContactNetwork : public vector< TNode<T> >
	{
	public:
	ContactNetwork(list<T *> *_p): branch_fabric_M(Matrix(3,3)), normal_fabric_M(Matrix(3,3)), N(0), packing(_p){}
	~ContactNetwork(){
		}

	Matrix branch_fabric_tensor(){
		double N=this->size();

		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++){
				branch_fabric_M(i,j)=0;
				}

		for(size_t i=0;  i<N; i++){
			branch_fabric_M+=this->at(i).branch_fabric_tensor();
			}
		branch_fabric_M=branch_fabric_M/N;
		return branch_fabric_M;
		}
	Matrix normal_fabric_tensor(){
		double N=this->size();

		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++){
				normal_fabric_M(i,j)=0;
				}

		for(size_t i=0;  i<N; i++){
			normal_fabric_M+=this->at(i).normal_fabric_tensor();
			}
		normal_fabric_M=normal_fabric_M/N;
		return normal_fabric_M;
		}

	void build(){
		N=packing->size();

		typename list<T *>::const_iterator it;
		for(it=packing->begin();  it!=packing->end(); it++){
			this->push_back(TNode<T>(*it));
			}

		typename list<T *>::const_iterator it1, it2;
		size_t i=0, j=0;
		for(it1=packing->begin(), i=0;  it1!=packing->end(); it1++, ++i){
			ERROR(*it1!=this->at(i).p, "In consistency in contact network");
			for(it2=it1, ++it2, j=i+1;  it2!=packing->end(); it2++, ++j){
				ShapeContact ovs((*it1)->shape,(*it2)->shape);
				CInteraction::overlaps(&ovs, (*it1)->shape, (*it2)->shape ); 
				if(!ovs.empty()) {
					assert(ovs.size()==1);
					Contact temp=ovs.back();
					this->at(i).push_back(TContact<T>((*it1), (*it2), temp ));
					temp.n*=-1.0;
					this->at(j).push_back(TContact<T>((*it2), (*it1), temp ));
					ERROR(*it2!=this->at(j).p, "Inconsistency in contact network");
					}
				}
			}

	}

	double avg_contact_number()const{
		double total=0;
		for(unsigned int i=0; i<this->size(); i++){
			total+=this->at(i).size();
			}
		return total/(double)this->size();
		
		}

	void print_branch_vectors(ostream &out,const vec3d &box_c1=vec3d(0,0,0),const vec3d &box_c2=vec3d(1,1,1)){
		typename ContactNetwork<T>::iterator it1;
		for(it1=this->begin();  it1!=this->end(); it1++){
			it1->print_branch_vectors(out,box_c1,box_c2);
			}
		}
	void print_eigen(ostream &out ,const vec3d &c1=vec3d(0,0,0),const vec3d &c2=vec3d(1,1,1)){
		//CEigSys eigsys(branch_fabric_tensor());
		//eigsys.print_eigen_vals(cout);
		typename ContactNetwork<T>::iterator it1;
		size_t i=0;
		for(it1=this->begin(), i=0;  it1!=this->end(); it1++, ++i){
			vec3d x=it1->p->shape->Xc;
			if(it1->p->is_shadow) continue;
			if( x(0)<c1(0) or 
			x(0)>c2(0) or 
			x(1)<c1(1) or 
			x(1)>c2(1) or 
			x(2)<c1(2) or 
			x(2)>c2(2) ) continue;
			CEigSys eigsys1((*it1).branch_fabric_tensor());
			CEigSys eigsys2((*it1).normal_fabric_tensor());
			eigsys1.print_eigens(cout);
			eigsys2.print_eigens(cout);
			cout<<endl;
			}
		}


	void print(ostream &out)const{
		//assert(elems);
		assert(N);
		cerr<< N <<endl;
		typename list<T *>::const_iterator it1;
		size_t i=0;
		for(it1=packing->begin(), i=0;  it1!=packing->end(); it1++, ++i){
			out<<"2  ";
			out<<(*it1)->shape->Xc<<"  0.008"<<endl;;
		for(size_t j=0; j<this->at(i).size(); j++){
			//if(this->[i].at(j).p1->shape->Xc(1)>0.12)continue;
			//out<<"2  ";
			//out<<this->[i].at(j).x<<"  0.005  0 0 50000"<<endl;;

			out<<"3  ";
			out<<this->at(i).at(j).p1->shape->Xc<<"  0.001  ";
			out<<this->at(i).at(j).x<<"  0.001"<<endl;;

			}
			}
	}
		

	Matrix branch_fabric_M;//fabric tensor of branch vectors (center to contact points)
	Matrix normal_fabric_M;//fabric tensor of contact normals
	size_t N;
	private:
	list<T *> *packing;
	};
#endif /* CONTACT_NETWORK_H */
