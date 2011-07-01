#ifndef CONTACT_NETWORK_H
#define CONTACT_NETWORK_H 
#include<vector>
#include<list>
#include"matrix.h"
#include"eigsys.h"

using namespace std;

template<class T>
class TContact : public BasicContact{ 
	public: 
	TContact(T *_p1, T *_p2, vec _x, vec _n):BasicContact(_x, _n), p1(_p1), p2(_p2),l(_x-_p1->shape->Xc){} ; 
	TContact(T *_p1, T *_p2, const BasicContact &bc):BasicContact(bc), p1(_p1), p2(_p2),l(bc.x-_p1->shape->Xc){} ; 
	T * p1, *p2; 
	vec l;
	}; 

template<class T>
class TNode : public vector<TContact<T> >{
	public:
	TNode():F(Matrix(3,3)){};
	Matrix &cal_fabric_tensor()
	{
		typename TNode<T>::iterator it;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				F(i,j)=0;

		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				for(it=this->begin(); it!=this->end(); it++)
					F(i,j)+=(*it).l(i)*(*it).l(j);
		return F;
	}

	Matrix F;
	};

template<class T>
class ContactNetwork : public vector< TNode<T> >
	{
	public:
	ContactNetwork(list<T *> *_p): F(Matrix(3,3)), N(0), packing(_p){}
	~ContactNetwork(){
		}

	Matrix cal_fabric_tensor(){
		double N=this->size();

		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				F(i,j)=0;

		for(size_t i=0;  i<N; i++){
			F+=this->at(i).cal_fabric_tensor();
			}
		F=F/N;
		return F;
		}
	void build(){
		typename vector<T *>::const_iterator it;
		N=packing->size();

		for(size_t i=0;  i<N; i++){
			this->push_back(TNode<T>());
			}

		typename list<T *>::const_iterator it1, it2;
		size_t i=0, j=0;
		for(it1=packing->begin(), i=0;  it1!=packing->end(); it1++, ++i){
			for(it2=it1, ++it2, j=i+1;  it2!=packing->end(); it2++, ++j){
				ShapeContact ovs((*it1)->shape,(*it2)->shape);
				CInteraction::overlaps(&ovs, (*it1)->shape, (*it2)->shape ); 
				if(!ovs.empty()) {
					assert(ovs.size()==1);
					this->at(i).push_back(TContact<T>((*it1), (*it2), static_cast<BasicContact>(ovs.back()) ));
					this->at(j).push_back(TContact<T>((*it2), (*it1), static_cast<BasicContact>(ovs.back()) ));
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

	void print_eigen(ostream &out){
		CEigSys eigsys(cal_fabric_tensor());
		eigsys.print_eigen_vals(cout);
		typename ContactNetwork<T>::iterator it1;
		size_t i=0;
		for(it1=this->begin(), i=0;  it1!=this->end(); it1++, ++i){
			CEigSys eigsys((*it1).F);
			eigsys.print_eigen_vals(cout);
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
		

	Matrix F;//fabric tensor
	size_t N;
	private:
	list<T *> *packing;
	};
#endif /* CONTACT_NETWORK_H */
