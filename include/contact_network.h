#ifndef CONTACT_NETWORK_H
#define CONTACT_NETWORK_H 
#include<vector>
#include<list>
#include"matrix.h"

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
	ContactNetwork(vector<T *> *_p): F(Matrix(3,3)), N(0), packing(_p){}
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

		for(size_t i=0;  i<N; i++){
			for(size_t j=i+1;  j<N; j++){
				ShapeContact ovs(packing->at(i)->shape,packing->at(i)->shape);
				CInteraction::overlaps(&ovs, packing->at(i)->shape, packing->at(j)->shape ); 
				if(!ovs.empty()) {
					assert(ovs.size()==1);
					this->at(i).push_back(TContact<T>(packing->at(i), packing->at(j), static_cast<BasicContact>(ovs.back()) ));
					this->at(j).push_back(TContact<T>(packing->at(j), packing->at(i), static_cast<BasicContact>(ovs.back()) ));
					}
				}
			}

		/*
		for(size_t i=0;  i<NContactNetwork; i++){
		for(size_t j=0; j<this->[i].size(); j++){
			this->[i].at(j).l
			}
		*/
	}

	void print(ostream &out)const{
		//assert(elems);
		assert(N);
		cerr<< N <<endl;
		for(size_t i=0; i<N; i++){
			out<<"2  ";
			out<<packing->at(i)->shape->Xc<<"  0.008"<<endl;;
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
	vector<T *> *packing;
	};
#endif /* CONTACT_NETWORK_H */
