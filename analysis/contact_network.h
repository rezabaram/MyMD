#ifndef CONTACT_NETWORK_H
#define CONTACT_NETWORK_H 

template<class T>
class ContactNetwork 
	{
	class Contact : public BasicContact{ 
		public: 
		Contact(T *_p1, T *_p2, vec _x, vec _n):BasicContact(_x, _n), p1(_p1), p2(_p2){} ; 
		Contact(T *_p1, T *_p2, const BasicContact &bc):BasicContact(bc), p1(_p1), p2(_p2){} ; 
		T * p1, *p2; 
		}; 
	public:
	ContactNetwork(vector<T *> *_p):N(0), packing(_p){}
	~ContactNetwork(){
		//delete [] elems;
		}

	void build(){
		typename vector<T *>::const_iterator it;
		N=packing->size();
		//elems= new vector<Contact>[N];
		for(size_t i=0;  i<N; i++){
			elems.push_back(vector<Contact>());
			}
		//assert(elems);
		for(size_t i=0;  i<N; i++){
			for(size_t j=i+1;  j<N; j++){
				ShapeContact ovs(packing->at(i),packing->at(i));
				//ovs.back();
				CInteraction::overlaps(&ovs, packing->at(i), packing->at(j) ); 
				if(!ovs.empty()) {
					assert(ovs.size()==1);
					elems[i].push_back(Contact(packing->at(i), packing->at(j), static_cast<BasicContact>(ovs.back()) ));
					elems[j].push_back(Contact(packing->at(j), packing->at(i), static_cast<BasicContact>(ovs.back()) ));
					}
				}
			}

		/*
		for(size_t i=0;  i<NContactNetwork; i++){
		for(size_t j=0; j<elems[i].size(); j++){
			elems[i].at(j).l
			}
		*/
	}

	void print(ostream &out)const{
		//assert(elems);
		assert(N);
		cerr<< N <<endl;
		for(size_t i=0; i<N; i++){
			out<<"2  ";
			out<<packing->at(i)->Xc<<"  0.008"<<endl;;
		for(size_t j=0; j<elems[i].size(); j++){
			//if(elems[i].at(j).p1->Xc(1)>0.12)continue;
			//out<<"2  ";
			//out<<elems[i].at(j).x<<"  0.005  0 0 50000"<<endl;;

			out<<"3  ";
			out<<elems[i].at(j).p1->Xc<<"  0.001  ";
			out<<elems[i].at(j).x<<"  0.001"<<endl;;

			}
			}
	}
		

	size_t N;
	private:
	vector < vector<Contact> > elems;
	vector<T *> *packing;
	};
#endif /* CONTACT_NETWORK_H */
