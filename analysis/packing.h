#ifndef PACKING_H
#define PACKING_H 
#include"../exception.h"
#include"../grid.h"
#include"../ellips_contact.h"
#include"../MersenneTwister.h"
#include"../interaction.h"



vec randomVec(const vec &x1, const vec &x2){
	return vec( x1(0)+(x2(0)-x1(0))*rgen(), x1(1)+(x2(1)-x1(1))*rgen(), x1(2)+(x2(2)-x1(2))*rgen());
	}

//member classes for contact network
template<class T>
class ContactNetwork;

template < class T>
class CPacking : public vector <T *>
	{

	public:
	CPacking():contacts(ContactNetwork<T>(this)){
		}
	~CPacking(){
		typename vector<T *>::iterator it;
		for(it=this->begin(); it!=this->end(); it++)
			delete (*it);
		}

	void print(std::ostream& out, bool raster=false)const;
	void printRaster3D(std::ostream& out)const{print(out, true);}
	void parse(string infilename);
	void parse(istream &inputFile);
	double packFraction(const vec &x1, const vec &x2, unsigned long N=10000);
	void output_contact_network(ostream &out);
	double totalVolume();
	void BuildContactNetwork();

	ContactNetwork<T> contacts;
	
 	private:
	};

template < class T>
double CPacking<T>::totalVolume(){
	double v=0;
	typename CPacking<T>::const_iterator it;
	for(it=this->begin(); it!=this->end(); it++){
		v+=(**it).vol();
		}
	return v;
	}

//calculating packing fraction inside cube with x1, x2 az opposite corners
template < class T>
double CPacking<T>::packFraction(const vec &x1, const vec &x2, unsigned long N){
TRY
	double n=0, nt=0;
	vec x;
	double d=1;
	typename CPacking<T>::const_iterator it;

	for(unsigned long i=0; i<N; i++){
		++nt;
		x=randomVec(x1,x2);
		for(it=this->begin(); it!=this->end(); it++){
			d=(**it)(x);
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
void CPacking<T>::print(std::ostream& out, bool raster)const{
		typename CPacking::const_iterator it;
		for(it=this->begin(); it!=this->end(); it++){
			if(raster) (*it)->printRaster3D(out);
			else (*it)->print(out);
			out<<endl;
			}
		}

template < class T>
void CPacking<T>::parse(string infilename) {
	ifstream inputFile(infilename.c_str());
	if(!inputFile.good())
	{
	cerr << "Unable to open input file: " << infilename << endl;
	return;
	}
	parse(inputFile);
	inputFile.close();
	}

template < class T>
void CPacking<T>::parse(istream &inputFile) {

	string line;

	//Parse the line
	while(getline(inputFile,line))
	{

	//Insert the line string into a stream
	stringstream ss(line);

	int id;
	ss>>id;
	if(id==14){
		T *shape=new T();
		shape->parse(ss);
		this->push_back(shape);
		}

	//else if(id==14){
		//GeomObjectBase *shape=new GeomObject<tsphere>();
		//shape->parse(ss);
		//push_back(shape);
		//}
	else{
		WARNING("Skip reading shape with id: "<< id);
		continue;
		}
		

	//Read up to second whitespace
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

/*
void print_connectivity(vector<CCircle> &packing, const char * name){
  // !!! C++ style Numbering starting with 0
  ofstream outconnect(name);
  outconnect<<packing.size()<<endl;
  for(int i=0; i<packing.size(); i++){
    outconnect<< i<<"  "<< (double)(packing.at(i).color)*packing.at(i).radius
              <<"  "<<packing.at(i).neighbours.size()<<"  ";
    for(int j=0; j<packing.at(i).neighbours.size(); j++){
      outconnect<<packing.at(i).neighbours.at(j)<<"  ";
    }
    outconnect<<endl;
  }
}

*/

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
	ContactNetwork(CPacking<T> *_p):N(0), elems(NULL), packing(_p){}
	~ContactNetwork(){
		delete [] elems;
		}

	void build(){
		typename CPacking<T>::const_iterator it;
		N=packing->size();
		elems= new vector<Contact>[N];
		assert(elems);
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
		assert(elems);
		assert(N);
		cerr<< N <<endl;
		for(size_t i=0; i<N; i++){
			out<<"2  ";
			out<<packing->at(i)->Xc<<"  0.008  0 0 50000"<<endl;;
		for(size_t j=0; j<elems[i].size(); j++){
			//if(elems[i].at(j).p1->Xc(1)>0.12)continue;
			//out<<"2  ";
			//out<<elems[i].at(j).x<<"  0.005  0 0 50000"<<endl;;

			out<<"5  ";
			out<<elems[i].at(j).p1->Xc<<"  0.001  ";
			out<<elems[i].at(j).x<<"  0.005  50000 50000 50000"<<endl;;

			}
			}
	}
		

	size_t N;
	private:
	vector<Contact> *elems;
	CPacking<T> *packing;
	};
#endif /* PACKING_H */
