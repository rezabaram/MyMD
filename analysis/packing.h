#ifndef PACKING_H
#define PACKING_H 
#include"../grid.h"
#include"../ellips_contact.h"
#include"../MersenneTwister.h"
#include"../interaction.h"


vec randomVec(const vec &x1, const vec &x2){
	return vec( x1(0)+(x2(0)-x1(0))*rgen(), x1(1)+(x2(1)-x1(1))*rgen(), x1(2)+(x2(2)-x1(2))*rgen());
	}

template < class T>
class CPacking : public vector <T *>
	{

	//member classes for contact network
	class Contact : public BasicContact{ public: Contact(T *_p, vec _x, vec _n):BasicContact(_x, _n), p(_p){} ; T * p; }; 
	class ContactNetwork : public vector < vector <Contact> >{};

	public:
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
	double totalVolume();
	void BuildContactNetwork();

	ContactNetwork contacts;
	
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
		contacts.clear();
		typename CPacking::const_iterator it1, it2;
		for(it1=this->begin();  it1<this->end(); it1++){
				contacts.push_back(vector<Contact>() );
			for(it2=it1+1;  it2<this->end(); it2++){
				ShapeContact ovs(*it1, *it2);
				//ovs.back();
				CInteraction::overlaps(&ovs, *it1, *it2); 
				if(!ovs.empty()) {
						assert(ovs.size()==1);
						contacts.back()->push_back(Contact(*ovs.back()));
					}
		
				}
			}
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
#endif /* PACKING_H */
