#ifndef PACKING_H
#define PACKING_H 
#include"exception.h"
#include"common.h"
#include"shapes.h"
#include"grid.h"
#include"ellips_contact.h"
#include"MersenneTwister.h"
#include"interaction.h"
#include"contact_network.h"




//member classes for contact network

template < class T>
class CPacking : public vector <T *>
	{

	public:
	CPacking():contacts(ContactNetwork<T>(this))
		{
		maxr=0.0;
		grid_built=false;
		}
	~CPacking(){
		typename vector<T *>::iterator it;
		if(grid_built)delete grid;
		for(it=this->begin(); it!=this->end(); it++)
			delete (*it);
		}

	void print(std::ostream& out, bool raster=false)const;
	void printRaster3D(std::ostream& out)const{print(out, true);}
	void printEuler(std::ostream& out)const;
	void parse(string infilename);
	void parse(istream &inputFile);
	double packFraction(const vec &x1, const vec &x2, unsigned long N=10000);
	void output_contact_network(ostream &out);
	double totalVolume();
	void BuildContactNetwork();
	void BuildGrid(){
		if(grid_built)return;
		grid=(new CRecGrid<T>(this, vec(0.0, 0.0, 0.0), vec(1.5, 1.5, 1.5), 1.1*maxr));
		grid->build();
		grid_built=true;
		};

	ContactNetwork<T> contacts;
	
	CRecGrid<T> * grid;
 	private:
	double maxr;
	bool grid_built;
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
	BuildGrid();
	double n=0, nt=0;
	vec x;
	//vec l=x2-x1;
	//double boxvolume=fabs(l(0)*l(1)*l(2));
	double d=1;
	typename CNode3D<T>::const_iterator it;

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
		CEllipsoid shape;
		T *p=new T(shape);
		p->shape->parse(ss);
		this->push_back(p);
		if(p->shape->radius>maxr)maxr=p->shape->radius;
		}

	//else if(id==14){
		//GeomObjectBase *shape=new GeomObject<tsphere>();
		//shape->parse(ss);
		//push_back(shape);
		//}
	else{
		//WARNING("Skip reading shape with id: "<< id);
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

#endif /* PACKING_H */
