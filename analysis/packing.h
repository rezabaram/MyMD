#ifndef PACKING_H
#define PACKING_H 
#include"../ellipsoid.h"
#include"../MersenneTwister.h"


vec randomVec(const vec &x1, const vec &x2){
	return vec( x1(0)+(x2(0)-x1(0))*rgen(), x1(1)+(x2(1)-x1(1))*rgen(), x1(2)+(x2(2)-x1(2))*rgen());
	}

class CPacking : public vector <GeomObjectBase *>
	{
	public:
	~CPacking(){
		CPacking::iterator it;
		for(it=begin(); it!=end(); it++)
			delete (*it);
		}

	void print(std::ostream& out, bool raster=false)const;
	void printRaster3D(std::ostream& out)const{print(out, true);}
	void parse(string infilename);
	void parse(istream &inputFile);
	double packFraction(const vec &x1, const vec &x2, unsigned long N=10000);
	double totalVolume();

 	private:
	};

double CPacking::totalVolume(){
	double v=0;
	CPacking::const_iterator it;
	for(it=this->begin(); it!=this->end(); it++){
		v+=(**it).vol();
		}
	return v;
	}
//calculating packing fraction inside cube with x1, x2 az opposite corners
double CPacking::packFraction(const vec &x1, const vec &x2, unsigned long N){
TRY
	double n=0, nt=0;
	vec x;
	double d=1;
	CPacking::const_iterator it;

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
void CPacking::print(std::ostream& out, bool raster)const{
		CPacking::const_iterator it;
		for(it=this->begin(); it!=this->end(); it++){
			if(raster) (*it)->printRaster3D(out);
			else (*it)->print(out);
			out<<endl;
			}
		}

void CPacking::parse(string infilename) {
	ifstream inputFile(infilename.c_str());
	if(!inputFile.good())
	{
	cerr << "Unable to open input file: " << infilename << endl;
	return;
	}
	parse(inputFile);
	inputFile.close();
	}

void CPacking::parse(istream &inputFile) {

	string line;

	//Parse the line
	while(getline(inputFile,line))
	{

	//Insert the line string into a stream
	stringstream ss(line);

	int id;
	ss>>id;
	if(id==14){
		GeomObjectBase *shape=new GeomObject<tellipsoid>();
		shape->parse(ss);
		push_back(shape);
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
#endif /* PACKING_H */
