#ifndef SIZE_DIST_H
#define SIZE_DIST_H 
#include<iostream>
#include<string>


#include<tr1/random>
std::tr1::ranlux64_base_01 eng0;


class CBaseDistribution
	{
	public:
	CBaseDistribution(string _name):name(_name), max_value(0){}
	virtual double get()=0;
	virtual void parse(istream &in)=0;
	virtual void print(ostream &out)const=0;
	double get_max(){return max_value;}
	string name;
	protected:
	double max_value;
 	private:
	};

//a primitive factory 
class CSizeDistribution
	{
	public:
	CSizeDistribution():p_dist(NULL){}
	template<class T>
	CSizeDistribution(const T &dist):p_dist(new T(dist)){}
	~CSizeDistribution(){if(!p_dist) delete p_dist;};

	friend istream & operator>>(istream &in, CSizeDistribution &dist);
	friend std::ostream & operator<< (std::ostream &out, const CSizeDistribution &dist);

	void print(ostream &out)const{
		p_dist->print(out);
		}
	double get(){
		p_dist->get();
		}

	string get_name(){
		return p_dist->name;
		}

 	private:
	CBaseDistribution *p_dist;
	};


class CMonoDist : public CBaseDistribution
	{
	public:
	CMonoDist(double _r=0):CBaseDistribution("mono"), r(_r){}
	double get(){return r;}
	void parse(istream &in){
		in>>r;
		max_value=r;
		}

	void print(ostream &out)const{
		out<<name<<"  "<<r;
		}

	friend istream & operator>>(istream &in, CMonoDist &dist);
	friend ostream & operator<< (ostream &out, const CMonoDist &dist);

 	private:
	double r;
	};
istream & operator>>(istream &in, CMonoDist &dist){
	dist.parse(in);
	}
ostream & operator<< (ostream &out, const CMonoDist &dist){
	dist.print(out);
	}


class CUniformDist: public CBaseDistribution
	{
	public:
	CUniformDist(double _min=0.5, double _max=1):CBaseDistribution("uniform"), min(_min), max(_max)
		,unif(tr1::uniform_real<double> (min, max))
		{
		}
	double get(){
		return unif(eng0);
		}
	void parse(istream &in){
		in>>min>>max;
		unif=tr1::uniform_real<double> (min, max);
		max_value=max;
		}

	void print(ostream &out)const{
		out<<name<<"  "<<min<<"  "<<max;
		}

	friend istream & operator>>(istream &in, CUniformDist &dist);
	friend ostream & operator<< (ostream &out, const CUniformDist &dist);

 	private:
	tr1::uniform_real<double> unif;
	double min, max;
	};
istream & operator>>(istream &in, CUniformDist &dist){
	dist.parse(in);
	}
ostream & operator<< (ostream &out, const CUniformDist &dist){
	dist.print(out);
	}


class CReadDist: public CBaseDistribution
	{
	public:
	CReadDist():CBaseDistribution("read") {}
	double get(){
		return values.at(unif(eng0));
		}
	void parse(istream &in){
		in>>filename;
		read(filename);
		}
	
	void read(string file){ 
		ifstream ifile(file.c_str());
		string line;
		//Parse the line
		double r;
		while(getline(ifile,line)){
			stringstream ss(line);
			ss>>r;
			values.push_back(r);
			max_value=max(max_value,r);
			}
		unif=tr1::uniform_int<int> (0, values.size()-1);
		}

	void print(ostream &out)const{
		out<<name<<"  from  "<<filename;
		}

	friend istream & operator>>(istream &in, CUniformDist &dist);
	friend ostream & operator<< (ostream &out, const CUniformDist &dist);

 	private:
	tr1::uniform_int<int> unif;
	vector<int> values;
	string filename;
	};
istream & operator>>(istream &in, CReadDist &dist){
	dist.parse(in);
	}
ostream & operator<< (ostream &out, const CReadDist &dist){
	dist.print(out);
	}



istream & operator>>(istream &in, CSizeDistribution &dist){
	string name;
	in>>name;
	if(name=="mono") { 
		dist.p_dist=new CMonoDist();
		dist.p_dist->parse(in);
		}
	else if(name=="uniform") { 
		dist.p_dist=new CUniformDist();
		dist.p_dist->parse(in);
		}
	else if(name=="read") { 
		dist.p_dist=new CReadDist();
		dist.p_dist->parse(in);
		}
	else{
		ERROR(1,"distribution"+name+"not defined!");
		}
}
ostream & operator<<(ostream &out, const CSizeDistribution &dist){
	dist.print(out);
	}
#endif /* SIZE_DIST_H */
