// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef SIZE_DIST_H
#define SIZE_DIST_H 
#include<iostream>
#include<string>


#include<random>
// TR1's ranlux64_base_01 was a floating-point subtract-with-carry engine that
// C++11 dropped.  ranlux48_base is the standard's equivalent: same 48-word lag
// and 5-word short lag.  NOTE this changes the random stream, so packings
// drawn from a SizeDistribution differ from a pre-C++17 build -- the physics
// is the same, the realisation is not.  See docs/PORTING-NOTES.md.
std::ranlux48_base eng0;


class CBaseDistribution
	{
	public:
	CBaseDistribution(string _name):name(_name), max_value(0){}
	// Virtual, because CSizeDistribution holds these by base pointer and
	// deletes through it.  Without this, `delete p_dist` was undefined
	// behaviour -- latent only because the destructor was inverted and never
	// actually deleted anything (see below).
	virtual ~CBaseDistribution(){}
	// Deep copy, so CSizeDistribution can have real value semantics.  It used
	// to share the pointer, which is why the inverted destructor could not
	// simply be fixed: flipping it turned a leak into a double free.
	virtual CBaseDistribution *clone() const = 0;
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
	// Rule of three.  The compiler-generated copy shared p_dist, so the old
	// destructor had been written inverted -- `if(!p_dist) delete p_dist` --
	// which leaked rather than double-freeing.  With a deep copy the obvious
	// destructor is correct.
	CSizeDistribution(const CSizeDistribution &other)
		:p_dist(other.p_dist ? other.p_dist->clone() : NULL){}
	CSizeDistribution &operator=(const CSizeDistribution &other){
		if(this!=&other){
			CBaseDistribution *copy =
				other.p_dist ? other.p_dist->clone() : NULL;
			delete p_dist;
			p_dist=copy;
			}
		return *this;
		}
	~CSizeDistribution(){ delete p_dist; };

	friend istream & operator>>(istream &in, CSizeDistribution &dist);
	friend std::ostream & operator<< (std::ostream &out, const CSizeDistribution &dist);

	void print(ostream &out)const{
		p_dist->print(out);
		}
	double get(){
		ERROR(!p_dist, "size distribution used before it was set");
		return p_dist->get();
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
	CBaseDistribution *clone() const { return new CMonoDist(*this); }
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
	return in;
	}
ostream & operator<< (ostream &out, const CMonoDist &dist){
	dist.print(out);
	return out;
	}


class CUniformDist: public CBaseDistribution
	{
	public:
	CUniformDist(double _min=0.5, double _max=1):CBaseDistribution("uniform"), min(_min), max(_max)
		{
		unif=std::uniform_real_distribution<double> (min, max);
		}
	CBaseDistribution *clone() const { return new CUniformDist(*this); }
	double get(){
		return unif(eng0);
		}
	void parse(istream &in){
		in>>min>>max;
		unif=std::uniform_real_distribution<double> (min, max);
		max_value=max;
		}

	void print(ostream &out)const{
		out<<name<<"  "<<min<<"  "<<max;
		}

	friend istream & operator>>(istream &in, CUniformDist &dist);
	friend ostream & operator<< (ostream &out, const CUniformDist &dist);

 	private:
	// min/max are declared and initialised before unif: the constructor used to
	// build unif from min and max in the member-initialiser list, but members
	// are initialised in *declaration* order, so unif was being constructed
	// from uninitialised bounds.
	double min, max;
	std::uniform_real_distribution<double> unif;
	};
istream & operator>>(istream &in, CUniformDist &dist){
	dist.parse(in);
	return in;
	}
ostream & operator<< (ostream &out, const CUniformDist &dist){
	dist.print(out);
	return out;
	}


class CReadDist: public CBaseDistribution
	{
	public:
	CReadDist():CBaseDistribution("read") {}
	CBaseDistribution *clone() const { return new CReadDist(*this); }
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
		unif=std::uniform_int_distribution<int> (0, values.size()-1);
		}

	void print(ostream &out)const{
		out<<name<<"  from  "<<filename;
		}

	friend istream & operator>>(istream &in, CUniformDist &dist);
	friend ostream & operator<< (ostream &out, const CUniformDist &dist);

 	private:
	std::uniform_int_distribution<int> unif;
	vector<double> values;
	string filename;
	};
istream & operator>>(istream &in, CReadDist &dist){
	dist.parse(in);
	return in;
	}
ostream & operator<< (ostream &out, const CReadDist &dist){
	dist.print(out);
	return out;
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
	return in;
}
ostream & operator<<(ostream &out, const CSizeDistribution &dist){
	dist.print(out);
	return out;
}
#endif /* SIZE_DIST_H */
