// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


/*By: 	Reza Baram
CopyRight: 
	Do whatever you want with this code. You can even replace 
	my name with yours. 

Disclaimer:
	There is no guarantee whatsoever that this code will work 
	perfectly as it is intended to. So, use is at your own
	risk.
*/


#ifndef CCONFIG_H
#define CCONFIG_H
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include"exception.h"


using namespace std;
const string comm="#";
class CParamBase {
	public:
		typedef enum {Normal, Default} out_type;
		
		CParamBase(string _name):name(_name){}
		virtual ~CParamBase();
		const string name;
		virtual void parse(istream &in=cin)=0;
		virtual void print(ostream &out=cout, CParamBase::out_type def=Normal)const=0;
	
 	private:
	};


template <class T>
class CParam : public CParamBase {
	public:
		CParam(string _name, T _value):CParamBase(_name), value(_value), def_value(_value){};
		T get() const {
			return value;
			};
		T get_default() const {
			return def_value;
			};
		void set(T _value){
			value=_value;
			}
		void parse(istream &in){
			in>>value;//param.name<<":   "<<param.value;
			}
		void print(ostream &out, CParamBase::out_type def=Normal)const{
			if(def==CParamBase::Default) out<<def_value;

			else out<<value;
			}
		template <class U>
		friend istream &operator>>(istream &in, CParam<U> &param);

		template <class U>
		friend ostream &operator<<(ostream &out, CParam<U> &param);
 	private:
		CParam(){};
		T value;
		const T def_value;
	};

template <class T>
istream &operator>>(istream &in, CParam<T> &param){
			in>>param.value;//param.name<<":   "<<param.value;
			return in;
			}
template <class T>
ostream &operator<<(ostream &out, CParam<T> &param){
			out<<param.name<<":   "<<param.value;
			return out;
			}


// the class CBaseConfig is a singleton 
class CBaseConfig : public CParamBase {
	public:
	CBaseConfig(const string name="main_config"):CParamBase(name){
		}; 
	virtual ~CBaseConfig();                                 
	//static CBaseConfig& Instance();
	bool isValidParam(string fname)const;
	bool isValidSection(string fname)const;
	virtual void parse(istream &input){
		};
	void print(ostream &out=cout, CParamBase::out_type def=Normal)const;

	virtual void define_parameters(){};
	void validate_param(string name)const{
		ERROR( !isValidParam(name), name+" is not a valid parameter or keyword.");
		}
	void validate_section(string name)const{
		ERROR( !isValidSection(name), name+" is not a valid section name.");
		}

	template<class T>
	T get_param(string name, CParamBase::out_type def=CParamBase::Normal)const {
		//const CParamBase  *pp=params[name];
		validate_param(name);
		CParam<T> * const p = static_cast< CParam<T>* > (params.find(name)->second );//i dont know if there is a better solution, i like this tough

		if(def==CParamBase::Default) return p->get_default();
		else return p->get();
		};

	CBaseConfig& get_section(string name){
		validate_section(name);
		return *(sections.find(name)->second);
	}
	template<class T>
	void add_param(string name, T value){
		try{
		CParam<T> *p=new CParam<T>(name,value);
		params[name]=p;
		param_name_dict.push_back(name);
		}catch(...){
			ERROR(1, "Failed creating new paramter.");
			}
		};

	CBaseConfig &add_section(string name){
		try{
		CBaseConfig *c=new CBaseConfig(name);
		sections[name]=c;
		section_name_dict.push_back(name);
		return *(sections.find(name)->second);
		}catch(...){
			ERROR(1, "Failed creating new paramter section.");
			}
		};

	map<string, CParamBase *> params;
	map<string, CBaseConfig *> sections;

	protected:
	vector<string> param_name_dict;
	vector<string> section_name_dict;

	CBaseConfig(const CBaseConfig&);                 // Prevent copy-construction

	private:
	static CBaseConfig *pConfig;
	CBaseConfig& operator=(const CBaseConfig&);      // Prevent assignment
};

#define  paramsDouble config.get_param<double> // just some abbreviations 
#define  paramsStr config.get_param<string> 

//---------------------------- definitions ------------------------------------------
CParamBase::~CParamBase(){}
CBaseConfig::~CBaseConfig(){
	map<string, CParamBase *>::iterator it;
	for(it=params.begin(); it!=params.end(); ++it) delete (*it).second;
	}                                 


/*
CBaseConfig * CBaseConfig::pConfig=NULL;
CBaseConfig &config=CBaseConfig::Instance();
CBaseConfig& CBaseConfig::Instance() {
	static CBaseConfig config;
	return config;
	}
*/


void CBaseConfig::print(ostream &out, CParamBase::out_type def)const{
	map<string, CParamBase *>::const_iterator it;
	for(it=params.begin(); it != params.end(); it++){
		out<<it->first<<"  ";
		it->second->print(out);
		out<<endl;
		}
}



bool CBaseConfig::isValidParam(string vname) const{
	vector<string>::const_iterator it;
	for( it=param_name_dict.begin(); it!=param_name_dict.end(); ++it){
		if(vname==(*it))return true;
		}
	return false;
	}

bool CBaseConfig::isValidSection(string vname) const{
	vector<string>::const_iterator it;
	for( it=section_name_dict.begin(); it!=section_name_dict.end(); ++it){
		if(vname==(*it))return true;
		}
	return false;
	}


#endif
