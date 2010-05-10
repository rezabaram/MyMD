#ifndef EXCEPTION_H
#define EXCEPTION_H 
#include<iostream>
#include<sstream>
using namespace std;

class CException
	{
	string message;
	string file;
	string function;
	long line;
	CException(){};
	public:
	CException(const string &m="Unknown", const string &f="UnknownFile", const string &func="UnknownFunction", long i=0):message(m), file(f), function(func), line(i) {};
	void Report(){
		cerr<<"Error: "<<file<<":"<<line<<":"<<function<<": "<< message <<endl;
		}
	};

//ERROR is defined for convenience 
#define ERROR(b, m)  if(b){throw CException(m, __FILE__, __PRETTY_FUNCTION__, __LINE__);}
#define WARNING(x)  std::cerr<<"Warning: "<<__FILE__<<":"<<__LINE__<<": "<<x<<std::endl;
#endif /* EXCEPTION_H */
