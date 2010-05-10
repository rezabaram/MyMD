#ifndef EXCEPTION_H
#define EXCEPTION_H 
#include<iostream>
#include<sstream>

class CException
	{
	std::string message;
	std::string file;
	std::string function;
	long line;
	CException(){};
	public:
	CException(const std::string &m="Unknown", const std::string &f="UnknownFile", const std::string &func="UnknownFunction", long i=0):message(m), file(f), function(func), line(i) {};
	void Report(){
		std::cerr<<"Error: "<<file<<":"<<line<<": "<<function<<": "<< message <<std::endl;
		}
	};

//ERROR is defined for convenience 
#define ERROR(b, m)  if(b){throw CException(m, __FILE__, __PRETTY_FUNCTION__, __LINE__);}
#define WARNING(x)  std::cerr<<"Warning: "<<__FILE__<<":"<<__LINE__<<": "<<x<<std::endl;
#endif /* EXCEPTION_H */
