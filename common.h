#ifndef COMMON_H
#define COMMON_H 

#define ERROR(x)  std::cerr<<"Error: In file " __FILE__<<" line "<<__LINE__<<":  "<<x<<std::endl;
//#include"matrix.h"
#include<assert.h>
#include<memory>
#include<limits>
#include<sstream>
#include<ostream>
#include<vector>
#include"vector.h"
#include"log.h"

#include"config.h"
extern CConfig &config; // don't forget "&" or you get a vicious bug, which took me one day to find

typedef vec3d vec;

double aG[]={0.0, 0, -10.0};
vec G(aG);

void define_parameters()
{
	config.add_param<vec>("Gravity", G);
	config.add_param<double>("outDt", 0.02);
	config.add_param<double>("stiffness1", 5000000.0); 
	config.add_param<double>("stiffness2", 1000000.0); 
	config.add_param<double>("density", 10000.0); 
	config.add_param<double>("particleSize", 0.05); 
}
#endif /* COMMON_H */
