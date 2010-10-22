#ifndef MAIN_H
#define MAIN_H 

#include"common.h"

extern CConfig &config; // don't forget "&" or you get a vicious bug, which took me one day to find

void define_parameters()
{
	config.add_param<vec>("Gravity", G);
	config.add_param<double>("outDt", 0.02);
	config.add_param<double>("stiffness", 5.0e+05); 
	config.add_param<double>("damping", 1); 
	config.add_param<double>("fluiddampping", 0.05); 
	config.add_param<double>("friction", 0); 
	config.add_param<double>("cohesion", 0); 
	config.add_param<double>("density", 10000.0); 
	config.add_param<double>("particleSize", 0.05); 
	config.add_param<double>("timeStep", 0.00001); 
	config.add_param<double>("maxTime", 4.0); 
	config.add_param<double>("verletfactor", 0.5); 
	config.add_param<size_t>("nParticle", 5); 
	config.add_param<string>("particleShape", "sphere"); 
	config.add_param<double>("e", 0.5); 
	config.add_param<double>("asphericity", -0.5); 
	config.add_param<double>("asphericityWidth", 0.1); 
}
#endif /* MAIN_H */