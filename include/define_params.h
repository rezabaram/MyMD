#ifndef DEFINE_PARAMS_H
#define DEFINE_PARAMS_H 
#include"params.h"
#include"vec.h"
#include<string>
using namespace std;

void define_parameters()
{
       config.add_param<vec>("Gravity", vec(0.0, 0.0, -10.0));

	//controlling the output
       config.add_param<double>("outStart", 0.00);
       config.add_param<double>("outEnd", 1000.00);
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
       config.add_param<string>("boundary", "wall"); 
       config.add_param<string>("radii", "radii.dat"); 
}

#endif /* DEFINE_PARAMS_H */
