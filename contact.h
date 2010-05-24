#ifndef CONTACT_H
#define CONTACT_H 
#include"vec.h"
#include"plane.h"
#include<iostream>

class Contact{
	public:
	Contact(const vec &_x, const vec &_dx , bool pers=false):x(_x), dx(_dx), persist(pers){
		}
	vec x, dx;
	bool persist;
	};

#endif /* CONTACT_H */
