#ifndef CONTACT_H
#define CONTACT_H 
#include"vec.h"
#include"plane.h"
#include<iostream>

class Contact{
	public:
	Contact(const vec &_x, const vec &_n , double dd, bool pers=false):x(_x), n(_n), dx_n(dd), persist(pers){
		}
	vec x, n;
	double dx_n;
	bool persist;
	};

#endif /* CONTACT_H */
