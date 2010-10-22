#ifndef CONTACT_H
#define CONTACT_H 
#include<iostream>
#include"vec.h"
#include"plane.h"

class BasicContact{
	public:
	BasicContact(const vec &_x, const vec &_n):x(_x), n(_n){
		}
	vec x, n;
	};

class Contact : public BasicContact{
	public:
	Contact(const vec &_x, const vec &_n , double dd, bool pers=false):BasicContact(_x, _n), dx_n(dd), persist(pers){
		}
	double dx_n;
	bool persist;
	};

#endif /* CONTACT_H */
