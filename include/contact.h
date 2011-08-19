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
	Contact(const vec &_x, const vec &_n , double dd, const void * _p1, const void * _p2):
		BasicContact(_x, _n), dx_n(dd), p1(_p1), p2(_p2),
		static_friction_on(false)
		{
		
		}
	double dx_n;
	const void *p1, *p2;
	vec x1, x2;//used to implement static friction
	bool static_friction_on;
	};

#endif /* CONTACT_H */
