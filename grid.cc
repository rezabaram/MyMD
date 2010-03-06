#include "grid.h"
std::ostream& operator<< (std::ostream &os, const CCoord &a)
	{
	os<<a.print_str()<<endl;
	return os;
	} ;
