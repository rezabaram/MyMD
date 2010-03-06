#ifndef BOX_H
#define BOX_H 


class CBox {
	public:
	CBox(vec3d x1=vec3d(std::numeric_limits<double>::max()), vec3d x2=vec3d(0.0)):
		corner(x1), L(x2) {};

	vec3d top()const{return corner+L;}

	vec3d corner, L;
	vector<GeomObject *> elems;
 	private:
	};

#endif /* BOX_H */
