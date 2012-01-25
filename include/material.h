#ifndef MATERIAL_H
#define MATERIAL_H 

class CMaterial
	{
	public:
	CMaterial(){
		stiffness=5e+5; friction=0;
		static_friction=0; friction_threshold=100;
		cohesion=0;
		density=1;
		color=" 1";
		}
	~CMaterial(){}
	double stiffness, damping, friction, static_friction, friction_threshold, cohesion;
	double density;
	string color;
 	private:
	};
#endif /* MATERIAL_H */
