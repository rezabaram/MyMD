// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


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
