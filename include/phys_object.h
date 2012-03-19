// This file is a part of Molecular Dynamics code for 
// simulating ellipsoidal packing. The author cannot 
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram 


#ifndef PHYS_OBJECT_H
#define PHYS_OBJECT_H 

class PhysObject
	{
	public:
	PhysObject():frozen(false){};

	bool frozen;
 	private:
	};
#endif /* PHYS_OBJECT_H */
