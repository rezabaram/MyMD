#ifndef BOX_H
#define BOX_H 
#include<limits>
#include"geombase.h"
#include"phys_object.h"

typedef enum{wall, soft, periodic} BoundaryType;

class CBox: public GeomObjectBase
	{
	public:
	CBox(vec corner=vec(std::numeric_limits<double>::max()), vec _L=vec(0.0), string _btype="wall"):
		GeomObjectBase(corner+_L/0.5, tbox), corner(corner), L(_L), nFaces(6),
		btype(_btype),
		u0(vec(1.0,0.0,0.0)), u1(vec(0.0,1.0,0.0)), u2(vec(0.0,0.0,1.0))
		 {
		identifier=6;
		face =  new CPlane * [nFaces];
		face[0]=new CPlane (corner,u0);
		face[1]=new CPlane (corner,u1);
		face[2]=new CPlane (corner,u2);

		face[3]=new CPlane (corner+L,-u0);
		face[4]=new CPlane (corner+L,-u1);


		if(btype=="periodic_x")
			{
			face[0]->solid=false;
			face[3]->solid=false;
			}
		else if(btype=="periodic_xy")
			{
			face[0]->solid=false;
			face[3]->solid=false;
			face[1]->solid=false;
			face[4]->solid=false;
			}
		else if(btype=="periodic_xyz")
			{
			face[0]->solid=false;
			face[3]->solid=false;
			face[1]->solid=false;
			face[4]->solid=false;
			face[2]->solid=false;
			face[5]=new CPlane (corner+L,-u2);
			face[5]->solid=false;
			}
		else if(btype=="wall"){}
		else{ERROR(1, "Boundary condition not defined");}
		};
	virtual ~CBox(){
		for(int i=0; i<5; i++){
			delete face[i];
			}
		if(btype=="periodic_xyz")delete face[5];
		delete [] face;
		}

	void shift(const vec &x){corner+=x;}
	void rotate(const vec &n, double alpha){}//FIXME 
	void scale(double scale){L*=scale;};
	void print(std::ostream &out)const{//FIXME this is ad-hoc
		out<< identifier<< "   "<<corner<<"  "<<corner+vec(1,0,0)<<"  "<<corner+vec(1,1,0)<<endl;
		out<< identifier<< "   "<<corner<<"  "<<corner+vec(0,0,1)<<"  "<<corner+vec(0,1,1)<<endl;
		out<< identifier<< "   "<<corner+L<<"  "<<corner+L+vec(0,0,1)<<"  "<<corner+vec(0,1,1)<<endl;
		out<< identifier<< "   "<<corner+L<<"  "<<corner+L+vec(0,0,1)<<"  "<<corner+vec(1,0,1)<<endl;
		out<< identifier<< "   "<<corner+L<<"  "<<corner+L+vec(1,0,0)<<"  "<<corner+L+vec(1,1,0)<<endl;
		}

	double vol()const{return L(0)*L(1)*L(2);}
	double I(vec n){
		ERROR(true, "Not implemented."); //FIXME
		return 0;}

	void parse(std::istream &in){
		in>>identifier;
		in>>corner>>L;
		}

	void scale(double s0, double s1, double s2){L(0)*=s0; L(1)*=s1; L(2)*=s2;}

	vec top()const{return corner+L;}
	vec corner, L;
	const size_t nFaces;
	CPlane ** face;
	string btype;
	protected:
	vec u0, u1, u2;
 	private:
	CBox();
	};

class BoxContainer : public CBox, public PhysObject
	{
	public:
	BoxContainer(vec corner=vec(std::numeric_limits<double>::max()), vec _L=vec(0.0), string _btype="wall"):
	CBox(corner, _L, _btype)
		{ 
		
		//assert(corner<L); FIXME make this work
		}	

 	private:
	};


#endif /* BOX_H */
