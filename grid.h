#ifndef GRID_H
#define GRID_H
#include"common.h"
#include"particle.h"
#include"box.h"

enum {
	Periodic=1,
};

enum EDirction
	{ left, right,
	back, forth,
	up, down
	};

//give names to the neighbours not to have to deal with their indices
enum {
	self=0, neigh_N,
	neigh_S, neigh_E,
	neigh_W, neigh_NE,
	neigh_NW, neigh_SE,
	neigh_SW, 
	neigh_U,
	neigh_UN, neigh_US,
	neigh_UE, neigh_UW,
	neigh_UNE, neigh_UNW,
	neigh_USE, neigh_USW,
	neigh_D,
	neigh_DN, neigh_DS,
	neigh_DE, neigh_DW,
	neigh_DNE, neigh_DNW,
	neigh_DSE, neigh_DSW,
};


//an ad hoc class for int coordinates. Maybe I will make vec3d a template to simply include this class in it.
class CCoord : public CVector<3, int>
	{
	public:
	CCoord(){};
	CCoord(const int i, const int j, const int k)
		:CVector<3,int>()
		{
		(*this)(0)=i; (*this)(1)=j; (*this)(2)=k;
		}
	};


//a class to define an interval (window) in 3D
class CWindow 
	{
      public:
	CWindow(){};
	explicit CWindow(const CCoord &_lleft, const CCoord &_uright):lleft(_lleft),uright(_uright){ }
	explicit CWindow(const int i1, const int j1, const int k1,  const int i2, const int j2, const int k2):lleft(CCoord(i1,j1, k1)),uright(CCoord(i2,j2, k2))
		{
		assert(i2>=i1 && j2>=j1 && k2>=k1);
		 }

	void trim(const CWindow &win) //something like croping
		{
		//to warn about unhandled cases
		//assert(win.uright(0)(0) >= lleft.x(0) ) ;
		//assert(win.lleft(0)(0) <= uright.x(0) ) ;
		//assert(win.lleft(0)(1) <= uright.x(1) ) ;
		//assert(win.uright(0)(1) >= lleft.x(1) ) ;

		if(win.lleft(0) > lleft(0) ) lleft(0)=win.lleft(0);
		if(win.lleft(1) > lleft(1) ) lleft(1)=win.lleft(1);
		if(win.lleft(2) > lleft(2) ) lleft(2)=win.lleft(2);
		if(win.uright(0) < uright(0) ) uright(0)=win.uright(0);
		if(win.uright(1) < uright(1) ) uright(1)=win.uright(1);
		if(win.uright(2) < uright(2) ) uright(2)=win.uright(2);
		//assert(uright(0)(0)>=lleft.x(0) && uright.x(1)>=lleft.x(1));
		}
	CCoord lleft, uright;//the name is misleading, these are actually two opposite corners of a cube
	private:
};

typedef std::vector<CParticle *> Container;

//this class is for nodes which keeps the information of their neighborhood 
//the node is unaware where in the grid it is. this information is kept in grid class (CRecGrid).
class CNode3D : public Container 
	{
	public:
	CNode3D():neighs(std::vector<CNode3D*>(27,(CNode3D *)NULL)){}
	~CNode3D(){};
	
	std::vector<CNode3D*> neighs;
	//vector<CNode3D*> extend_neighs;
	private:
	CNode3D(const CNode3D &);
};

//and the Grid class
class CRecGrid
	{
	public:
	CRecGrid(const vec &_corner, const vec & _L, double _d);
	~CRecGrid();
	CNode3D &operator ()(int i, int j, int k); //returns the ith node
	CNode3D *node(int i, int j, int k); //returns the ith node

	CNode3D* top_filled_node(double x, double y);
	inline CNode3D* top_node(int i, int y);
	void add(CParticle *part);
	void add(CParticle *part, const CWindow &win);
	CNode3D* which(const vec3d &p) {
		return node (floor((p-corner)(0)/dL(0)), floor((p-corner)(1)/dL(1)), floor((p-corner)(2)/dL(2)));
		}

	void print(std::ostream &out=std::cout){
		out<< "print not yet implemented for this class" <<std::endl;
        	}

	vec corner;
	vec L, dL;//dimensions of the grid L, and cells dL 
	CCoord N;//number of cells in each direction

	//private:
	CWindow window;
	void setup_neighs();
	CNode3D *nodes;
	int *top_nodes;//keeps the k's of the topest node
};


CRecGrid::CRecGrid(const vec & _corner, const vec & _L, double _d):
	corner(_corner), L(_L)
	{
	
	if(_d>L(0))_d=L(0);
	if(_d>L(1))_d=L(1);
	if(_d>L(2))_d=L(2);

	N(0)=(unsigned int)(L(0)/_d+1e-12);
	N(1)=(unsigned int)(L(1)/_d+1e-12);
	N(2)=(unsigned int)(L(2)/_d+1e-12);
	dL(0)=L(0)/(double)(N(0));
	dL(1)=L(1)/(double)(N(1));
	dL(2)=L(1)/(double)(N(2));

	nodes=(CNode3D *)NULL;
	top_nodes=(int *)NULL;
	if(N(0)==0||N(1)==0){ERROR("No grid will be created.") return;}
	nodes=new CNode3D[N(0)*N(1)*N(2)];
	assert(nodes);
	top_nodes=new int[N(0)*N(1)];
	if(nodes==NULL || top_nodes==NULL){ERROR("No grid will be created.") return;}

	setup_neighs();
	for(int i=0; i<N(0)*N(1); i++){
                top_nodes[i]=NULL;
                }
	window=CWindow(0,0,0, N(0)-1,N(1)-1, N(2)-1);
}

CRecGrid::~CRecGrid()
	{
	if(nodes!=NULL) delete [] nodes;
	if(top_nodes!=NULL) delete [] top_nodes;
	}

CNode3D & CRecGrid::operator ()(int i, int j, int k)
	{
	return *node(i,j, k);
	}

void CRecGrid::add(CParticle *cir,const CWindow &win)
	{
	assert(cir);
	for(int i=win.lleft(0); i<=win.uright(0); i++){
	for(int j=win.lleft(1); j<=win.uright(1); j++){
	for(int k=win.lleft(2); k<=win.uright(2); k++){
		node(i,j,k)->push_back(cir);
		//if(k>top_nodes[j*N(0)+i])top_nodes[j*N(0)+i]=k;
		}
		}
		}
	}

void CRecGrid::add(CParticle *part)
	{
	assert(part);
	vec3d R=part->x(0)-corner;
	double r=part->size;
	CWindow win(floor((R(0)-r)/dL(0)), floor((R(1)-r)/dL(1)), floor((R(2)-r)/dL(2)),
	            floor((R(0)+r)/dL(0)), floor((R(1)+r)/dL(1)), floor((R(2)+r)/dL(2))
		   );
	win.trim(window);
	add(part,win);
	}

inline CNode3D * CRecGrid::node(int i, int j, int k)
	{
	if(0  /*flags&Periodic*/){//FIX ME
		while (i<0)i+=N(0);
		while (i>N(0)-1)i-=N(0);
		
		while (j<0)j+=N(1);
		while (j>N(1)-1)j-=N(1);
		
		}
	else {
		if( i<0)return NULL;
		if( i>N(0)-1)return NULL;

		if( j<0)return NULL;
		if( j>N(1)-1)return NULL;

		}

//no periodic in z direction
	if( k<0)return NULL;
	if( k>N(2)-1)return NULL;

//origin is considered at lower left
	assert(nodes);
	return (&nodes[k*(N(0)*N(1))+j*N(0)+i]);
}


void CRecGrid::setup_neighs()
	{
	for(int i=0; i<N(0); i++){
		for(int j=0; j<N(1); j++){
		for(int k=0; k<N(2); k++){
			//if(node(i,j)==NULL)continue;
			node(i,j,k)->neighs.at(self)=node(i,j, k);
			node(i,j,k)->neighs.at(neigh_N)=node(i,j+1, k);
			node(i,j,k)->neighs.at(neigh_S)=node(i,j-1, k);
			node(i,j,k)->neighs.at(neigh_E)=node(i+1,j, k);
			node(i,j,k)->neighs.at(neigh_W)=node(i-1,j, k);
			node(i,j,k)->neighs.at(neigh_NE)=node(i+1,j+1, k);
			node(i,j,k)->neighs.at(neigh_NW)=node(i-1,j+1, k);
			node(i,j,k)->neighs.at(neigh_SE)=node(i+1,j-1, k);
			node(i,j,k)->neighs.at(neigh_SW)=node(i-1,j-1, k);

			node(i,j,k)->neighs.at(neigh_U)=node(i,j, k+1);
			node(i,j,k)->neighs.at(neigh_UN)=node(i,j+1, k+1);
			node(i,j,k)->neighs.at(neigh_US)=node(i,j-1, k+1);
			node(i,j,k)->neighs.at(neigh_UE)=node(i+1,j, k+1);
			node(i,j,k)->neighs.at(neigh_UW)=node(i-1,j, k+1);
			node(i,j,k)->neighs.at(neigh_UNE)=node(i+1,j+1, k+1);
			node(i,j,k)->neighs.at(neigh_UNW)=node(i-1,j+1, k+1);
			node(i,j,k)->neighs.at(neigh_USE)=node(i+1,j-1, k+1);
			node(i,j,k)->neighs.at(neigh_USW)=node(i-1,j-1, k+1);

			node(i,j,k)->neighs.at(neigh_D)=node(i,j, k-1);
			node(i,j,k)->neighs.at(neigh_DN)=node(i,j+1, k-1);
			node(i,j,k)->neighs.at(neigh_DS)=node(i,j-1, k-1);
			node(i,j,k)->neighs.at(neigh_DE)=node(i+1,j, k-1);
			node(i,j,k)->neighs.at(neigh_DW)=node(i-1,j, k-1);
			node(i,j,k)->neighs.at(neigh_DNE)=node(i+1,j+1, k-1);
			node(i,j,k)->neighs.at(neigh_DNW)=node(i-1,j+1, k-1);
			node(i,j,k)->neighs.at(neigh_DSE)=node(i+1,j-1, k-1);
			node(i,j,k)->neighs.at(neigh_DSW)=node(i-1,j-1, k-1);

			}
		}
		}
	}

CNode3D * CRecGrid::top_node(int i, int j)
        {
        assert(top_nodes);
        return node(i, j, top_nodes[j*N(0)+i]);
        }


CNode3D* CRecGrid::top_filled_node(double x, double y)
	{
	return top_node(floor((x-corner(0))/dL(0)), floor((y-corner(1))/dL(1)));
	}
	
#endif
