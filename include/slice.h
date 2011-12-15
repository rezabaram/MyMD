#ifndef SLICE_H
#define SLICE_H 
#include<iostream>
#include<fstream>
#include<string>

#include"exception.h"
using namespace std;

class CPixel
{
	public:
	CPixel():totalvalue(0),norm_factor(1){}
	void reset(){totalvalue=0;norm_factor=1;}
	float value(){return totalvalue/norm_factor;}
	void add(float a){totalvalue+=a; norm_factor++;}
 	private:
	float totalvalue;
	float norm_factor;
};

class CSlice
{
	public:
	~CSlice(){
		delete [] grid;
	}
	CSlice(const vec2d &_corner=vec2d(0.0,0.0), const vec2d & _L=vec2d(1.0,1.0), double _d=0.002): corner(_corner), L(_L){
		ERROR(_d<0.0000001,"Invalid grid size either negative or too small.");
		if(_d>L(0))_d=L(0);
		if(_d>L(1))_d=L(1);

		nx=(size_t)(L(0)/_d+1e-12);
		ny=(size_t)(L(1)/_d+1e-12);
		if(nx<1)nx=1;
		if(ny<1)ny=1;

		dx=L(0)/(double)(nx);
		dy=L(1)/(double)(ny);
		allocate();
	}

	void reset(){
		for (unsigned int i = 0;i < nx*ny;i++) 
			grid[i].reset(); 
	}
	void allocate(){
		grid=new CPixel[nx*ny];
	}
	void add(const vec2d &p, float a){
		size_t i= (size_t)floor((p-corner)(0)/dx);
		size_t j= (size_t)floor((p-corner)(1)/dy);
		if(i>=nx or j>=ny)return;
		(*this)(i,j).add(a);
	}
	
	CPixel &which(const vec2d &p){
		size_t i= (size_t)floor((p-corner)(0)/dx);
		size_t j= (size_t)floor((p-corner)(1)/dy);
		return (*this)(i,j);
	}

	CPixel &operator()(size_t i, size_t j){
		ERROR(i>=nx or j>=ny, "Index out of bound");
		return grid[i*ny+j];
	}
	void export_ppm(string filename){
		ofstream out(filename.c_str());
		out<<"P3"<<endl;
		out<<"# created by Reza Baram"<<endl;
		out<< nx <<"\t"<< ny<<endl;
		out<<"255" <<endl;
		for (unsigned int j = 0;j < ny;j++) { 
			for (unsigned int i = 0;i < nx;i++) { 
				double value=(*this)(nx-i-1,ny-j-1).value();
				out<< (int)(value*255) <<" "<< (int)(value*255) <<" "<<(int)(value*255)<<"  "; ;
			}
			out<<endl;
		}
	}

	size_t nx, ny;
	float dx, dy;
	vec2d  corner, L;
 	private:
	CPixel *grid;
};

#endif /* SLICE_H */
