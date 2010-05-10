#ifndef COMPOSITE_H
#define COMPOSITE_H 
#include"geombase.h"
template<>
class GeomObject<tcomposite>: public GeomObjectBase{
	public:	
	GeomObject<tcomposite> (const GeomObject<tcomposite> & p):GeomObjectBase(p.Xc,tcomposite){
		for(indexType i=0; i<p.elems.size(); i++){
			elems.push_back(new CSphere(*(p.elems.at(i))));
			}
		radius=p.radius;
		Xc=p.Xc;
		}
		
	GeomObject<tcomposite> (const vec &v, double r):GeomObjectBase(v,tcomposite){

		radius=r;
		CSphere *s1=NULL;
		CSphere *s2=NULL, *s3;
		s1=new CSphere(vec(-2*r/3,0.0,0.0), r/3);
		s1->identifier=1;
		s2=new CSphere(vec(0.0), 2*r/3);
		s2->identifier=2;
		s3=new CSphere(vec(2*r/3, 0.0, 0.0), r/3);
		s3->identifier=1;

		ERROR(s1==NULL || s2==NULL, "Improper memory allocation");

		elems.push_back(s1);
		elems.push_back(s2);
		elems.push_back(s3);

		s3=new CSphere(vec(0.0,2*r/3, 0.0), r/3);
		s3->identifier=1;
		elems.push_back(s3);

		s3=new CSphere(vec(0.0,-2*r/3, 0.0), r/3);
		s3->identifier=1;
		elems.push_back(s3);

		s3=new CSphere(vec(0.0,0.0,2*r/3), r/3);
		s3->identifier=1;
		elems.push_back(s3);

		s3=new CSphere(vec(0.0,0.0,-2*r/3), r/3);
		s3->identifier=1;
		elems.push_back(s3);
		//elems.push_back(s2);
		//s2=new CSphere(vec(2,-r/2.), r);
		//elems.push_back(s2);
		//s2=new CSphere(vec(1,-r/2.), r);
		//elems.push_back(s2);
		//s2=new CSphere(vec(1,+r/2.), r);
		//elems.push_back(s2);

		Xc=v;
		for(indexType i=0; i<elems.size(); i++){
			elems.at(i)->Xc=Xc+elems.at(i)->Xc0;
			}
		}

	
	~GeomObject<tcomposite>(){
		for(indexType i=0; i<elems.size(); i++){
			delete elems.at(i);
			elems.at(i)=NULL;
			}
		}

	void moveto(const vec &v){
		vec dx = v-Xc;
		for(indexType i=0; i<elems.size(); i++){
			elems.at(i)->Xc+=dx;
			}
		Xc+=dx;
		}

	double vol(){
		double v=0;
		for(indexType i=0; i<elems.size(); i++){
			v+=elems.at(i)->vol();
			}
		return v;
		}

	double I(vec n){//FIXME this works only when the elements are on the axes
		n.normalize();
		double II=0;
		for(indexType i=0; i<elems.size(); i++){
			II+=elems.at(i)->I(n)+elems.at(i)->vol()*(elems.at(i)->Xc0-(elems.at(i)->Xc0*n)*n).abs2();
			}
		return II;
		}

	double min(size_t j){
		double miny=10;
		for(indexType i=0; i<elems.size(); i++){
			if((elems.at(i)->Xc)(j)-elems.at(i)->radius < miny) miny=((elems.at(i)->Xc)(j)-elems.at(i)->radius);
			}
		return miny;
		};
	void rotateTo(const Quaternion &q){
		for(indexType i=0; i<elems.size(); i++){
			elems.at(i)->Xc=Xc+q.rotate(elems.at(i)->Xc0);
			}
		};

	void rotate(const vec& n , double alpha){//FIXME
		ERROR(true, "check this");
		//q.setRotation(n, alpha);
		for(indexType i=0; i<elems.size(); i++){
			//elems.at(i)->Xc0=q.rotate(elems.at(i)->Xc0);
			}
		};

	void shift(const vec& v){Xc+=v;};
	void scale(double scale){};//FIXME

	void print(std::ostream &out)const{
		for(indexType i=0; i<elems.size()-1; i++){
			elems.at(i)->print(out);
			out<< endl;
			}
			elems.back()->print(out);
		}
	void parse(std::istream &in){ERROR(true, "check this.");}
	
	vector<CSphere *> elems;
	private:
	GeomObject<tcomposite> ();
	};

#endif /* COMPOSITE_H */
