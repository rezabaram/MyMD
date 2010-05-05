#ifndef DFREEDOM_H
#define DFREEDOM_H 

template<indexType depth=5>//number of derivatives, i.e., r, v, a, j, etc ...
class CDFreedom {
	public:
	vec &operator()(indexType i){return x[i];}
	const vec &operator()(indexType i)const{return x[i];}
	template<short order>
	void gear_predict(double dt){
		if(depth<order){
			ERROR("Number of derivatives cannot be less than order of gear integration.");
			exit(1);
			}
		double ddt=1.0;
		for(int i=0; i<order-1; i++){
			ddt=1.0;
			double fac=0;
			for(int j=i+1; j<order; j++){//tailor exapansion
				fac++;
				ddt*=dt/fac;
				x[i]+=ddt*x[j];
				//ERROR(i<<"  "<<j<<"  "<<ddt)
				}
			}
		}

	template<short order>
	void gear_correct(double dt, const vec &dA){
		if(order<4 || order>6){
			ERROR("Only gear schemes of order 4, 5, and 6 are implemented.");
			exit(1);
			}
		if(depth<order){
			ERROR("Number of derivatives cannot be less than order of gear integration.");
			exit(1);
			}
		static const double c4[] ={ 1.0/6.0,  5.0/6.0,  1.0, 1.0/3.0};//forth order
		static const double c5[] ={ 19.0/90.0,  3.0/4.0,  1.0, 1.0/2.0, 1./12.0};//fifth order
		static const double c6[] ={ 3.0/16.0,  251.0/360.0,  1.0, 11.0/18.0, 1./6.0, 1./60.0};//6th order
	

		double ddt=dt*dt/2.0;
		for(int i=0; i<order; i++){
			if(order==4)x[i]+=dA*(ddt*c4[i]);
			if(order==5)x[i]+=dA*(ddt*c5[i]);
			if(order==6)x[i]+=dA*(ddt*c6[i]);
			ddt/=dt;
			ddt*=(double)(i+1);
			}
		}
 	private:
	vec x[depth];
	};

#endif /* DFREEDOM_H */
