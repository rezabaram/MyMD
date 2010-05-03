#ifndef EIGEN_H
#define EIGEN_H 
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "common.h"

void mysort(gsl_vector_complex *eval, int ind[] , int n){

	bool swapped;

	double val1, val2;
	do{
		swapped = false;
		for(int i=0; i<n-1; ++i){
                         val1=GSL_REAL(gsl_vector_complex_get (eval, ind[i]));
                         val2=GSL_REAL(gsl_vector_complex_get (eval, ind[i+1]));
		      if (val1> val2){
			swap( ind[i], ind[i+1] );
			swapped = true;
			}
			}
		--n;
	}while (swapped);


}

void eigens(Matrix &M, vector<double> &eigenvals, vector<vec3d<double> > &eigenvecs)
     {
     


	assert(M.RowNo()==M.ColNo());
	size_t N=M.RowNo();
	assert(N==4);

	double *data=new double[N*N];
	assert(data);
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			data[i*N+j]=M(i,j);
			}
		}
	//vector< complex >


       gsl_matrix_view m = gsl_matrix_view_array (data, N, N);
       gsl_vector_complex *eval = gsl_vector_complex_alloc (N);
       gsl_matrix_complex *evec = gsl_matrix_complex_alloc (N, N);
     
       gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (N);
       gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
     
       gsl_eigen_nonsymmv_free (w);
       //gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	

///     sorting the eigenvalues 
	int ind[N];//for holding sorted indeces 
	for(int i=0; i<N; ++i){
		ind[i]=i;
		}
	mysort(eval, ind, N);
///
       
       { 
         for (int i = 0; i < N; i++)
           {
	     
             gsl_complex eval_i = gsl_vector_complex_get (eval, ind[i]);
	     if(abs(GSL_IMAG(eval_i)) > epsilon) continue; //only real eigen values 
             gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, ind[i]);
			eigenvals.push_back((GSL_REAL(eval_i), GSL_IMAG(eval_i))); 
			double w=-GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 3));
			assert(fabs(w)>epsilon);
			eigenvecs.push_back(
				vec3d<double > (
					GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 0))/w,
					GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 1))/w,
					GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 2))/w
					)
				); 

			//eigenvecs.push_back(vec3d<complex<double> > (
			 	//complex<double> (GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 0)), GSL_IMAG(gsl_vector_complex_get (&evec_i.vector, 0))), 
			 	//complex<double> (GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 1)), GSL_IMAG(gsl_vector_complex_get (&evec_i.vector, 1))), 
			 	//complex<double> (GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 2)), GSL_IMAG(gsl_vector_complex_get (&evec_i.vector, 2))))
				//); 
     
             //printf ("\neigenvalue = %g\n", GSL_REAL(eval_i));
             //printf ("eigenvector = \n");
             //gsl_vector_complex_fprintf (stdout, &evec_i.vector, "%g");
           }
       }
     
       gsl_vector_complex_free (eval);
       gsl_matrix_complex_free (evec);

	
	delete [] data;
     }

#endif /* EIGEN_H */
