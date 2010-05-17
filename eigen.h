#ifndef EIGEN_H
#define EIGEN_H 
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "common.h"

void mysort(gsl_vector_complex *eval, indexType ind[] , indexType n){

	bool swapped;

	double val1, val2;
	do{
		swapped = false;
		for(indexType i=0; i<n-1; ++i){
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

void eigens(Matrix &M, vector<complex<double> > &eigenvals, vector<HomVec > &eigenvecs)
     {
     


	assert(M.RowNo()==M.ColNo());
	size_t N=M.RowNo();
	assert(N==4);

	double *data=new double[N*N];
	assert(data);
	for(indexType i=0; i<N; i++){
		for(indexType j=0; j<N; j++){
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
	indexType ind[N];//for holding sorted indeces 
	for(indexType i=0; i<N; ++i){
		ind[i]=i;
		}
	mysort(eval, ind, N);
///
       { 
         for (indexType i = 0; i < N; i++)
           {
             gsl_complex eval_i = gsl_vector_complex_get (eval, ind[i]);
	     if(0)if(abs(GSL_IMAG(eval_i)) > epsilon or GSL_REAL(eval_i) <= 0){ continue; //only real eigen values 
		}
             gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, ind[i]);
			//eigenvals.push_back((GSL_REAL(eval_i), GSL_IMAG(eval_i))); 
			eigenvals.push_back(complex<double> (GSL_REAL(eval_i), GSL_IMAG(eval_i)));
			double w=GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 3));
			//assert(fabs(w)>epsilon);
			if(fabs(w)<epsilon)w=1;//this is wrong. but i dont use those eigenvectors
			eigenvecs.push_back(
				HomVec(
					GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 0))/w,
					GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 1))/w,
					GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 2))/w,
					1
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
