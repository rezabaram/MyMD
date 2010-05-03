#ifndef EIGEN_H
#define EIGEN_H 
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "common.h"


 void eigens(Matrix &M)
     {
     
	assert(M.RowNo()==M.ColNo());
	size_t N=M.RowNo();
	double *data=new double(N*N);
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			data[i*N+j]=M(i,j);
			}
		}

       gsl_matrix_view m = gsl_matrix_view_array (data, N, N);
       gsl_vector_complex *eval = gsl_vector_complex_alloc (N);
       gsl_matrix_complex *evec = gsl_matrix_complex_alloc (N, N);
     
       gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (N);
       gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
     
       gsl_eigen_nonsymmv_free (w);
       gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
       
       { int i;
     
         for (i = 0; i < N; i++)
           {
             gsl_complex eval_i 
                = gsl_vector_complex_get (eval, i);
             gsl_vector_complex_view evec_i 
                = gsl_matrix_complex_column (evec, i);
     
             printf ("eigenvalue = %g\n", GSL_REAL(eval_i));
             printf ("eigenvector = \n");
             gsl_vector_complex_fprintf (stdout, &evec_i.vector, "%g");
           }
       }
     
       gsl_vector_complex_free (eval);
       gsl_matrix_complex_free (evec);
     }

#endif /* EIGEN_H */
