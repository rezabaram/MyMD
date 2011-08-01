#ifndef EIGEN_H
#define EIGEN_H 
#include <stdio.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "matrix.h"
#include "vec.h"

using namespace std;

//FIXME not goood. put things in a better place
using namespace math;
typedef matrix<double> Matrix;
//////////

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

	static double data[4*4];
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
         for (indexType i = 0; i < N; i++) {
		gsl_complex eval_i = gsl_vector_complex_get (eval, ind[i]);
		gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, ind[i]);
		//eigenvals.push_back((GSL_REAL(eval_i), GSL_IMAG(eval_i))); 

		eigenvals.push_back(complex<double> (GSL_REAL(eval_i), GSL_IMAG(eval_i)));
		eigenvecs.push_back(
			//convert to my format
			HomVec(
				GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 0)),
				GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 1)),
				GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 2)),
				GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 3))
				)
			); 
           	}
     
       gsl_vector_complex_free (eval);
       gsl_matrix_complex_free (evec);
	
     }
//for 3x3 matrix
void eigens(const Matrix &M, vector<double> &eigenvals, vector<vec3d> &eigenvecs)
     {

    	assert(M.RowNo()==M.ColNo());
	size_t N=M.RowNo();
	assert(N==3);

	static double data[3*3];
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
	//mysort(eval, ind, N);
///
         for (indexType i = 0; i < N; i++) {
		gsl_complex eval_i = gsl_vector_complex_get (eval, ind[i]);
		gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, ind[i]);
		//eigenvals.push_back((GSL_REAL(eval_i), GSL_IMAG(eval_i))); 

		eigenvals.push_back(GSL_REAL(eval_i));
		eigenvecs.push_back(
			//convert to my format
			vec3d(
				GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 0)),
				GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 1)),
				GSL_REAL(gsl_vector_complex_get (&evec_i.vector, 2))
				)
			); 
           	}
     
       gsl_vector_complex_free (eval);
       gsl_matrix_complex_free (evec);
	
     }
#endif /* EIGEN_H */
