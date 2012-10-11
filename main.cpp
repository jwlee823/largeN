#include <iostream>
#include <gsl/gsl_matrix.h> // if GSL is not installed, you would get error messages
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

int main() {
//int main(int argc, char *argv) {
  int Nrep = 2; //dimension of the representation
  int Ngen = Nrep*Nrep-1; //number of the SU(N) generators
  int i, j, k;
  gsl_matrix_complex * tf[Ngen]; //generators of SU(N) in the fundamental representation
  gsl_matrix_complex * ta[Ngen]; //generators of SU(N) in the adjoint representation
  int tmp = 0;
  gsl_complex z;
  GSL_SET_COMPLEX (&z, 0.0, 1.0); //define the pure imaginary number

//////////////////////////////////////////////////////////////////////////////////////////////////
//                 Fundamental representation of SU(N)
//////////////////////////////////////////////////////////////////////////////////////////////////

//Off-diagonal matrices
    for (i = 0; i < Nrep-1; i++)
      for (j = i+1; j < Nrep; j++) {
        tf[tmp] = gsl_matrix_complex_calloc (Nrep, Nrep);
        gsl_matrix_complex_set (tf[tmp], i, j, gsl_complex_rect (1/M_SQRT2,0.0));
        gsl_matrix_complex_set (tf[tmp], j, i, gsl_complex_rect (1/M_SQRT2,0.0));
	tmp += 1;
      }

    for (i = 0; i < Nrep-1; i++)
      for (j = i+1; j < Nrep; j++) {
        tf[tmp] = gsl_matrix_complex_calloc (Nrep, Nrep);
        gsl_matrix_complex_set (tf[tmp], i, j, gsl_complex_div_real (z, -1.0*M_SQRT2));
        gsl_matrix_complex_set (tf[tmp], j, i, gsl_complex_div_real (z,  M_SQRT2));
	tmp += 1;
      }

//Diagonal matrices
    for (k = 1; k < Nrep; k++) {
      tf[tmp] = gsl_matrix_complex_calloc (Nrep, Nrep);
      for (i = 0; i < k; i++) {
	 gsl_matrix_complex_set (tf[tmp], i, i, gsl_complex_rect (1.0/sqrt (k*(k+1.0)),0.0));
     }
      gsl_matrix_complex_set (tf[tmp], k, k, gsl_complex_rect (-k/sqrt (k*(k+1.0)),0.0));
      tmp += 1;
    }

//Print out the generators of SU(N) in the fundamental rep

  for (k = 0; k < tmp; k++){
    for (i = 0; i < Nrep; i++) 
      for (j = 0; j < Nrep; j++)
        printf ("Tf(%d,%d,%d) = %10g%10g\n",k, i, j, GSL_REAL(gsl_matrix_complex_get (tf[k], i, j)), GSL_IMAG(gsl_matrix_complex_get (tf[k], i,j)));

  }

//////////////////////////////////////////////////////////////////////////////////////////////////
//                 Adjoint representation of SU(N)
//////////////////////////////////////////////////////////////////////////////////////////////////

  gsl_complex alpha, beta;
  GSL_SET_COMPLEX (&alpha, 1.0, 0.0);
  GSL_SET_COMPLEX (&beta, 0.0, 0.0); 
  int a, b, c;
  gsl_matrix_complex * Mtmp;
  gsl_matrix_complex * Ntmp;

  for (k = 0; k < Ngen; k++) {
    ta[k] = gsl_matrix_complex_calloc (Ngen, Ngen);
  } 

  for (a = 0; a < Ngen; a++) {
    for (b = 0; b < Ngen; b++) {
      for (c = 0; c < Ngen; c++) {

        Mtmp = gsl_matrix_complex_calloc (Nrep, Nrep);
        Ntmp = gsl_matrix_complex_calloc (Nrep, Nrep);

        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, alpha, tf[a], tf[b], beta, Mtmp);
        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, alpha, Mtmp, tf[c], beta, Ntmp);

        gsl_complex ftmp = gsl_complex_rect (0.0, 0.0);;

	for (i = 0; i < Nrep; i++){
          ftmp = gsl_complex_add (ftmp, gsl_matrix_complex_get (Ntmp, i, i));
        }

        Mtmp = gsl_matrix_complex_calloc (Nrep, Nrep);
        Ntmp = gsl_matrix_complex_calloc (Nrep, Nrep);

        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, alpha, tf[b], tf[a], beta, Mtmp);
        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, alpha, Mtmp, tf[c], beta, Ntmp);

        for (i = 0; i < Nrep; i++){
          ftmp = gsl_complex_sub (ftmp, gsl_matrix_complex_get (Ntmp, i, i));
        }

	gsl_matrix_complex_set (ta[a], b, c, gsl_complex_mul_real (ftmp, -1.0/M_SQRT2));

      }
    }
  }

//Print out the generators of SU(N) in the adjoint rep

  for (k = 0; k < Ngen; k++)
    for (i = 0; i < Ngen; i++) 
      for (j = 0; j < Ngen; j++)
        printf ("ta(%d,%d,%d) = %10g%10g\n", k, i, j, GSL_REAL(gsl_matrix_complex_get (ta[k], i, j)), GSL_IMAG(gsl_matrix_complex_get (ta[k], i,j)));

//gsl_matrix_free
  for (k = 0; k < tmp; k++){
    gsl_matrix_complex_free (tf[k]);
    gsl_matrix_complex_free (ta[k]);
  }

  gsl_matrix_complex_free (Mtmp);
  gsl_matrix_complex_free (Ntmp);

  return 0;

}
