#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <complex.h>


double complex scalar_product(double complex* a, double complex* b,int N);
/* takes two complex vectors and the lenght of the vectors
and retruns one complex number */

double norm(double complex* a,int N);
/*analogous but only computes the norm*/

/*Function for acting with matrix on vector*/
void multiplyvec(double complex **mat, double complex *vec, double complex *res);

/* Function for creating a semi definite positive random Matrix*/
void multAtimesv( double complex* in, double complex* out);

/* function for random vector */
void randvec(double complex* vec, int M);

#endif // !LinearAlgebra_h
