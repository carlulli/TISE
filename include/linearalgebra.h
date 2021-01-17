#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <complex.h>



double complex scalar_product(double complex* a, double complex* b,int N);
/* takes two complex vectors and the lenght of the vectors
and retruns one complex number */



double norm(double complex* a,int N);
/*analogous but only computes the norm*/



#endif // !LinearAlgebra_h
