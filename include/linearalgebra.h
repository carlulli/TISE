#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <complex.h>

#define dcompl double complex

dcompl scalar_product(dcompl* a, dcompl* b,int N);
/* takes two complex vectors and the lenght of the vectors
and retruns one complex number */



double norm(dcompl* a,int N);
/*analogous but only computes the norm*/



#endif // !LinearAlgebra_h
