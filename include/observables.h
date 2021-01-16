#include <complex.h>
#include <math.h>

/*******************************************************************************
Get functions:
That return the actual values that are calculated.
The functions that calculate the values are inside observables.c
Do they need to be referenced here, if the are not called in the main?
*******************************************************************************/
double get_avgx(double complex *psi);
double get_deltax(double complex *psi);
double get_avgp(double complex *psi);
double get_deltap(double complex *psi);
double get_avgsqr_p(double complex *psi);

double Normal(double complex *psi);
