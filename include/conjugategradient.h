#include <complex.h>


void conj_grad ( double complex b[], double complex x[], void(*pfunc)(double complex *, double complex*));
void set_residue(double res);
double get_residue();
