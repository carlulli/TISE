#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include <complex.h>

// test
void conj_grad ( double complex b[], double complex x[], void(*pfunc)(double complex *, double complex*));

void set_res(char *argv[]);

double get_res();


#endif
