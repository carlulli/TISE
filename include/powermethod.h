#ifndef POWERMETHOD_H
#define POWERMETHOD_H

#include <complex.h>

double power_method(double tol,
  void (*M)(double complex* in,double complex* out),
  double complex *out_ev);


#endif
