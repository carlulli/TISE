#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


/* Calculates the complex scalar product of two complex vectors */
double complex scalar_product(double complex* a, double complex* b, int N) {
	double complex sum = { 0 };
	for (int i = 0; i < N; i++) {
		sum += conj(a[i]) * b[i];
	}
	return sum;
}

/* Calculates the norm */
double norm(double complex* a, int N) {

	double sum = 0;
	for (int i = 0; i < N ; i++) {
  	sum += creal(a[i])*creal(a[i]) + cimag(a[i])*cimag(a[i]);
}
return sqrt(sum);
}
