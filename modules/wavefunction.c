#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "linearalgebra.h"

void set_random_wavefunction(double complex* psi,int N) {

  // srand(time(NULL)); is called in beginning of the main.c

  for(int  i = 0; i < N; i++) {
    psi[i] = 1.0*rand()/RAND_MAX + 1.0*rand()/RAND_MAX*I;
  }
  double psi_norm = norm(psi,N);
  for(int  i = 0; i < N; i++) {
    psi[i]  = psi[i]/psi_norm;

  }
}
