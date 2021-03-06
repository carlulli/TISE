#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>


#include "linearalgebra.h"
#include "geometry.h"


static int N = 0;

/*
Initializing eigenvalue and eigenvector outside funtion
to make them available for other methods ?!
- Where does the pointer have to be initialized?
*/

double power_method(double tol,void (*M)(double complex*,
  double complex*),
  double complex *out_ev) {
  srand(time(NULL));
  N = get_N();
  double complex w[N];
  double complex z[N];
  double complex delta[N];
  double mu;
  double delta_norm = 2*tol;
  int iter = 0;
  /*
  generate random starting vector w[N]
  Maybe we should initialize a random sinusiod or even sin + icos?
  */
  for(int n = 0; n < N; n++) {
    w[n] = 1.0 * rand()/RAND_MAX +
    1.0 * rand()/RAND_MAX*I;
  }

  while(delta_norm > tol){
    M(w,z);
    mu = norm(z,N);
    //printf("mu = %.20f\t  \n",mu);

    for (int n = 0; n < N; n++) {
      delta[n] = z[n] - mu*w[n];
    }
    // testing  printf("\n %f",creal(w[n]));

    for (int n = 0; n < N; n++) {
    w[n] = z[n]/mu;
    }
    delta_norm = norm(delta,N);
    iter++;

    if(iter > 100) {
      printf("The algorithm could not converge, please try again with a higher tolerance \n");
      exit(-1);
    }
    }

    for (int n = 0; n < N; n++) {
      out_ev[n] = w[n];
    }
    return mu;
}
