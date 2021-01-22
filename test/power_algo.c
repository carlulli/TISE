#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>


#include "powermethod.h"
#include "geometry.h"
#include "linearalgebra.h"
#include "assert.h"


double complex* M;
static int N = 5;

/* Creates the needed hermitian def pos matrix of NxN dimension */
void initialize_M() {
  M = malloc(sizeof(double complex)*N*N);
  double complex A[N][N];
  double complex sum;
/*  srand(time(NULL)); is called in beginning of main */
  /* creates random A */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A[i][j] = 1.0*rand()/RAND_MAX + 1.0*rand()/RAND_MAX*I;
    }
  }
  /* calculate M semi definite positive as M = (A^t)*A  */
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      sum = 0;
      for(int k = 0; k < N; k++){
        sum += conj(A[k][i])*A[k][j];
      }
      M[i*N + j] = sum;
    }
  }
}
/* applies M to w such that w (in[]) --> Mw  = z (out[])    */
void Mvect(double complex* in, double complex* out) {
  double complex sum;
  for(int i = 0; i < N; i++) {
    sum = 0;
    for(int j = 0; j < N; j++) {
      sum += M[i*N + j]*in[j];
    }
    out[i] = sum;
  }

}


#define _FILE_NAME_ "test/power_algo.c"
// takes parameters N dimension of the vector and tolerance tol for the testing
int main(int argc, char* argv[]) {

  assert(argc==3,_FILE_NAME_,"main"," For this routine input the following: [N] [TOLERANCE]");
  //use geometry.h to take the parameter
  set_params(argc,argv);
  N = get_N();
  double tol = atof(argv[2]);
  double complex out_ev[N];
  double mu;
  void (*ptr)(double complex*,double complex*) = Mvect;

  initialize_M();
  mu = power_method(tol,ptr,out_ev);
  printf("final mu = %f \n",mu);

  /* Here we study the distance between the two the analytical result and the approximated eigenvalue/eigenvector */

  double complex Mw_out[N],muw_out[N],delta[N];

  for(int i = 0; i < N; i++) {
    muw_out[i]  = mu*out_ev[i];
  }

  Mvect(out_ev,Mw_out);

  for(int i = 0; i < N; i++) {

    delta[i] = Mw_out[i] - muw_out[i];
  }
  printf("Note that the parameter for the tolerance shouldn't be smaller than 1e-14\n");
  printf("Computing || Mw - muw|| = %.12e < %.12e\nIf the output is smaller than the tolerance input, the test is correct!\n",norm(delta,N), tol);


  free(M);
  return 0;
}
