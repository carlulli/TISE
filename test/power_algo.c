#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>


#include "powermethod.h"
#include "geometry.h"


//
double complex* M;
static int N = 5;

//Creates the needed hermitian def pos matrix of NxN dimension
void initialize_M() {
  M = malloc(sizeof(double complex)*N*N);
  double complex A[N][N];
  double complex sum;
  srand(time(NULL));
  /* create random A */
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
        sum += conj(A[i][k])*A[k][j];
      }
      M[i*N + j] = sum;
      printf("M[%d][%d] = %f %f\n",i,j,creal(sum) ,cimag(sum));
    }
  }
}
/* applies M to w such that w (in[])
 --> Mw  = z (out[])       */
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
// takes parameters N dimension of the vector and tolerance tol for the testing
int main(int argc, char* argv[]) {
  //use geometry.h to take the parameter
  set_params(argc,argv);
  N = get_N();
  double complex out_ev[N];
  double tol = atof(argv[2]);
  void (*ptr)(double complex*,double complex*) = Mvect;

  initialize_M();
  double mu = power_method(tol,ptr,out_ev);
  free(M);
  return 0;
}
