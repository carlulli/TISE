#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "linearalgebra.h"
#include "geometry.h"
#include "assert.h"
#include "hamiltonian.h"
#include "wavefunction.h"

#define _FILE_NAME_ "test/linearity.c"


int main (int argc, char* argv[]) {

  assert(argc==5,_FILE_NAME_,"main","[N] [MASS] [SET_POTENTIAL] [TOLERANCE]");

  /* Here the parameters of the system are set */
  double mass = atof(argv[2]);
  double tol = atof(argv[4]);
  set_params(argc,argv);     //geometry library
  set_kinetic_params(mass);  //hamiltonian
  set_potential(atoi(argv[3])); //hamiltonian


  int N = get_N();   //array size
  double complex psi[N],theta[N]; // needed wavefunctions
  double complex out_psi[N],out_theta[N],out_left[N],out_right[N];
  double complex in_left[N],delta[N];
  double complex alpha,beta;

  /*gives random wavefunction*/
  set_random_wavefunction(psi,N);
  set_random_wavefunction(theta,N);

  /* test printf to check good functions */
  for(int i = 0; i < N; i++) {
    printf("psi[%d]  = %f %fi\n",i, creal(psi[i]),cimag(psi[i]));
  }
  printf("norm psi = %f \nnorm theta = %f\n", norm(psi,N),norm(theta,N));

  /* generate two random coefficients */
  srand(time(NULL));
  alpha = 1.0*rand()/RAND_MAX + 1.0*rand()/RAND_MAX*I;
  beta = 1.0*rand()/RAND_MAX + 1.0*rand()/RAND_MAX*I;


  for(int i = 0; i < N; i++) {
    in_left[i] = alpha*psi[i] + beta*theta[i];
  }
  H(in_left,out_left);//computes the left hand side

  H(psi,out_psi);
  H(theta,out_theta);
  for(int i = 0; i < N; i++) {
    out_right[i] = alpha*out_psi[i] + beta*out_theta[i];
    delta[i] = out_left[i] - out_right[i];
  }

  printf("The program calculates the distance between the two sides of the Linearity equation :\n ||H(a_psi + b_theta) - a(H_psi) + b(H_theta) || distance between the two sides : \n difference = %.12e\nDesidered tolerance %.12e \n", norm(delta,N), tol);
  return 0;
}
