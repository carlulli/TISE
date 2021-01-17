#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include "linearalgebra.h"
#include "geometry.h"
#include "assert.h"
#include "hamiltonian.h"
#include "wavefunction.h"

#define _FILE_NAME_ "test/new_hamiltonian/hermitian.c"

static int N = 0;


int main (int argc, char* argv[]) {
  assert(argc==3,_FILE_NAME_,"main","Usage:
  analytic N\nwhere N is the number of lattice
  points\n  And the second digit the assigned potential");
  double mass = atof(argv[2]);
  set_params(argc,argv);
  set_kinetic_params(mass);
  set_potential(atoi(argv[3]));
  N = get_N();

  //We now wish to study hermitianity
  printf("This program studies the Hermicity of the hamiltonian through it's action on a randomly generated wave function \n");
  double complex psi[N],left[N],right[N];
  set_random_wavefunction(psi);
  for(int i = 0; i < N; i++) {
    printf("psi[%d] = %f %fi\n",i,creal(psi[i]),cimag(psi[i]));
  }
  printf("norm of psi = %f", norm(psi,N));
  H(psi,left);
  H(psi,right);
  double complex delta[1];
  delta[0] = scalar_product(left,psi,N) - scalar_product(psi,right,N);
  printf(" ||(Hpsi,psi) - (psi,Hpsi)|| = %.12e\n",norm(delta,1));




  return 0;
}
