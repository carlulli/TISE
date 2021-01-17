#include <stdio.h>
#include <stdlib.h>


#include "geometry.h"
#include "hamiltonian.h"
#include "powermethod.h"
#include "observables.h"
#include "linearalgebra.h"
#include "conjugategradient.h"

/*******************************************************************************
Wrapper function that calculates the inverse of the hamiltonian apllied to v.
It is necessary because conjugategradient() can not be given to the power method
directly, as it needs more parameters than the function declared in the power method.

H_inv:
  - takes 2 parameters:
      - input and output vectors (actually their pointers)
  - calls the conjugate gradient with these vectors and
    the hamiltonian that has to be set before
*******************************************************************************/
//void H_defpos(double complex *in, double complex *out);
//void (*pH_defpos)(double complex *in, double complex *out);
 //pH_defpos = &H_defpos;
void H_inv(double complex *in, double complex *out){
  conj_grad(in, out, &H_defpos);
}
//void (*pH_inv)(double complex *in, double complex *out);
//pH_inv = &H_inv;
/*******************************************************************************
calling main with parameters argv:
  argv[1] = N
  argv[2] = m
  argv[3] = number to decide on the potential
  argv[4] = tolerance (for powermethod)
  argv[5] = residue (for conjugate gradient)
therefore argc should be 6
*******************************************************************************/
int main(int argc, char* argv[]) {
  if (argc == 6) {
  printf(
    "STARTING your TISE calculations with desired input parameters: \n"
    "NUM = %s\n" "mass = %s\n" "potential number = %s\n"
    "tolrance for powermethod = %s\n" "residue for conjugate gradient = %s\n",
    argv[1], argv[2], argv[3], argv[4], argv[5]
  );
  FILE *fp3;

  fp3 = fopen("output/input_parameters.txt", "w");
  fprintf(fp3,
    "Input parameters you chose your simulation.\n"
    "NUM\t%s\n" "mass\t%s\n" "potential number\t%s\n" "tolrance\t%s\n" "residue\t%s\n",
    argv[1], argv[2], argv[3], argv[4], argv[5]
  );

  fclose(fp3);
  }
  else {
    printf("\n[ main.c | input parameters ] ERROR. Necessary number of input parameters 6!\n"
    "You have %d!\n" "Remember: The executable Name is the first parameter.\n", argc);
    exit(-1);
  }


  /*****************************************************************************
  set hamiltonian with the parameters mass m, potential
  *****************************************************************************/
  set_params(argc, argv);
  int N = get_N();
  print_N();
  double m = atof(argv[2]);
  set_res(argv);
  set_kinetic_params(m);
  set_potential(atoi(argv[3]));
  print_hamiltonian_info();
  set_res(argv);
  get_res();
  /****************************************************************************
  calculating the largest eigenvalue and eigenvector of the inverse of H
  power method to calculate largest eigenvalue is called with:
    - tolerance: double tol
    - output vector (eigenstate): double complex *out_ev
    - function that calculates Mv: void M (double complex *in, double complex*out)
        this is gonna be the H_inv declared above, where conjugategradient is called
   conjugate gradient takes input, output and hamiltonian as parameters:
        - input vector: b = w (the starting vector from the power method
                              (new in every iteration))
        - output vector: H^-1b = z (z from power method theor. z=Mv)
    power method output:
     - int mu = lowest eigenvalue = groundstate energy
     - double v[N] = eigenstate corresponding to groundstate energy
  ****************************************************************************/
  double complex eigenstate[N];
  for(int i = 0; i < N; i++)
  {
      eigenstate[i] = 1.0 * rand()/RAND_MAX + 1.0 * rand()/RAND_MAX *I;
  }
  double eigenvalue;
  /* if tol is handled inside power method module, it needs a set_tol() and
  does not need to be passed the tolerance
  set_tol(); */
  double tol = atof(argv[4]);
  printf("Powermethod tolerance =\t%f\n", tol);
/* power_method(tol, &H_inv(), eigenstate) where H_inv() also gets N  */
  eigenvalue = power_method(tol, &H_inv, eigenstate);

  /* getting the actual ev of H
     needs a get_k() from the hamiltonian
     I am not to sure about when we need to substract or add
  */
  int k = get_minV();
  if (k < 0) {
    eigenvalue += k;
  }
  else {
    eigenvalue -= k;
  }
  printf("ANNOUNCMENT: The eigenvalue is: %f\n", eigenvalue);

  /****************************************************************************
  Calculate the observables
  avg_x, delta_x, avg_p, delta_p
  The get_ functions should probably be defined in a seperate module
  Alternatively they can be defined in here.
  ****************************************************************************/
  double avg_x, delta_x, avg_p, delta_p;
  avg_x = get_avgx(eigenstate);
  delta_x = get_deltax(eigenstate);
  avg_p = get_avgp(eigenstate);
  delta_p = get_deltap(eigenstate);

  printf("[main.c | ] NOTICE: If the position width = %f is close to N/2 = %d, you are getting finite volume effects!\n"
  "\t If so, you should be simulating with larger N!\n", delta_x, N/2);

  /****************************************************************************
  create file with eigenstate values
  later a seperate file with eigenvalue and the other observables is created
  ATTENTION: maybe the last \n might make problems when reading in the file again
  ****************************************************************************/
  FILE *fp1;

  fp1 = fopen("output/output_ES.txt", "w");
  fprintf(fp1, "n\tREAL(eigenstate[n])\tIMAG(eigenstate[n])");
  for (int i=0; i<N; i++) {
    fprintf(fp1, "%f\t%f\n", creal(eigenstate[i]), cimag(eigenstate[i]));
  }
  fclose(fp1);


  /****************************************************************************
  create file with eigenvalue and the other observables
  ****************************************************************************/
  FILE *fp2;

  fp2 = fopen("output/output_observ.txt", "w");
  fprintf(fp2,
    "eigenvalue\t%f\n" "<x^>\t%f\n" "delta_n\t%f\n" "<p^>\t%f\n" "delta_n\t%f",
    eigenvalue, avg_x, delta_x, avg_p, delta_p
  );

  fclose(fp2);


  return 0;
}
