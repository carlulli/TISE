#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "linearalgebra.h"
#include "geometry.h"
#include "assert.h"
#include "hamiltonian.h"



#define dcompl double complex

#define _FILE_NAME_ "test/hamiltonian/analytic.c"

static int N = 0;



void set_sinwave(double complex *psi,int k)
{
   for(int n=0;n<N;n++)
      psi[n]=sin((M_PI*((n+1)*k))/(N+1));
}

int main (int argc, char *argv[]) {
  assert(argc== 2,_FILE_NAME_,"main","Usage: analytic N\nwhere N is the number of lattice points");
  /* HERE THE MASS IS HARDCODED */
  double mass = 2.3513;
  set_params(argc,argv);
  N = get_N();
  //set_N(N);//put the parameter into geometry.c
  set_kinetic_params(mass);//define parameters in hamiltonian
  set_zero_potential();//sets a potential of 0
  print_hamiltonian_info();
  printf("\nThis programs checks the eigenvalues and eigenvectors of the Hamiltonian\n"
       "with V=0. In this case one can calculate analytically that the eigenvectors\n"
       "are labeled by an integer index k. The k-th eigenvector is given by\n"
       "   psi(n) = sin( pi*(n+1)*k/(N+1) )     for n=0,1,2,...,N-1\n"
       "and the corresponding eigenvalue is given by\n"
       "   E = 2/m * sin( pi*k/(2*(N+1)) )^2\n\n"
       "For each k=0,1,2,...,N this program prints\n"
       "   maxdev = max_n | H.psi(n) - E psi(n) |\n"
       "The test passes is all these numbers are < 1e-15\n\n");

       dcompl in[N],out[N];
       double ev,dev,maxdev;

       for(int k=0;k<=N;k++)
       {
          set_sinwave(in,k);
          H(in,out);
          ev=2.*sin((M_PI*k)/(2*(N+1)))*sin((M_PI*k)/(2*(N+1)))/mass;
          printf("Average state energy %f\n",creal(average_state_energy(out)));
          maxdev=0.0;
          for(int n=0;n<N;n++)
          {
             dev=cabs(out[n]-ev*in[n]);
             if(dev>maxdev) maxdev=dev;
          }
          printf("k= %d\tmaxdev= %.2e\n",k,maxdev);
       }
       printf("\n");



return 0;

}
