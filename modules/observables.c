#include <complex.h>
#include <math.h>
#include <stdio.h>

/*******************************************************************************
geometry.h declares the get_N() function
*******************************************************************************/
#include "geometry.h"

/*******************************************************************************
Calculates the discretized normalization of a complex vector
(usually some state psi) by:

N = scalar_product(adj(psi), psi) = sum_{n=0}^{N-1} cabs(psi(n))^2
*******************************************************************************/
double Normal(double complex *psi) {
    int N = get_N();
    double sum = 0;
    for (int n=0; n<N; n++) {
      sum += creal(psi[n])*creal(psi[n]) + cimag(psi[n])*cimag(psi[n]);
    }
    return sum;
}

/*******************************************************************************
Calculates the average position in lattice units of a complex vector by:

<x> = 1/N sum_{n=0}^{N-1} n*cabs(psi(n))^2
*******************************************************************************/
double avgx(double complex *psi) {
  int N = get_N();
  double sum = 0;
  for (int n=0; n<N; n++) {
    sum += n * ( creal(psi[n])*creal(psi[n]) + cimag(psi[n])*cimag(psi[n]) );
  }
/*  printf("[observables.c | avgx()] sum =\t%f\t Norm =\t%f\n" , sum, Normal(psi)); */
  return sum/Normal(psi);
}

/*******************************************************************************
Calculates the (dimensionaless) width of a complex vector by:

deltax = sqrt( 1/N * (sum_{n=0}^{N-1} n^2*cabs(psi(n))^2) - <x>^2 )

Question: does the 1/N only belong to the sum?
*******************************************************************************/
double deltax(double complex *psi) {
  int N = get_N();
  double sum = 0;
  for (int n=0; n<N; n++) {
/*    sum += n*n * ( creal(psi[n])*creal(psi[n]) + cimag(psi[n])*cimag(psi[n]) ); */
    sum+= ( creal(psi[n])*creal(psi[n]) +  cimag(psi[n])*cimag(psi[n]) )
          * (n - avgx(psi)) * (n - avgx(psi)) ;
    /*printf("[observables.c | deltax] sum[%d] =\t%f\n", n, sum);*/
  }

  /*double test = (1/Normal(psi)) * sum - avgx(psi) * avgx(psi);
  printf("[observables.c | deltax] test =\t%f\n", test);
  double test2 = (1/Normal(psi)) * sum;
  printf("[observables.c | deltax] test2 =\t%f\n", test2);
  double test4 = (1/Normal(psi));
  printf("[observables.c | deltax] test4 =\t%f\n", test4);
  printf("[observables.c | deltax] sum =\t%f\n", sum);
  double testNormal = Normal(psi);
  printf("[observables.c | deltax] Normal =\t%f\n", testNormal);
  double test3 = avgx(psi) * avgx(psi);
  printf("[observables.c | deltax] test3 =\t%f\n", test3);*/

/*  return sqrt( (1/Normal(psi)) * sum - avgx(psi) * avgx(psi) ); */
/*  return sqrt( 1/Normal(psi) * sum ); */
  return sqrt( sum / Normal(psi) );
}

/*******************************************************************************
Calculates the average momentum of a complex vector by:

<p> = 1/N * sum_{n=0}^{N-1} (-I) * scalar_product(conj(psi(n)), derivative(psi(n)))

with the symmetric derivative:

derivative(psi(n)) = 1/2 * (psi(n+1) - psi(n-1))

ATTENTION: here we use Dirichlet Boundary Conditions,
            setting psi((N-1)+1)=psi(N)=0 and psi(0-1)=psi(-1)=0

Possible systematic ERROR: if psi is not a eigenfunction of the momentum operator
          the values at sum are actually complex number and not real
          BUT this might not raise an error, even if it is wrong

Readbility: is it better to write the if, else if, else like done her
            or do you prefere it as done in the function below
*******************************************************************************/
double avgp(double complex *psi) {
  int N = get_N();
  double sum=0.0;
  for (int n=0; n<N; n++) {
    if (n==0) { sum +=  (-I) * conj(psi[n]) * (psi[n+1])/2; }
    else if (n==N-1) { sum +=  (-I) * conj(psi[n]) * ((-1)*psi[n-1])/2; }
    else {  sum +=  (-I) * conj(psi[n]) * ((psi[n+1]-psi[n-1]))/2; }
  /*  printf("[observables.c | avgp()] sum[%d] =\t%f\n", n, sum); */
  }
/*  printf("[observables.c | avgp()] sum =\t%f\n", sum); */
  return (sum) / Normal(psi);
}

/*******************************************************************************
Calculates the average squared momentum of a complex vector by:

<p^2> = 1/N * sum_{n=0}^{N-1} cabs(derivative(psi(n)))^2

Idea how to calculate cabs(dervative(psi(n))):
1. derivative(psi(n)) = 1/2 * (psi(n+1) - psi(n-1)) =: z
2. cabs(z) = conj(z) * z

ATTENTION: here we use Dirichlet Boundary Conditions,
            setting psi((N-1)+1)=psi(N)=0 and psi(0-1)=psi(-1)=0
*******************************************************************************/
double avgsqr_p(double complex *psi) {
  int N = get_N();
  double sum = 0.0;
  for (int n=0; n<N; n++) {
    if (n==0) {
      sum += conj( (psi[n+1])/2 ) * (psi[n+1])/2 ;
    }
    else if (n==N-1) {
      sum += conj( (-1)*psi[n-1]/2 ) * (-1) * psi[n-1]/2 ;
    }
    else {
      sum += conj( ( psi[n+1]-psi[n-1] )/2 ) * ( psi[n+1]-psi[n-1] )/2 ;
    }
/*    printf("[observables.c | avgsqr_p()] sum=\t%f\n", sum); */
  }
  // double test = sum/Normal(psi);
/*  printf("[observables.c | avgsqr_p()] sum/Normalization=\t%f\n", test); */
  return (sum) / Normal(psi);
}

/*******************************************************************************
Calculates the momentum width of a complex vector by:

deltap = sqrt( <p^2> - <p>^2 )
*******************************************************************************/
double deltap(double complex *psi) {
/*  printf("[observables.c | deltap)()] avgp=\t%f\n"
  "[observables.c | deltap()] avgsqr_p=\t%f\n",
  avgp(psi), avgsqr_p(psi) ); */
  return sqrt(avgsqr_p(psi) - avgp(psi) * avgp(psi) );
}

/*******************************************************************************
Get functions:
That return the actual values that are calculated.
*******************************************************************************/
double get_avgx(double complex *psi) {
  return avgx(psi);
}

double get_deltax(double complex *psi) {
  return deltax(psi);
}

double get_avgp(double complex *psi) {
  return avgp(psi);
}

double get_deltap(double complex *psi) {
  return deltap(psi);
}

double get_avgsqr_p(double complex *psi) {
  return avgsqr_p(psi);
}
