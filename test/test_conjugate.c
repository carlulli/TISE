#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "conjugategradient.h"
#include "geometry.h"
#include "linearalgebra.h"
#include "assert.h"


/******************************************************************************
test program:  compiles with gcc -Wall -o ${TDIR}/test_conjugate.o -I ${IDIR}
${MDIR}/geometry.c ${MDIR}/assert.c ${MDIR}/linearalgebra.c
 ${MDIR}/conjugategradient.c ${TDIR}/test_conjugate.c -lm   --> see compile.sh


 generates random Matrix B
 creates random positive semidefinite Matrix A = B_transpose * B



******************************************************************************/

static int N;
double accuracy;
int kfac;


// /*Function for acting with matrix on vector*/
// void multiplyvec(double complex **mat, double complex *vec, double complex *res)
// {
//     int i, k;
//     for (i = 0; i < N; i++)
//     {
//         res[i] = 0;
//         for (k = 0; k < N; k++)
//             res[i] += mat[i][k] * vec[k];
//     }
//     return;
// }
// /* Function for creating a semi definite positive random Matrix*/
// void multAtimesv( double complex* in, double complex* out)
// {
//     static double complex **A=NULL;
//     if(A==NULL)
//     {
//         /* Allocate the memory for the matrix A */
//         A=malloc(N*sizeof(double complex*));
//         for(int i=0;i<N;i++)
//             A[i]=malloc(N*sizeof(double complex));
//
//         /* Allocate the memory for the matrix C */
//         double complex **C=NULL;
//         C=malloc(N*sizeof(double complex*));
//         for(int i=0;i<N;i++)
//         {
//             C[i]=malloc(N*sizeof(double complex));
//         }
//
//
//         /* Set C equal to a random matrix */
//         for(int i = 0; i < N; i++)
//         {
//             for(int k = 0; k < N; k++)
//                 C[i][k] = 1.0/N * rand()/RAND_MAX ;
//         }
//
//         /* Set A equal to Cdag.C */
//         for (int i = 0; i < N ; i++)
//         {
//             for (int j = 0; j < N; j++)
//             {
//                 double complex sum = 0;
//                 for( int k = 0; k<N; k++)
//                     sum += conj(C[k][i])*C[k][j];
//
//                 A[i][j] = sum;
//             }
//         }
//
//         /* Free the memory for C. */
//         for(int i=0;i<N;i++) free(C[i]);
//         free(C);
//     }
//     /* Applying the created Matrix on the input vector*/
//     multiplyvec(A,in,out);
// }
//
//
// /**************************************************************************
//  function for random vector
// ***************************************************************************/
// void randvec(double complex vec[N], int N)
// {
//     srand(time(NULL));
//     for(int i = 0; i < N; i++)
//     {
//         vec[i] = 1.0 * rand()/RAND_MAX + 1.0 * rand()/RAND_MAX *I;
//     }
// }



#define _FILE_NAME_ "test/test_conjugate.c"
// takes parameters N dimension of the vector and accuracy (=residue) for the testing
int main(int argc, char* argv[]) {

  assert(argc==4,_FILE_NAME_,"main"," For this routine input the following: [N] [residue] [kfac] \n \
  kfac being the integer value of the factor of the identity matrix, which is going to be added to our random matrix.");
  //use geometry.h to take the parameter
  set_params(argc,argv);
  N = get_N();
  accuracy = atof(argv[2]);
  kfac = atof(argv[3]);
  double complex x[N];

  double complex vec[N];
  double complex b[N] ;
  double complex c[N] ;
  double complex result[N];
  double deviation;
  double differ;

  srand(time(NULL));

  /*
Printing introduction to the program, explaining what it does
  */

printf("\n Welcome to the conjugate gradient test module. It creates  a positive semi definite Matrix A, and a random vector vec.\
Then it applies the Matrix A on the random Vector vec. The resulting vector is fed as input parameter to the conjugate gradient function.\
It provides us with the solution x for A * x = b. Then the matrix A is applied to this vecor x. " );


/************************************************************************************
creates and prints random vector
************************************************************************************/
    randvec(vec,N);
    printf("\n  random vector: vec\n");
    for(int i = 0; i < N; i++)
    {
        printf("%f, <--real and imaginary--> %f \n",creal(vec[i]),cimag(vec[i]));
    }
/************************************************************************************
pointer for matrix * vector multiplication
************************************************************************************/
    void multAtimesv( double complex* in, double complex* out);

    void (*pmultAtimesv)( double complex* in, double complex* out);
    pmultAtimesv = &multAtimesv;
/************************************************************************************
multiply our positive semi def random Matrix A with a random vector vec  A* vec = b
************************************************************************************/
    multAtimesv(vec,b);
    printf("\n  vector b from A*vec = b\n");
    for(int i = 0; i < N; i++)
    {
        printf("%f, <--real and imaginary--> %f \n",creal(b[i]),cimag(b[i]));
    }
/************************************************************************************
conjugate gradient that should solve A*X = b
************************************************************************************/
    conj_grad(b,x,pmultAtimesv);
/************************************************************************************
print out x
************************************************************************************/
    printf("\n  vextor x from conj grad A *x = b xxxx\n");
    for(int i = 0; i < N; i++)
    {
        printf("%f, <--real and imaginary--> %f \n",creal(x[i]),cimag(x[i]));
    }
/************************************************************************************
multiply A*x to see if its similar to b
************************************************************************************/
    multAtimesv(x,c);
    printf("\n  c für A*x = c cccc\n");
    for(int i = 0; i < N; i++)
    {
        printf("%f, <--real and imaginary--> %f \n",creal(c[i]),cimag(c[i]));
    }

    for(int i = 0; i < N; i++)
    {
        result[i] = c[i]- b[i];
    }
    deviation = norm(result,N);
    differ = accuracy - deviation;
    printf("\n The deviation: abs(A*x - b) is %.20f .\n The accuracy parameter, which should be achieved, was %.20f \n\n.\
     The difference (accuracy minus deviation) is %.20f \n \
     If this value is positive our desired accuracy was achieved. \n ", deviation, accuracy, differ);




    return 0;
}
