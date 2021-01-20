#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include "geometry.h"


/* Calculates the complex scalar product of two complex vectors */
double complex scalar_product(double complex* a, double complex* b, int N) {
	double complex sum = { 0 };
	for (int i = 0; i < N; i++) {
		sum += conj(a[i]) * b[i];
	}
	return sum;
}

/* Calculates the norm */
double norm(double complex* a, int N) {

	double sum = 0;
	for (int i = 0; i < N ; i++) {
  	sum += creal(a[i])*creal(a[i]) + cimag(a[i])*cimag(a[i]);
}
return sqrt(sum);
}


/*Function for acting with matrix on vector*/
void multiplyvec(double complex **mat, double complex *vec, double complex *res)
{
    int i, k;
		int N = get_N();
    for (i = 0; i < N; i++)
    {
        res[i] = 0;
        for (k = 0; k < N; k++)
            res[i] += mat[i][k] * vec[k];
    }
    return;
}
/* Function for creating a semi definite positive random Matrix*/
void multAtimesv( double complex* in, double complex* out)
{
	 	int N = get_N();
		static double complex **A=NULL;
    if(A==NULL)
    {
        /* Allocate the memory for the matrix A */
        A=malloc(N*sizeof(double complex*));
        for(int i=0;i<N;i++)
            A[i]=malloc(N*sizeof(double complex));

        /* Allocate the memory for the matrix C */
        double complex **C=NULL;
        C=malloc(N*sizeof(double complex*));
        for(int i=0;i<N;i++)
        {
            C[i]=malloc(N*sizeof(double complex));
        }


        /* Set C equal to a random matrix */
        for(int i = 0; i < N; i++)
        {
            for(int k = 0; k < N; k++)
                C[i][k] = 1.0/N * rand()/RAND_MAX ;
        }

        /* Set A equal to Cdag.C */
        for (int i = 0; i < N ; i++)
        {
            for (int j = 0; j < N; j++)
            {
                double complex sum = 0;
                for( int k = 0; k<N; k++)
                    sum += conj(C[k][i])*C[k][j];

                A[i][j] = sum;
            }
        }

        /* Free the memory for C. */
        for(int i=0;i<N;i++) free(C[i]);
        free(C);
    }
    /* Applying the created Matrix on the input vector*/
    multiplyvec(A,in,out);
}

/**************************************************************************
 function for random vector
***************************************************************************/
void randvec(double complex* vec, int M)
{
		//srand(time(NULL));
    for(int i = 0; i < M; i++)
    {
        vec[i] = 1.0 * rand()/RAND_MAX + 1.0 * rand()/RAND_MAX *I;
    }
}
