#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>



/******************************************************************************
test program:  compiles with gcc -Wall test_conjugate.c -o test_conjugate -lm

 generates random Matrix B
 creates random positive semidefinite Matrix A = B_transpose * B



******************************************************************************/

/**************************************************************************
 !!!!  here N has to be selected in all three lines
***************************************************************************/

const int N = 30;
double complex x[30] = {0};


int getN()
{
    return N;
}

void multiplyvec(double complex **mat, double complex *vec, double complex *res)
{
    int N=getN();

    int i, k;
    for (i = 0; i < N; i++)
    {
        res[i] = 0;
        for (k = 0; k < N; k++)
            res[i] += mat[i][k] * vec[k];
    }
    return;
}


void multAtimesv( double complex* in, double complex* out)
{
    /* A static variable inside a function is inizialized only the first time
       the function is called, then its value is remembered every time the
       function is called again */

    static double complex **A=NULL;
    int N=getN();

    /* Notice that this block is executed only the first time the function is
       called */
    if(A==NULL)
    {
        /* Allocate the memory for the matrix A */
        A=malloc(N*sizeof(double complex*));
        for(int i=0;i<N;i++)
            A[i]=malloc(N*sizeof(double complex));
        /* Notice: if there is only one instruction in the for loop, you do not
           need the brackets { ... } */

        /* Allocate the memory for the matrix C */
        double complex **C=NULL;
        C=malloc(N*sizeof(double complex*));
        for(int i=0;i<N;i++)
        {
            C[i]=malloc(N*sizeof(double complex));
        }

        /* This is a logic error: the seen should be set only once in the whole
           program, typically at the beginning of the main function, and *NOT*
           everytime you need a random number */
        srand(time(NULL));

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

        /* Free the memory for C. Notice that you should not free A, since you
          will need it after */
        for(int i=0;i<N;i++) free(C[i]);
        free(C);
    }

    multiplyvec(A,in,out);
}


/**************************************************************************
 function for random vector
***************************************************************************/
void randvec(double complex vec[N], int N)
{
    srand(time(NULL));
    for(int i = 0; i < N; i++)
    {
        vec[i] = 1.0 * rand()/RAND_MAX + 1.0 * rand()/RAND_MAX *I;
    }
}

/**************************************************************************
  function for A_dagger* A  (given A)
***************************************************************************/


/************************************************************************************
function for dot product used in Conj gradient
************************************************************************************/

double complex dot_product(double complex a[], double complex b[], int length)
{
    double complex result = 0.0 + 0 * I;
    for (int i = 0; i< length; i++)
    {
        result += conj(a[i]) * b[i];
    }
    return result;

}
/**************************************************************************
 The conjugent gradient  for solving A* x = b
 , in our case we do not have A, but a method to get A*x
***************************************************************************/

void conj_grad ( double complex b[], double complex x[], void(*pfunc)(double complex *, double complex*))
{


   //int get_N();
    //N = get_N();

    int i;
    double   rsold;
    double   rsnew;
    double complex alpha;
    double complex r[N];
    double complex p[N];
    double complex Ap[N];
    double complex Ax[N] ;
    double accuracy = 1.e-3;

   // double complex dot_product(double complex a[], double complex b[], int length);

/************************************************************************
    1. calcultate A*x
***************************************************************************/

    pfunc(x, Ax);

/********************************************************
     2. initialize r = b-Ax0
***************************************/

    for ( i = 0; i < N; i++ )
    {
        r[i] = b[i]-Ax[i];
    }

/************************************************************************
    3.  set p = r
***************************************************************************/

    for ( i = 0; i < N; i++ )
    {
        p[i] = r[i];
    }
/************************************************************************
    4. rsold = r' * r
***************************************************************************/

    rsold = dot_product(r,r,N);

/************************************************************************
   This is where the Loop starts
***************************************************************************/

    for ( int it = 1; it < 1000; it++ )
    {
/************************************************************************
         5. Ap = A * p   done with pfunc
***************************************************************************/
        pfunc(p, Ap);

/************************************************************************
    alpha = rsold / (p'*Ap)
***************************************************************************/                                                              //alpha = dot_product(p,r,N) / ( dot_product(p,Ap,N)  ) ;// uses N  /option 1

        alpha = rsold / ( dot_product(p,Ap,N)  ) ;
/************************************************************************
    set x = x + alpha * p
***************************************************************************/
        for ( i = 0; i < N; i++ )
        {
            x[i] = x[i] + alpha * p[i];
        }
/************************************************************************
   r = r - alpha * Ap
***************************************************************************/

        for ( i = 0; i < N; i++ )
        {
            r[i] = r[i] - alpha * Ap[i];
        }

/************************************************************************
   rsnew = r' * r
***************************************************************************/                                                                                                                                                //rsnew = dot_product(r,Ap,N);    // uses N   option 1   change
                                                                                                                                        //rsold = dot_product(p,Ap,N);    // option 1             chan
        rsnew = creal(dot_product(r,r,N));
/************************************************************************
   print rsnew for testing
***************************************************************************/

        printf("\n rsnold %f rsnew %f accuracy %f,  \n",rsold,rsnew,accuracy);

/************************************************************************
   break condition
***************************************************************************/

        if (sqrt(rsnew) < accuracy)     // uses accuracy
        {
            break;
        }

/************************************************************************
  p = r + (rsnew / rsold) * p
***************************************************************************/

        for ( i = 0; i < N; i++ )
        {

            p[i] = r[i] + (rsnew / rsold) * p[i];

        }

/************************************************************************
   prints loop number and can print x
***************************************************************************/

        printf("loop number %i",it);
        /*printf("\n x bei loop %i\n",it);
        for(int i = 0; i < N; i++)
        {
           // printf("%f, und imaginär %f \n",creal(x[i]),cimag(x[i]));
        } */

/************************************************************************
   rsold = rsnew
***************************************************************************/

        rsold = rsnew;

/************************************************************************
   for testing
***************************************************************************/

        /* printf("p%.5f \n",creal( p[0]));
         printf("p%.5f \n",creal( p[1]));
         printf("p%.5f \n",creal( p[2]));
         printf("p%.5f \n",cimag( p[0]));
         printf("p%.5f \n",cimag( p[1]));
         printf("p%.5f \n",cimag( p[2]));  */

    }



    return;

}











int main()
{

    double complex vec[N];
    double complex b[N] ;
    double complex c[N] ;


/************************************************************************************
creates and prints random vector
************************************************************************************/

    randvec(vec,N);
    printf("\n  random vector: vec\n");
    for(int i = 0; i < N; i++)
    {
        printf("%f, und imaginär %f \n",creal(vec[i]),cimag(vec[i]));
    }

//
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
        printf("%f, und imaginär %f \n",creal(b[i]),cimag(b[i]));
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
        printf("%f, und imaginär %f \n",creal(x[i]),cimag(x[i]));
    }
/************************************************************************************
multiply A*x to see if its similar to b
************************************************************************************/
    multAtimesv(x,c);
    printf("\n  c für A*x = c cccc\n");
    for(int i = 0; i < N; i++)
    {
        printf("%f, und imaginär %f \n",creal(c[i]),cimag(c[i]));
    }

    return 0;
}
