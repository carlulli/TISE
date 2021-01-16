#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "geometry.h"





double complex dot_product(double complex a[], double complex b[], int length)
{
    double complex result = 0.0 + 0 * I;
    for (int i = 0; i< length; i++)
    {
        result += conj(a[i]) * b[i];
    }
    return result;

}








void conj_grad ( double complex b[], double complex x[], void(*pfunc)(double complex *, double complex*))    // § need to change the parameters of pfunc depending on which pfunc
/*   for solving A* x = b  , in our case we do not have A, but a method to get A*x

*/
{

    int N = get_N();
    int i;    // for the loops
    double  rsold;
    double  rsnew;
    double complex alpha;
    double complex r[N];
    double complex p[N];
    double complex Ap[N];
    double complex Ax[N] ;
    double accuracy = 1.e-3;

    double complex dot_product(double complex a[], double complex b[], int length);

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

    for ( int it = 1; it <= N; it++ )
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
        rsnew = dot_product(r,r,N);
/************************************************************************
   print rsnew for testing
***************************************************************************/

        printf("\n rsnew %f,  \n",rsnew);

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
