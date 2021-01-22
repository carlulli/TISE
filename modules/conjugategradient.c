#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "geometry.h"
#include "hamiltonian.h"
#include "linearalgebra.h"



double accuracy;

/*******************************************************************************
set_res initialzes the residue value by assigning the 5th input parameter to res
get_res returns residue
residue is called accuracy inside conjugategradient.c
*******************************************************************************/
void set_res(char *argv[]) {
  double val = 0;
  val = atof(argv[5]);
  if (val > 0) { accuracy = val; }
  else {
    printf("[conjugategradient | set_residue()] ERROR: residue has to be > 0!");
    exit(-1);
  }
}

double get_res() {
  if (accuracy == 0) {
    printf("[conjugategradient | set_res()] ERROR: Residue not initialzed yet!");
    exit(-1);
  }
  return accuracy;
}


void conj_grad ( double complex b[], double complex x[], void(*pfunc)(double complex *, double complex*))
/*   for solving A* x = b  , in our case we do not have A, but a method to get A*x

*/
{

    int N = get_N();
    //double accuracy = get_res();
    int i;
    /************************************************************************
    testing
    **************************************************************************/
  /*
    for ( i = 0; i < N; i++ )
    {
      printf("\n conj grad kriegt\n b \t x \n %f \t%f \t %f \n ",creal(b[i]),cimag(b[i]),creal(x[i]) );
    }*/
/*******************************************************************************/
    double  rsold;
    double  rsnew;
    double complex alpha;
    double complex r[N];
    double complex p[N];
    double complex Ap[N];
    double complex Ax[N] ;
/************************************************************************
    1. calcultate A*x
***************************************************************************/
for(int n = 0; n < N; n++) {
  x[n] = 1.0 * rand()/RAND_MAX +
  1.0 * rand()/RAND_MAX*I;
}

    pfunc(x, Ax);

    /*******************************************************
    testing
    **********************************************************/

  /*  for(int i = 0; i < N; i++)
    {
        printf("\n%f, und imaginär %f ",creal(Ax[i]),cimag(Ax[i]));
    }*/
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

    rsold = scalar_product(r,r,N);

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

        alpha = rsold / ( scalar_product(p,Ap,N)  ) ;
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
        rsnew = scalar_product(r,r,N);
/************************************************************************
   print rsnew for testing
***************************************************************************/

      //  printf("[conjugategradient.c | ] INFO: rsnew\t%f  \n",rsnew);

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

        //printf("[conjugategradient.c | ]  loop number %i\n",it);

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
    if (sqrt(rsnew) > accuracy)     // uses accuracy
        {
            printf("Attention: The desired attention was not achieved, \n instead we only got to %.15f \n", sqrt(rsnew));
            exit(-1);
        }
    return;
}
