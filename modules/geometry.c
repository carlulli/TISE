/******************************************************************************
Module reads input values from main
sets parameters by assining input values to parameter variables
has get_params functions that return the parameters
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>

/*
defined inside c file: variables can only be accessed from this file
and returned from get_functions OR accessed from other files with extern keyword
static: can not be accessed with extern keyword (safest way)
*/
static int N=0;

/*
function set_params
takes input values
transform input char to int
assign input int to variables (global?):
NUM, larger than 1 and odd;
A;
eps (tolerance), around 10e-15;
 */
void set_params(int argc, char *argv[]){
  int val;
  val = atoi(argv[1]);
  if (val > 1 && val % 2) {
    N = val;
  } else {
    printf("[ geometry.c| set_params ] NUM has to be > 1 and an odd number!\n");

    exit(-1);
  }
}

/* This is another way to return N (to try out the adress/reference stuff)
void get_N(int *pN) {
    static int sN=0;

    if(N==0)
    {
      printf("[ geometry.c| get_N ] Error! Not not yet set!\n");
      exit(-1);
    }

    if(sN==0)
    {
      sN=N;
    }
    else
    {
      if((N!=sN))
      {
         printf("[ geometry.c| get_N ] Error! (N) has changed: (%d) -> (%d)\n",sN,N);
         exit(-1);
      }
    }
    (*pN)=N;
}*/

int get_N() {
    static int sN=0;

    if(N==0)
    {
      printf("[ geometry.c| get_N ] Error! Not not yet set!\n");
      exit(-1);
    }

    if(sN==0)
    {
      sN=N;
    }
    else
    {
      if((N!=sN))
      {
         printf("[ geometry.c| get_N ] Error! (N) has changed: (%d) -> (%d)\n",sN,N);
         exit(-1);
      }
    }
    return N;
}

/* prints the currently set parameters */
void print_N() {
  printf("The parameters set from the input are: " "NUM =\t%d \n",N );
}
