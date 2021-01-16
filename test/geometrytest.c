#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "geometry.h"

/* this was just to try the random thing out  */
static double drand()
{
   return (1.0*rand())/RAND_MAX;
}

int main(int argc, char *argv[]) {
/*
  printf(
    "argc: %d\n
    argv[0]: %s\n
    argv[1]: %d\n"
    , argc, argv[0], atoi(argv[1]));
*/

  if (argc < 2) {
    printf("You forgot your NUM parameter.\n"
    "Remeber: executable [NUM]!\n"
  "That is why the following zsh error is raised.\n\n");
  exit(-1);
    }

  /* When running the executable after compiling:
    argv[0] is the executable name
    argv[1] is the NUM parameter (an number after the exe name)
      NUM has to be bigger than 0 and odd (so the center is a lattice point)
  */
  set_params(argc, argv);
  int N = get_N();
  print_N();

  printf("Check if N from Module = N from main: ModuleN see before; thisN= %d\n", N);

  double rand;
  rand = drand();
  printf("Random number: %f\n", rand);

  return 0;
}
