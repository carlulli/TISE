/******************************************************************************
Module reads input values from main
sets parameters by assining input values to parameter variables
has get_params functions that return the parameters
******************************************************************************/
#ifndef GEOMETRY_H
#define GEOMETRY_H

void set_params(int argc, char *argv[]);
/*
takes input values
transform input char to int
assign input int to variables (global?):
NUM, larger than 1 and odd;
 */

/*void get_N(int *pN);*/
int get_N();
 /* returns N */

void print_N();
/* prints the currently set parameters */

#endif /* GEOMETRY_H */
