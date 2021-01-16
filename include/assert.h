#ifndef ASSERT_H
#define ASSERT_H

#include <stdarg.h>

/*******************************************************************************
If test fails, this function prints an error message and exit the execution with
return value -1.

This is what is printed in the standard output:
   ERROR [filename] function
   message
where the message is built out of the
   const char* format,...
parameters as in the printf function.
*******************************************************************************/
void assert(int test,const char* filename,const char* functionname,const char* format,...);


#endif /* ASSERT_H */
