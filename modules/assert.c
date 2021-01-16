#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


void assert(int test,const char* filename,const char* functionname,const char* format,...)
{
   if(test) return;

   va_list arglist;
   printf("\nERROR [%s] %s\n",filename,functionname);
   va_start(arglist,format);
   vprintf(format,arglist);
   va_end(arglist);
   printf("\n\n");
   exit(-1);
}
