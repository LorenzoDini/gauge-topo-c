#ifndef DEBUG_RNG_C
#define DEBUG_RNG_C

#include<stdio.h>

#include"../Macro/macro.h"
#include"random.h"

int main(void)
   {
   int i;
   int seme;

   printf("\n");

   seme=1;
   initrand(seme);
   printf("seed=%d\n", seme);
   for(i=0; i<5; i++)
      {
      printf("  random[%d]=%.16g\n", i, casuale());
      }

   printf("\n");

   seme=2;
   initrand(seme);
   printf("seed=%d\n", seme);
   for(i=0; i<5; i++)
      {
      printf("  random[%d]=%.16g\n", i, casuale());
      }

   printf("\n");

   seme=1;
   initrand(seme);
   printf("seed=%d\n", seme);
   for(i=0; i<5; i++)
      {
      printf("  random[%d]=%.16g\n", i, casuale());
      }

   printf("\n");

   seme=0;
   initrand(seme);
   printf("seed=%d\n", seme);
   for(i=0; i<5; i++)
      {
      printf("  random[%d]=%g\n", i, casuale());
      }

   printf("\n");



   return 0;
   }

#endif
