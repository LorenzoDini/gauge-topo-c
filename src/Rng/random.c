#ifndef RANDOM_C
#define RANDOM_C

#include<math.h>
#include<time.h>

#include"./dSFMT/dSFMT.h"
#include"random.h"
#include"../Macro/macro.h"

/* random number in (0,1) */
double casuale(void)
 {
 return dsfmt_gv_genrand_open_open();
 }


/* initialize random generator */
void initrand(unsigned long s)
  {
  if(s==0)
    {
    dsfmt_gv_init_gen_rand(time(NULL));
    }
  else
    {
    dsfmt_gv_init_gen_rand(s);
    }
  }

/* normal gaussian random number generator (polar method, knuth vol 2, p. 117) */
double gauss1()
   {
   double v1, v2, s, ris;

   do
     {
     v1=1.0-2.0*casuale();
     v2=1.0-2.0*casuale();
     s=v1*v1+v2*v2;
     }
   while(s >= 1);

   ris=v1*sqrt(-2*log(s)/s);
    return ris;
   }

#endif
