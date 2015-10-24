/* 
  DEBUG FOR G2 GROUP
*/

#include<math.h>
#include<stdio.h>

#include"../GParam/gparam.h"
#include"g2.h"
#include"g2_check.h"
#include"g2_un.h"
#include"g2_upd.h"
#include"../Macro/macro.h"
#include"../Rng/random.h"

int main(void)
   {
   int count;
   int seme=0;
   double n, cc, energy;
   GParam param;

   G2 M, N, L;

   /* inizialize random number generator */
   initrand(seme);

   /* fix a value for d_beta */
   param.d_beta=2.3;

   printf("\n****************************\n");
   printf("PROGRAM FOR DEBUG OF THE G2 GROUP\n");
   printf("****************************\n\n");
  
   printf("Precision: double");

   printf("\n\n");
   printf("VERIFY THAT FUNDAMENTAL ROTATIONS ARE IN G2\n\n");

   printf("  matrix D1  ....");
   D1(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D2  ....");
   D2(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D3  ....");
   D3(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D4  ....");
   D4(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D5  ....");
   D5(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D6  ....");
   D6(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D7  ....");
   D7(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D8  ....");
   D8(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D9  ....");
   D9(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D10 ....");
   D10(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D11 ....");
   D11(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D12 ....");
   D12(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D13 ....");
   D13(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   printf("  matrix D14 ....");
   D14(&M, casuale());
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);



   printf("\n\n");
   printf("VERIFY THAT THE RANDOM MATRIX IS IN G2\n\n");
   printf("  random matrix ....");
   rand_matrix_G2(&M);
   if(check_G2(&M, &n, &cc)) printf("    OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);



   printf("\n\n");
   printf("VERIFY THAT UPDATE G2->G2\n\n");

   rand_matrix_G2(&M);
   rand_matrix_G2(&N);
   rand_matrix_G2(&L);
   plus_equal_G2(&N, &L);    /*  M in G2, N no   (M=link, N=staple) */ 
   single_heatbath_G2(&M, &N, &param);
   printf("  Heatbath ...");
   if(check_G2(&M, &n, &cc)) printf("  OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);

   rand_matrix_G2(&M);
   rand_matrix_G2(&N);
   rand_matrix_G2(&L);
   plus_equal_G2(&N, &L);    /*  M in G2, N no   (M=link, N=staple) */
   single_overrelaxation_G2(&M, &N, &param);
   printf("  Overrelaxation ...");
   if(check_G2(&M, &n, &cc)) printf("  OK\n");
   else printf("    ERROR!!!!!!!!!!!   n=%.16g  c=%.16g\n", n, cc);



   printf("\n\n");
   printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
   rand_matrix_G2(&M);
   rand_matrix_G2(&N);
   rand_matrix_G2(&L);
   plus_equal_G2(&N, &L);    /*  M in G2, N no   (M=link, N=staple) */

   times_G2(&L, &M, &N);
   energy=retr_G2(&L);  /* initial energy */

   single_overrelaxation_G2(&M, &N, &param);

   times_G2(&L, &M, &N);
   energy-=retr_G2(&L);  /* initial - final */

   if(fabs(energy)<MIN_VALUE) printf("  OK\n");
   else printf("  ERROR!!!!!!!!!!!   DeltaE=%g\n", energy);



   printf("\n\n");
   printf("VERIFT PROJECTION ON G2\n");
   printf("    limit precision: %g\n", MIN_VALUE);
   count=0;    
   rand_matrix_G2(&M);
   while(check_G2(&M, &n, &cc))
        {
        rand_matrix_G2(&N);
        times_equal_G2(&M, &N);
        count++;
        } 
   printf("    product of rand_matrix_G2() is out of G2 after %d iterations\n", count);
   printf("    n=unitary condition\t\tcc=cubic condition\n");  
   printf("    n=%g  cc=%g\n", n, cc);
   printf("    projection ...");
   unitarize_G2(&M);
   if(check_G2(&M, &n, &cc)) printf("  OK\n");
   else printf("  ERROR!!!!!!!!!!!\n");
   printf("    n=%g  cc=%g\n", n, cc);

   printf("\n\n");
   printf("VERIFY THAT COOLING DECREASES THE ENERGY ...");
   rand_matrix_G2(&M);
   rand_matrix_G2(&N);
   rand_matrix_G2(&L);
   plus_equal_G2(&N, &L);    /*  M in G2, N no   (M=link, N=staple) */

   times_G2(&L, &M, &N);
   energy=retr_G2(&L);  /* initial energy */

   cool_G2(&M,&N);

   times_G2(&L, &M, &N);
   energy-=retr_G2(&L);  /* initial -final */

   if(energy<MIN_VALUE)  printf("  OK\n");
   else printf("  ERROR!!!!!!!!!!!   DeltaE=%g\n", energy);
   printf("\n\n");

   return 0;
   }
