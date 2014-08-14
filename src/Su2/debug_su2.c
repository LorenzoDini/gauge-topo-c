#ifndef DEBUG_SU2_C
#define DEBUG_SU2_C

#include"math.h"

#include"../Const/const.h"
#include"../Macro/macro.h"
#include"../Rng/random.h"
#include"su2.h"
#include"su2_upd.h"

int main(void)
  {
  int seme=0;
  double energy;
  Const param;
  
  Su2 M, N, L, T, mI;
   
  /* initialize random seed */
  initrand(seme);

  /* fix a value for d_beta */
  param.d_beta=2.3;
    
  printf("\n*******************************\n");
  printf("PROGRAM FOR THE DEBUG OF SU(2)\n");
  printf("*******************************\n\n");

  printf("Precision: double");
 

 
  printf("\n\n");
  printf("VERIFY THAT THE RANDOM MATRIX IS IN SU(2)\n\n");
  printf("  random matrix ....");
  rand_matrix_Su2(&M);
  one_Su2(&mI);
  times_equal_real_Su2(&mI, -1.0);  /* mI=-1.0 */
  times_dag2_Su2(&T, &M, &M);
  plus_equal_Su2(&T, &mI);
  if(sqrtdet_Su2(&T) <=MIN_VALUE) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }



  printf("\n\n");
  printf("VERIFY THAT UPDATE SU(2)->SU(2)\n\n");
  rand_matrix_Su2(&M);
  rand_matrix_Su2(&N);
  rand_matrix_Su2(&L);
  plus_equal_Su2(&N, &L); /* N+=L,  M in SU(2), N no   (M=link, N=staple) */

  /* heatbath */
  single_heatbath_Su2(&M, &N, &param);
  printf("  Heatbath ...");
  times_dag2_Su2(&T, &M, &M);
  plus_equal_Su2(&T, &mI);
  if(sqrtdet_Su2(&T) <=MIN_VALUE) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }

  /* overrelaxation */
  rand_matrix_Su2(&M);
  rand_matrix_Su2(&N);
  rand_matrix_Su2(&L);
  plus_equal_Su2(&N, &L); /* N+=L,  M in SU(2), N no   (M=link, N=staple) */
  single_overrelaxation_Su2(&M, &N, &param);
  printf("  Overrelaxation ...");
  times_dag2_Su2(&T, &M, &M);
  plus_equal_Su2(&T, &mI);
  if(sqrtdet_Su2(&T) <=MIN_VALUE) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }



  printf("\n\n");
  printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
  rand_matrix_Su2(&M);
  rand_matrix_Su2(&N);
  rand_matrix_Su2(&L);
  plus_equal_Su2(&N, &L); /* N+=L,  M in SU(2), N no   (M=link, N=staple) */

  times_Su2(&T, &M, &N);  /* T=M*N */
  energy=retr_Su2(&T);    /* initial energy */
  single_overrelaxation_Su2(&M, &N, &param);
  times_Su2(&T, &M, &N);  /* T=M*N */
  energy-=retr_Su2(&T);    /* -=final energy */
  if(fabs(energy)<MIN_VALUE) 
    {
    printf("  OK\n");
    }
  else 
    {
    printf("  ERROR!!!!!!!!!!!   DeltaE=%g\n", energy);
    }



  printf("\n\n");
  printf("VERIFY THAT COOLING DECREASES THE ENERGY ...");
  rand_matrix_Su2(&M);
  rand_matrix_Su2(&N);
  rand_matrix_Su2(&L);
  plus_equal_Su2(&N, &L); /* N+=L,  M in SU(2), N no   (M=link, N=staple) */

  times_Su2(&T, &M, &N);  /* T=M*N */
  energy=retr_Su2(&T);    /* initial energy */

  cool_Su2(&M, &N);

  times_Su2(&T, &M, &N);  /* T=M*N */
  energy-=retr_Su2(&T);    /* -=final energy */
  if(energy<MIN_VALUE) 

    {
    printf("  OK\n");
    }
  else 
    {
    printf("  ERROR!!!!!!!!!!!   DeltaE=%g\n", energy);
    }

  printf("\n\n");

  return 0;
  }

#endif
