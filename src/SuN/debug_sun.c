#ifndef DEBUG_SUN_C
#define DEBUG_SUN_C

#include"math.h"

#include"../Const/const.h"
#include"../Macro/macro.h"
#include"../Rng/random.h"
#include"sun.h"
#include"sun_upd.h"


int main(void)
  {
  int seme=0;
  double energy;
  Const param;
  
  SuN M, N, L, T;
   
  /* initialize random seed */
  initrand(seme);

  /* fix a value for d_beta */
  param.d_beta=2.3;
    
  printf("\n*******************************\n");
  printf("PROGRAM FOR THE DEBUG OF SU(N)\n");
  printf("*******************************\n\n");

  printf("Precision: double");
 

 
  printf("\n\n");
  printf("VERIFY THAT THE RANDOM MATRIX IS IN SU(N)\n\n");
  printf("  random matrix ....");
  rand_matrix_SuN(&M);
  if(scheck_SuN(&M) == 0) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }

  printf("\n\n");
  printf("VERIFY THAT UPDATE SU(N)->SU(N)\n\n");
  rand_matrix_SuN(&M);
  rand_matrix_SuN(&N);
  rand_matrix_SuN(&L);
  plus_equal_SuN(&N, &L); /* N+=L,  M in SU(N), N no   (M=link, N=staple) */

  /* heatbath */
  single_heatbath_SuN(&M, &N, &param);
  printf("  Heatbath ...");
  if(scheck_SuN(&M) == 0) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }

  /* overrelaxation */
  rand_matrix_SuN(&M);
  rand_matrix_SuN(&N);
  rand_matrix_SuN(&L);
  plus_equal_SuN(&N, &L); /* N+=L,  M in SU(N), N no   (M=link, N=staple) */
  single_overrelaxation_SuN(&M, &N, &param);
  printf("  Overrelaxation ...");
  if(scheck_SuN(&M) == 0) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }



  printf("\n\n");
  printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
  rand_matrix_SuN(&M);
  rand_matrix_SuN(&N);
  rand_matrix_SuN(&L);
  plus_equal_SuN(&N, &L); /* N+=L,  M in SU(N), N no   (M=link, N=staple) */

  times_SuN(&T, &M, &N);  /* T=M*N */
  energy=retr_SuN(&T);    /* initial energy */
  single_overrelaxation_SuN(&M, &N, &param);
  times_SuN(&T, &M, &N);  /* T=M*N */
  energy-=retr_SuN(&T);    /* -=final energy */
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
  rand_matrix_SuN(&M);
  rand_matrix_SuN(&N);
  rand_matrix_SuN(&L);
  plus_equal_SuN(&N, &L); /* N+=L,  M in SU(N), N no   (M=link, N=staple) */

  times_SuN(&T, &M, &N);  /* T=M*N */
  energy=retr_SuN(&T);    /* initial energy */

  cool_SuN(&M, &N);

  times_SuN(&T, &M, &N);  /* T=M*N */
  energy-=retr_SuN(&T);    /* -=final energy */
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
