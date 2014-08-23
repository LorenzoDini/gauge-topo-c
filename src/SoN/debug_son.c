#ifndef DEBUG_SUN_C
#define DEBUG_SUN_C

#include"math.h"

#include"../Const/const.h"
#include"../Macro/macro.h"
#include"../Rng/random.h"
#include"son.h"
#include"son_upd.h"


int main(void)
  {
  int seme=0;
  double energy;
  Const param;
  
  SoN M, N, L, T;
   
  /* initialize random seed */
  initrand(seme);

  /* fix a value for d_beta */
  param.d_beta=2.3;
    
  printf("\n*******************************\n");
  printf("PROGRAM FOR THE DEBUG OF SO(N)\n");
  printf("*******************************\n\n");

  printf("Precision: double\n");
  printf("N_c=%d\n", NCOLOR);

 
  printf("\n\n");
  printf("VERIFY THAT THE RANDOM MATRIX IS IN SO(N)\n\n");
  printf("  random matrix ....");
  rand_matrix_SoN(&M);
  if(scheck_SoN(&M) == 0) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }


  printf("\n\n");
  printf("VERIFY THAT UPDATE SO(N)->SO(N)\n\n");
  rand_matrix_SoN(&M);
  rand_matrix_SoN(&N);
  rand_matrix_SoN(&L);
  plus_equal_SoN(&N, &L); /* N+=L,  M in SO(N), N no   (M=link, N=staple) */

  /* heatbath */
  single_heatbath_SoN(&M, &N, &param);
  printf("  Heatbath ...");
  if(scheck_SoN(&M) == 0) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }

  /* overrelaxation */
  rand_matrix_SoN(&M);
  rand_matrix_SoN(&N);
  rand_matrix_SoN(&L);
  plus_equal_SoN(&N, &L); /* N+=L,  M in SO(N), N no   (M=link, N=staple) */
  single_overrelaxation_SoN(&M, &N, &param);
  printf("  Overrelaxation ...");
  if(scheck_SoN(&M) == 0) 
    {
    printf("    OK\n");
    }
  else 
    {
    printf("    ERROR!!!!!!!!!!!\n");
    }



  printf("\n\n");
  printf("VERIFY THAT OVERRELAXATION DOES NOT CHANGE THE ENERGY ...");
  rand_matrix_SoN(&M);
  rand_matrix_SoN(&N);
  rand_matrix_SoN(&L);
  plus_equal_SoN(&N, &L); /* N+=L,  M in SO(N), N no   (M=link, N=staple) */

  times_SoN(&T, &M, &N);  /* T=M*N */
  energy=retr_SoN(&T);    /* initial energy */
  single_overrelaxation_SoN(&M, &N, &param);
  times_SoN(&T, &M, &N);  /* T=M*N */
  energy-=retr_SoN(&T);    /* -=final energy */
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
  rand_matrix_SoN(&M);
  rand_matrix_SoN(&N);
  rand_matrix_SoN(&L);
  plus_equal_SoN(&N, &L); /* N+=L,  M in SO(N), N no   (M=link, N=staple) */

  times_SoN(&T, &M, &N);  /* T=M*N */
  energy=retr_SoN(&T);    /* initial energy */

  cool_SoN(&M, &N);

  times_SoN(&T, &M, &N);  /* T=M*N */
  energy-=retr_SoN(&T);    /* -=final energy */
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
