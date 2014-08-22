#ifndef SON_H
#define SON_H

#include<stdio.h>

#include"../Macro/macro.h"

typedef struct SoN {
   #ifdef __GNUC__
     double comp[NCOLOR][NCOLOR] __attribute__ ((aligned (DOUBLE_ALIGN)));
   #else
     double comp[NCOLOR][NCOLOR];
   #endif
} SoN;

void one_SoN(SoN *__restrict__ A);         /* A=1 */
void zero_SoN(SoN *__restrict__ A);        /* A=0 */

void equal_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B);      /* A=B       */
void equal_dag_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B);  /* A=B^{dag} */

void plus_equal_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B);     /* A+=B       */
void plus_equal_dag_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B); /* A+=B^{dag} */

void minus_equal_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B);     /* A-=B       */
void minus_equal_dag_SoN(SoN *__restrict A, SoN const *__restrict__ const B); /* A-=B^{dag} */

void lin_comb_SoN(SoN *__restrict__ A, 
                  double b, 
                  SoN const *__restrict__ const B, 
                  double c, 
                  SoN const *__restrict__ const C);       /* A=b*B+c*C */
void lin_comb_dag1_SoN(SoN *__restrict__ A, 
                       double b, 
                       SoN const *__restrict__ const B, 
                       double c, 
                       SoN const *__restrict__ const C);  /* A=b*B^{dag}+c*C */
void lin_comb_dag2_SoN(SoN *__restrict__ A, 
                       double b, 
                       SoN const *__restrict__ const B,  
                       double c, 
                       SoN const *__restrict__ const C);  /* A=b*B+c*C^{dag} */
void lin_comb_dag12_SoN(SoN *__restrict__ A, 
                       double b, 
                       SoN const *__restrict__ const B, 
                       double c, 
                       SoN const *__restrict__ const C); /* A=b*B^{dag}+c*C^{dag} */

void times_equal_real_SoN(SoN *__restrict__ A, double r); /* A*=r */

void times_equal_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B);     /* A*=B       */
void times_equal_dag_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B); /* A*=B^{dag} */

void times_SoN(SoN *__restrict__ A,
               SoN const *__restrict__ const B, 
               SoN const *__restrict__ const C);       /* A=B*C             */
void times_dag1_SoN(SoN *__restrict__ A, 
                    SoN const *__restrict__ const B, 
                    SoN const *__restrict__ const C);  /* A=B^{dag}*C       */
void times_dag2_SoN(SoN *__restrict__ A, 
                    SoN const *__restrict__ const B,  
                    SoN const *__restrict__ const C);  /* A=B*C^{dag}       */
void times_dag12_SoN(SoN *__restrict__ A, 
                     SoN const *__restrict__ const B, 
                     SoN const *__restrict__ const C); /* A=B^{dag}*C^{dag} */

void rand_matrix_SoN(SoN *__restrict__ A);

double norm_SoN(SoN const *__restrict__ const A);
double retr_SoN(SoN const *__restrict__ const A);
double imtr_SoN(SoN const *__restrict__ const A);

void LU_SoN(SoN const *__restrict__ const A, SoN *__restrict__ ris, int *sign);
double det_SoN(SoN const *__restrict__ const A);
int scheck_SoN(SoN const *__restrict__ const A); /* gives 0 if the matrix is in SU(N) and 1 otherwise */
void unitarize_aux_SoN(SoN *__restrict__ A);
void unitarize_SoN(SoN *__restrict__ A);

void print_on_screen_SoN(SoN const *__restrict__ const A);
void print_on_file_SoN(FILE *fp, SoN const *__restrict__ const A);
void read_from_file_SoN(FILE *fp, SoN *__restrict__ A);


#endif
