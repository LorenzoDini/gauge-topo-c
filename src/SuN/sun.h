#ifndef SUN_H
#define SUN_H

#include<complex.h>
#include<stdio.h>

#include"../Macro/macro.h"

typedef struct SuN {
   double complex comp[NCOLOR][NCOLOR];
} SuN;

void one_SuN(SuN *__restrict__ A);         /* A=1 */
void zero_SuN(SuN *__restrict__ A);        /* A=0 */

void equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B);      /* A=B       */
void equal_dag_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B);  /* A=B^{dag} */

void plus_equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B);     /* A+=B       */
void plus_equal_dag_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B); /* A+=B^{dag} */

void minus_equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B);     /* A-=B       */
void minus_equal_dag_SuN(SuN *__restrict A, SuN const *__restrict__ const B); /* A-=B^{dag} */

void lin_comb_SuN(SuN *__restrict__ A, 
                  double b, 
                  SuN const *__restrict__ const B, 
                  double c, 
                  SuN const *__restrict__ const C);       /* A=b*B+c*C */
void lin_comb_dag1_SuN(SuN *__restrict__ A, 
                       double b, 
                       SuN const *__restrict__ const B, 
                       double c, 
                       SuN const *__restrict__ const C);  /* A=b*B^{dag}+c*C */
void lin_comb_dag2_SuN(SuN *__restrict__ A, 
                       double b, 
                       SuN const *__restrict__ const B,  
                       double c, 
                       SuN const *__restrict__ const C);  /* A=b*B+c*C^{dag} */
void lin_comb_dag12_SuN(SuN *__restrict__ A, 
                       double b, 
                       SuN const *__restrict__ const B, 
                       double c, 
                       SuN const *__restrict__ const C); /* A=b*B^{dag}+c*C^{dag} */

void times_equal_real_SuN(SuN *__restrict__ A, double r); /* A*=r */

void times_equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B);     /* A*=B       */
void times_equal_dag_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B); /* A*=B^{dag} */

void times_SuN(SuN *__restrict__ A,
               SuN const *__restrict__ const B, 
               SuN const *__restrict__ const C);       /* A=B*C             */
void times_dag1_SuN(SuN *__restrict__ A, 
                    SuN const *__restrict__ const B, 
                    SuN const *__restrict__ const C);  /* A=B^{dag}*C       */
void times_dag2_SuN(SuN *__restrict__ A, 
                    SuN const *__restrict__ const B,  
                    SuN const *__restrict__ const C);  /* A=B*C^{dag}       */
void times_dag12_SuN(SuN *__restrict__ A, 
                     SuN const *__restrict__ const B, 
                     SuN const *__restrict__ const C); /* A=B^{dag}*C^{dag} */

void rand_matrix_SuN(SuN *__restrict__ A);

double norm_SuN(SuN const *__restrict__ const A);
double retr_SuN(SuN const *__restrict__ const A);
double imtr_SuN(SuN const *__restrict__ const A);

void LU_SuN(SuN const *__restrict__ const A, SuN *__restrict__ ris, int *sign);
complex double det_SuN(SuN const *__restrict__ const A);
int scheck_SuN(SuN const *__restrict__ const A); /* gives 0 if the matrix is in SU(N) and 1 otherwise */
void unitarize_aux_SuN(SuN *__restrict__ A);
void unitarize_SuN(SuN *__restrict__ A);

void print_on_screen_SuN(SuN const *__restrict__ const A);
void print_on_file_SuN(FILE *fp, SuN const *__restrict__ const A);
void read_from_file_SuN(FILE *fp, SuN *__restrict__ A);


#endif
