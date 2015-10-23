#ifndef SUN_H
#define SUN_H

#include"../Macro/macro.h"

#include<complex.h>
#include<stdio.h>

typedef struct SuN {
   double complex comp[Ncolor][Ncolor];
} SuN;

void one_SuN(SuN *A);         /* A=1 */
void zero_SuN(SuN *A);        /* A=0 */

void equal_SuN(SuN *A, SuN const * const B);      /* A=B       */
void equal_dag_SuN(SuN *A, SuN const * const B);  /* A=B^{dag} */

void plus_equal_SuN(SuN *A, SuN const * const B);     /* A+=B       */
void plus_equal_dag_SuN(SuN *A, SuN const * const B); /* A+=B^{dag} */

void minus_equal_SuN(SuN *A, SuN const * const B);     /* A-=B       */
void minus_equal_dag_SuN(SuN *A, SuN const * const B); /* A-=B^{dag} */

void lin_comb_SuN(SuN *A, double b, SuN const * const B, double c, SuN const * const C);       /* A=b*B+c*C */
void lin_comb_dag1_SuN(SuN *A, double b, SuN const * const B, double c, SuN const * const C);  /* A=b*B^{dag}+c*C */
void lin_comb_dag2_SuN(SuN *A, double b, SuN const * const B, double c, SuN const * const C);  /* A=b*B+c*C^{dag} */
void lin_comb_dag12_SuN(SuN *A, double b, SuN const * const B, double c, SuN const * const C); /* A=b*B^{dag}+c*C^{dag} */

void times_equal_real_SuN(SuN *A, double r); /* A*=r */

void times_equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B);     /* A*=B       */
void times_equal_dag_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B); /* A*=B^{dag} */

void times_SuN(SuN * A, SuN const * const B, SuN const * const C);       /* A=B*C             */
void times_dag1_SuN(SuN *A, SuN const * const B, SuN const * const C);  /* A=B^{dag}*C       */
void times_dag2_SuN(SuN *A, SuN const * const B, SuN const * const C);  /* A=B*C^{dag}       */
void times_dag12_SuN(SuN *A, SuN const * const B, SuN const * const C); /* A=B^{dag}*C^{dag} */

void rand_matrix_SuN(SuN *A);

double norm_SuN(SuN const * const A);
double retr_SuN(SuN const * const A);
double imtr_SuN(SuN const * const A);

void LU_SuN(SuN const * const A, SuN *ris, int *sign);
complex double det_SuN(SuN const * const A);
int scheck_SuN(SuN const * const A); /* gives 0 if the matrix is in SU(N) and 1 otherwise */
void unitarize_aux_SuN(SuN *A);
void unitarize_SuN(SuN *A);

void print_on_screen_SuN(SuN const * const A);
void print_on_file_SuN(FILE *fp, SuN const * const A);
void print_on_binary_file_SuN(FILE *fp, SuN const * const A);
void read_from_file_SuN(FILE *fp, SuN *A);
void read_from_binary_file_SuN(FILE *fp, SuN *A);
void read_from_binary_file_swap_SuN(FILE *fp, SuN *A);


#endif
