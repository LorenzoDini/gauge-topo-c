#ifndef SU2_H
#define SU2_H

#include"../Macro/macro.h"

#include<stdio.h>

typedef struct Su2 {
   #ifdef __GNUC__
     double comp[4] __attribute__ ((aligned (DOUBLE_ALIGN)));
   #else
     double comp[4];
   #endif
} Su2;

void init_Su2(Su2 *__restrict__ A, double vec[4]);
void one_Su2(Su2 *__restrict__ A);         /* A=1 */
void zero_Su2(Su2 *__restrict__ A);        /* A=0 */

void equal_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B);      /* A=B       */
void equal_dag_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B);  /* A=B^{dag} */

void plus_equal_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B);     /* A+=B       */
void plus_equal_dag_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B); /* A+=B^{dag} */

void minus_equal_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B);     /* A-=B       */
void minus_equal_dag_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B); /* A-=B^{dag} */

void lin_comb_Su2(Su2 *__restrict__ A, double b, Su2 const *__restrict__ const B, double c, Su2 const *__restrict__ const C);       /* A=b*B+c*C */
void lin_comb_dag1_Su2(Su2 *__restrict__ A, double b, Su2 const *__restrict__ const B, double c, Su2 const *__restrict__ const C);  /* A=b*B^{dag}+c*C */
void lin_comb_dag2_Su2(Su2 *__restrict__ A, double b, Su2 const *__restrict__ const B, double c, Su2 const *__restrict__ const C);  /* A=b*B+c*C^{dag} */
void lin_comb_dag12_Su2(Su2 *__restrict__ A, double b, Su2 const *__restrict__ const B, double c, Su2 const *__restrict__ const C); /* A=b*B^{dag}+c*C^{dag} */

void times_equal_real_Su2(Su2 *__restrict__ A, double r); /* A*=r */

void times_equal_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B);     /* A*=B       */
void times_equal_dag_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B); /* A*=B^{dag} */

void times_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B, Su2 const *__restrict__ const C);       /* A=B*C             */
void times_dag1_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B, Su2 const *__restrict__ const C);  /* A=B^{dag}*C       */
void times_dag2_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B, Su2 const *__restrict__ const C);  /* A=B*C^{dag}       */
void times_dag12_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B, Su2 const *__restrict__ const C); /* A=B^{dag}*C^{dag} */

void rand_matrix_Su2(Su2 *__restrict__A);
void rand_matrix_p0_Su2(double p0, Su2 *__restrict__A);

double retr_Su2(Su2 const *__restrict__ const A);
double imtr_Su2(Su2 const *__restrict__ const A);
double sqrtdet_Su2(Su2 const *__restrict__ const A);

void unitarize_Su2(Su2 *__restrict__ A);

void print_on_screen_Su2(Su2 const *__restrict__ const A);
void print_on_file_Su2(FILE *fp, Su2 const *__restrict__ const A);
void read_from_file_Su2(FILE *fp, Su2 *__restrict__ A);


#endif
