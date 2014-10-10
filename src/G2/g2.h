#ifndef G2_H
#define G2_H

/* for definitions see
   Greensite, Langfeld, Olejnik, Reinhardt, Tok   Phys Rev D 75, p.034501 (2007) */

#include"../Macro/macro.h"

#include<stdio.h>

typedef struct G2 {
   #ifdef __GNUC__
     double comp[7][7] __attribute__ ((aligned (DOUBLE_ALIGN)));
   #else
     double comp[7][7];
   #endif
} G2;

void one_G2(G2 *__restrict__ A);         /* A=1 */
void zero_G2(G2 *__restrict__ A);        /* A=0 */

void equal_G2(G2 *__restrict__ A, G2 const *__restrict__ const B);      /* A=B       */
void equal_dag_G2(G2 *__restrict__ A, G2 const *__restrict__ const B);  /* A=B^{dag} */

void plus_equal_G2(G2 *__restrict__ A, G2 const *__restrict__ const B);     /* A+=B       */
void plus_equal_dag_G2(G2 *__restrict__ A, G2 const *__restrict__ const B); /* A+=B^{dag} */

void minus_equal_G2(G2 *__restrict__ A, G2 const *__restrict__ const B);     /* A-=B       */
void minus_equal_dag_G2(G2 *__restrict__ A, G2 const *__restrict__ const B); /* A-=B^{dag} */

void lin_comb_G2(G2 *__restrict__ A, double b, G2 const *__restrict__ const B, double c, G2 const *__restrict__ const C);       /* A=b*B+c*C */
void lin_comb_dag1_G2(G2 *__restrict__ A, double b, G2 const *__restrict__ const B, double c, G2 const *__restrict__ const C);  /* A=b*B^{dag}+c*C */
void lin_comb_dag2_G2(G2 *__restrict__A, double b, G2 const *__restrict__ const B, double c, G2 const *__restrict__ const C);  /* A=b*B+c*C^{dag} */
void lin_comb_dag12_G2(G2 *__restrict__A, double b, G2 const *__restrict__ const B, double c, G2 const *__restrict__ const C); /* A=b*B^{dag}+c*C^{dag} */

void times_equal_real_G2(G2 *__restrict__ A, double r); /* A*=r */

void times_equal_G2(G2 *__restrict__ A, G2 const *__restrict__ const B);     /* A*=B       */
void times_equal_4_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, int i0, int i1, int i2, int i3);
void times_equal_dag_G2(G2 *__restrict__ A, G2 const *__restrict__ const B); /* A*=B^{dag} */

void times_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, G2 const *__restrict__ const C);       /* A=B*C             */
void times_dag1_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, G2 const *__restrict__ const C);  /* A=B^{dag}*C       */
void times_dag2_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, G2 const *__restrict__ const C);  /* A=B*C^{dag}       */
void times_dag12_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, G2 const *__restrict__ const C); /* A=B^{dag}*C^{dag} */

/* rotation matrices */
void D1(G2 *__restrict__ A, double x);
void D2(G2 *__restrict__ A, double x);
void D3(G2 *__restrict__ A, double x);
void D4(G2 *__restrict__ A, double x);
void D5(G2 *__restrict__ A, double x);
void D6(G2 *__restrict__ A, double x);
void D7(G2 *__restrict__ A, double x);
void D8(G2 *__restrict__ A, double x);
void D9(G2 *__restrict__ A, double x);
void D10(G2 *__restrict__ A, double x);
void D11(G2 *__restrict__ A, double x);
void D12(G2 *__restrict__ A, double x);
void D13(G2 *__restrict__ A, double x);
void D14(G2 *__restrict__ A, double x);

void rand_matrix_G2(G2 *__restrict__ A);

double retr_G2(G2 const *__restrict__ const A);
double imtr_G2(G2 const *__restrict__ const A);
double norm_G2(G2 const *__restrict__ const A);

void print_on_screen_G2(G2 const *__restrict__ const A);
void print_on_file_G2(FILE *fp, G2 const *__restrict__ const A);
void read_from_file_G2(FILE *fp, G2 *__restrict__ A);

#endif
