#ifndef FUNCTION_POINTERS_H
#define FUNCTION_POINTERS_H

#include<stdio.h>

#include"../G2/g2.h"
#include"../G2/g2_check.h"
#include"../G2/g2_un.h"
#include"../G2/g2_upd.h"
#include"../Macro/macro.h"
#include"../SoN/son.h"
#include"../SoN/son_upd.h"
#include"../Su2/su2.h"
#include"../Su2/su2_upd.h"
#include"../SuN/sun.h"
#include"../SuN/sun_upd.h"



void (*one)(GAUGE_GROUP *A);         /* A=1 */
void (*zero)(GAUGE_GROUP *A);        /* A=0 */

void (*equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B);      /* A=B       */
void (*equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B);  /* A=B^{dag} */

void (*plus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B);     /* A+=B       */
void (*plus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B); /* A+=B^{dag} */

void (*minus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B);     /* A-=B       */
void (*minus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B); /* A-=B^{dag} */

void (*lin_comb)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C);       /* A=b*B+c*C */
void (*lin_comb_dag1)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C);  /* A=b*B^{dag}+c*C */
void (*lin_comb_dag2)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C);  /* A=b*B+c*C^{dag} */
void (*lin_comb_dag12)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C); /* A=b*B^{dag}+c*C^{dag} */

void (*times_equal_real)(GAUGE_GROUP *A, double r); /* A*=r */

void (*times_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B);     /* A*=B       */
void (*times_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B); /* A*=B^{dag} */

void (*times)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C);       /* A=B*C             */
void (*times_dag1)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C);  /* A=B^{dag}*C       */
void (*times_dag2)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C);  /* A=B*C^{dag}       */
void (*times_dag12)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C); /* A=B^{dag}*C^{dag} */

void (*rand_matrix)(GAUGE_GROUP *A);

double (*retr)(GAUGE_GROUP const * const A);
double (*imtr)(GAUGE_GROUP const * const A);

void (*unitarize)(GAUGE_GROUP *A);

void (*print_on_screen)(GAUGE_GROUP const * const A);
void (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A);
void (*read_from_file)(FILE *fp, GAUGE_GROUP *A);

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple, Const const * const param);
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple, Const const * const param);
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);



void init_function_pointers(void);

#endif
