#ifndef RANDOM_H
#define RANDOM_H

#include"../Macro/macro.h"

double casuale(void);           /* random number in (0,1) */
void initrand(unsigned long s); /* initialize random generator */
double gauss1();                /* normal gaussian random number */

#endif
