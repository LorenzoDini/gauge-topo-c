#ifndef SUN_AUX_H
#define SUN_AUX_H

#include"sun.h"
#include"../Su2/su2.h"

void ennetodue(SuN const * const in, int i, int j, double *xi, Su2 *u);
void duetoenne(Su2 const * const in, int i, int j, SuN *out);


#endif
