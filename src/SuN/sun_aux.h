#ifndef SUN_AUX_H
#define SUN_AUX_H

#include"sun.h"
#include"../Su2/su2.h"

void ennetodue(SuN const *__restrict__ const in, 
               int i, 
               int j, 
               double *xi, 
               Su2 *__restrict__ u);
void duetoenne(Su2 const *__restrict__ const in, 
               int i, 
               int j, 
               SuN *__restrict__ out);


#endif
