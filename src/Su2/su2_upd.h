#ifndef SU2_UPD_H
#define SU2_UPD_H

#include"../Macro/macro.h"

#include"../Const/const.h"
#include"su2.h"

void randheat_Su2(double k, double *__restrict__ out);
void single_heatbath_Su2(Su2 *__restrict__ link, Su2 const *__restrict__ const staple, Const const *__restrict__ const param);
void single_overrelaxation_Su2(Su2 *__restrict__ link, Su2 const *__restrict__ const staple, Const const *__restrict__ const param);
void cool_Su2(Su2 *__restrict__ link, Su2 const *__restrict__ const staple);

#endif
