#ifndef SUN_UPD_H
#define SUN_UPD_H

#include"../Const/const.h"
#include"../Su2/su2.h"
#include"sun.h"

void single_heatbath_SuN(SuN *link, SuN const * const staple, Const const * const param);
void single_overrelaxation_SuN(SuN *link, SuN const * const staple, Const const * const param);
void cool_SuN(SuN *link, SuN const * const staple);

#endif
