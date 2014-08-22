#ifndef SUN_UPD_H
#define SUN_UPD_H

#include"../Const/const.h"
#include"../Su2/su2.h"
#include"sun.h"

void single_heatbath_SuN(SuN *__restrict__ link,
                         SuN const *__restrict__ const staple, 
                         Const const *__restrict__ const param);
void single_overrelaxation_SuN(SuN *__restrict__ link, 
                               SuN const *__restrict__ const staple, 
                               Const const *__restrict__ const param);
void cool_SuN(SuN *__restrict__ link, 
              SuN const *__restrict__ const staple);

#endif
