#ifndef SON_UPD_H
#define SON_UPD_H

#include"../Const/const.h"
#include"son.h"

void randheat_U1(double *link, double const * const staple);
void single_heatbath_SoN(SoN *__restrict__ link,
                         SoN const *__restrict__ const staple, 
                         Const const *__restrict__ const param);
void single_overrelaxation_SoN(SoN *__restrict__ link, 
                               SoN const *__restrict__ const staple, 
                               Const const *__restrict__ const param);
void cool_SoN(SoN *__restrict__ link, 
              SoN const *__restrict__ const staple);

#endif
