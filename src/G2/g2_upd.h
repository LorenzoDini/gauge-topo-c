#ifndef G2_UPD_H
#define G2_UPD_H

#include"../Macro/macro.h"

#include"../Const/const.h"
#include"g2.h"


void single_heatbath_G2(G2 *link, G2 const * const staple, Const const * const param);
void single_overrelaxation_G2(G2 *link, G2 const * const staple, Const const * const param);
void cool_G2(G2 * link, G2 const * const staple);

#endif
