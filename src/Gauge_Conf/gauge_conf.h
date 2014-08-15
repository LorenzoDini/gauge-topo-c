#ifndef GAUGE_CONF_H
#define GAUGE_CONF_H

#include<stdio.h>

#include"../Const/const.h"
#include"../G2/g2.h"
#include"../Geometry/geometry.h"
#include"../Macro/macro.h"
#include"../Su2/su2.h"

typedef struct Gauge_Conf {
  Geometry geo;
  GAUGE_GROUP **lattice;
  } Gauge_Conf;

/* in gauge_conf_def.c */
int init_gauge_conf(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param);
void end_gauge_conf(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param);
void save_on_file(Gauge_Conf const *__restrict__ const GC, Const const *__restrict__ const param);
void init_gauge_conf_from_gauge_conf(Gauge_Conf *__restrict__ GC1, 
                                     Gauge_Conf const *__restrict__ const GC2, 
                                     Const const *__restrict__ const param); /* GC1=GC2 */

/* in gauge_conf_meas.c */
double plaquettep(Gauge_Conf const *__restrict__ const GC, long int r, int i, int j);
void plaquette(Gauge_Conf const *__restrict__ const GC, 
               Const const *__restrict__ const param, 
               double *__restrict__ plaqs, 
               double *__restrict__ plaqt);
void polyakow(Gauge_Conf const *__restrict__ const GC, 
              Const const *__restrict__ const param, 
              double *__restrict__ repoly, 
              double *__restrict__impoly);
double topcharge(Gauge_Conf const *__restrict__ const GC, Const const *__restrict__ const param);
double topchargedens(Gauge_Conf const *__restrict__ const GC, int r);
void topcharge_cooling(Gauge_Conf const *__restrict__ const GC, 
                       Const const *__restrict__ const param, 
                       double *__restrict__ charge, 
                       double *__restrict__ meanplaq); 
void topcharge_cooling2(Gauge_Conf const *__restrict__ const GC, 
                        Const const * const param, 
                        double *__restrict__ charge, 
                        double *__restrict__ meanplaq);


/* in gauge_conf_upd.c */
void calcstaples(Gauge_Conf const *__restrict__ const GC, int r, int i, GAUGE_GROUP *__restrict__ M);
void heatbath(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param, int r, int i);
void overrelaxation(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param, int r, int i);
void update(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param);
void cooling(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param, int n);
void cooling2(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param, int n);
void quadrifoglio(Gauge_Conf const *__restrict__ const GC, int r, int j, int i, GAUGE_GROUP *__restrict__ M);




#endif
