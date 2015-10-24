#ifndef GAUGE_CONF_H
#define GAUGE_CONF_H

#include"../Macro/macro.h"

#include<stdio.h>

#include"../GParam/gparam.h"
#include"../G2/g2.h"
#include"../Geometry/geometry.h"
#include"../Su2/su2.h"
#include"../SuN/sun.h"

typedef struct Gauge_Conf {
  Geometry geo;
  GAUGE_GROUP **lattice;
  } Gauge_Conf;

/* in gauge_conf_def.c */
int init_gauge_conf(Gauge_Conf *, GParam const * const param);
void end_gauge_conf(Gauge_Conf *GC, GParam const * const param);
void save_on_file(Gauge_Conf const * const GC, GParam const * const param);
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC1, 
                                     Gauge_Conf const * const GC2, 
                                     GParam const * const param); /* GC1=GC2 */
void compute_md5sum(char *res, Gauge_Conf const * const GC, GParam const * const param);

/* in gauge_conf_meas.c */
double plaquettep(Gauge_Conf const * const GC, long int r, int i, int j);
void plaquette(Gauge_Conf const * const GC, GParam const * const param, double *plaqs, double *plaqt);
void polyakow(Gauge_Conf const * const GC, GParam const * const param, double *repoly, double *impoly);
double topcharge(Gauge_Conf const * const GC, GParam const * const param);
double topchargedens(Gauge_Conf const * const GC, int r);
void topcharge_cooling(Gauge_Conf const * const GC, GParam const * const param, double *charge, double *meanplaq); 
void topcharge_cooling2(Gauge_Conf const * const GC, GParam const * const param, double * charge, double * meanplaq);


/* in gauge_conf_upd.c */
void calcstaples(Gauge_Conf const * const GC, int r, int i, GAUGE_GROUP *M);
void heatbath(Gauge_Conf *GC, GParam const * const param, int r, int i);
void overrelaxation(Gauge_Conf *GC, GParam const * const param, int r, int i);
void update(Gauge_Conf *GC, GParam const * const param);
void cooling(Gauge_Conf *GC, GParam const * const param, int n);
void cooling2(Gauge_Conf *GC, GParam const * const param, int n);
void quadrifoglio(Gauge_Conf const * const GC, int r, int j, int i, GAUGE_GROUP *M);




#endif
