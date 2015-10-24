#ifndef CONST_H
#define CONST_H

#include"../Macro/macro.h"

#include<time.h>

typedef struct GParam {
  /* lattice dimensions */
  int d_latox;
  int d_latoy;
  int d_latoz;
  int d_latot;

  /* simulation parameters */
  double d_beta;
  int d_campione;
  int d_thermal;
  int d_over;
  int d_measevery;

  /* smoothing */
  int d_nummeas;
  int d_cooling;

  /* initializations */
  int d_inizio;
  int d_saveconf;

  /* file names */
  char d_conf_file[50];
  char d_data_file[50];
  char d_err_file[50];

  /* random seed */
  int d_randseed;

  /* derived constants */
  int d_volume;
  double d_inv_vol;
  int d_space_vol;
  double d_inv_space_vol;
} GParam;


int readinput(char *in_file, GParam *param);
void print_parameters(GParam const * const param, time_t time_start, time_t time_end);

#endif
