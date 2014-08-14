#ifndef CONST_H
#define CONST_H

#include<time.h>

#include"../Macro/macro.h"

typedef struct Const {
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
  char conf_file[20];       
  char data_file[20];       
  char err_file[20];        

  /* random seed */
  int d_randseed;

  /* derived constants */
  int d_volume;
  double d_inv_vol;
  int d_space_vol;
  double d_inv_space_vol;
} Const;


int readinput(char *in_file, Const *param);
void print_parameters(Const const * const param, time_t time_start, time_t time_end);

#endif
