#ifndef GEOMETRY_H
#define GEOMETRY_H

#include"../Macro/macro.h"

#include"../GParam/gparam.h"

int lex_index(int t, int x, int y, int z, GParam const * const param);

typedef struct Geometry {
   int **d_nnp;
   int **d_nnm;
   int *d_timeslice;  
} Geometry;

void init_geometry(Geometry *geo, GParam const * const param);
void end_geometry(Geometry *geo, GParam const * const param);
int nnp(Geometry const * const geo, int r, int i);
int nnm(Geometry const * const geo, int r, int i);
int timeslice(Geometry const * const geo, int r);


#endif
