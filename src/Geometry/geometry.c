#ifndef GEOMETRY_C
#define GEOMETRY_C

#include<stdio.h>
#include<stdlib.h>

#include"geometry.h"
#include"../Const/const.h"
#include"../Macro/macro.h"


/* lexicographic index */
int lex_index(int t, int x, int y, int z, Const const * const param)
  {
  int ris = t + (param->d_latot)*x + (param->d_latot)*(param->d_latox)*y + (param->d_latot)*(param->d_latox)*(param->d_latoy)*z;
  
  return ris;
  }


void init_geometry(Geometry *geo, Const const * const param)
  {
  int t, tp, tm, x, xp, xm, y, yp, ym, z, zp, zm;
  int ris, risp, rism;
  int i;
  #if HAVE_POSIX_MEMALIGN == 1
    void *vp;
  #endif

  /* allocate memory */  
  #if HAVE_POSIX_MEMALIGN == 1
    t=posix_memalign(&vp, INT_ALIGN, param->d_volume*sizeof(int*));
    if(t!=0)
      {
      fprintf(stderr, "Problems in allocating the geometry!\n");
      }
    geo->d_nnp=vp;
    t=posix_memalign(&vp, INT_ALIGN, param->d_volume*sizeof(int*));
    if(t!=0)
      {
      fprintf(stderr, "Problems in allocating the geometry!\n");
      }
    geo->d_nnm=vp;
  #else
    geo->d_nnp = (int **) malloc(param->d_volume * sizeof(int *));
    if(geo->d_nnp==NULL)
      {
      fprintf(stderr, "Problems in allocating the geometry!\n");
      }
    geo->d_nnm = (int **) malloc(param->d_volume * sizeof(int *));
    if(geo->d_nnm==NULL)
      {
      fprintf(stderr, "Problems in allocating the geometry!\n");
      }
  #endif
  for(i=0; i<(param->d_volume); i++)
     {
     #if HAVE_POSIX_MEMALIGN == 1
       t=posix_memalign(&vp, INT_ALIGN, 4*sizeof(int));
       if(t!=0)
         {
         fprintf(stderr, "Problems in allocating the geometry!\n");
         }
       geo->d_nnp[i]=vp;
       t=posix_memalign(&vp, INT_ALIGN, 4*sizeof(int));
       if(t!=0)
         {
         fprintf(stderr, "Problems in allocating the geometry!\n");
         }
       geo->d_nnm[i]=vp;
     #else
       geo->d_nnp[i] = (int *) malloc(4 * sizeof(int)); 
       if(geo->d_nnp[i]==NULL)
         {
         fprintf(stderr, "Problems in allocating the geometry!\n");
         }
       geo->d_nnm[i] = (int *) malloc(4 * sizeof(int)); 
       if(geo->d_nnm[i]==NULL)
         {
         fprintf(stderr, "Problems in allocating the geometry!\n");
         }
     #endif
     }
  #if HAVE_POSIX_MEMALIGN == 1
     t=posix_memalign(&vp, INT_ALIGN, param->d_volume*sizeof(int));
     if(t!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry!\n");
       }
     geo->d_timeslice=vp;
  #else
     geo->d_timeslice = (int *) malloc(param->d_volume * sizeof(int));
     if(geo->d_timeslice==NULL)
      {
      fprintf(stderr, "Problems in allocating the geometry!\n");
      }
  #endif

  /* assign values */
  for(t=0; t<(param->d_latot); t++)
     {
     if(t==param->d_latot-1) tp=0;
     else tp=t+1;
     if(t==0) tm=param->d_latot-1;
     else tm=t-1;

     for(x=0; x<(param->d_latox); x++)
        {
        if(x==param->d_latox-1) xp=0;
        else xp=x+1;
        if(x==0) xm=param->d_latox-1;
        else xm=x-1;

        for(y=0; y<(param->d_latoy); y++)
           {
           if(y==param->d_latoy-1) yp=0;
           else yp=y+1;
           if(y==0) ym=param->d_latoy-1;
           else ym=y-1;

           for(z=0; z<(param->d_latoz); z++)
              {
              if(z==param->d_latoz-1) zp=0;
              else zp=z+1;
              if(z==0) zm=param->d_latoz-1;
              else zm=z-1;

              ris = lex_index(t, x, y, z, param);
              geo->d_timeslice[ris]=t;          

              rism = lex_index(tm, x, y, z, param);
              risp = lex_index(tp, x, y, z, param); 
              geo->d_nnp[ris][0]=risp;
              geo->d_nnm[ris][0]=rism;

              rism = lex_index(t, xm, y, z, param);
              risp = lex_index(t, xp, y, z, param);
              geo->d_nnp[ris][1]=risp;
              geo->d_nnm[ris][1]=rism;

              rism = lex_index(t, x, ym, z, param);
              risp = lex_index(t, x, yp, z, param);
              geo->d_nnp[ris][2]=risp;
              geo->d_nnm[ris][2]=rism;

              rism = lex_index(t, x, y, zm, param);
              risp = lex_index(t, x, y, zp, param);
              geo->d_nnp[ris][3]=risp;
              geo->d_nnm[ris][3]=rism;
              } 
           }
        }
     }
  }  

/* free memory */
void end_geometry(Geometry *geo, Const const * const param)
  {
  int i;

  for(i=0; i<param->d_volume; i++)
     {
     free(geo->d_nnp[i]);
     free(geo->d_nnm[i]);
     }
  free(geo->d_nnp);
  free(geo->d_nnm);
  free(geo->d_timeslice);
  }


int nnp(Geometry const * const geo, int r, int i)
  {
  return geo->d_nnp[r][i];
  }


int nnm(Geometry const * const geo, int r, int i)
  {
  return geo->d_nnm[r][i];
  }


int timeslice(Geometry const * const geo, int r)
  {
  return geo->d_timeslice[r];
  }


#endif
