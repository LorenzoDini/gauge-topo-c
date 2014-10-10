#ifndef GAUGE_CONF_DEF_C
#define GAUGE_CONF_DEF_C

#include"../Macro/macro.h"

#include<stdio.h>
#include<stdlib.h>

#include"../Const/const.h"
#include"../Geometry/geometry.h"
#include"gauge_conf.h"
#include"../Func_Point/function_pointers.h"
#include"../Su2/su2.h"

int init_gauge_conf(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param)
  {
  int i, j;
  #if HAVE_POSIX_MEMALIGN == 1 
  void *vp;
  #endif

  /* allocate lattice */
  #if HAVE_POSIX_MEMALIGN == 1 
    j=posix_memalign(&vp, DOUBLE_ALIGN, param->d_volume * sizeof(GAUGE_GROUP *));
    if(j!=0)
      {
      fprintf(stderr, "Problems in allocating the lattice!\n");
      return 1;    
      }
    GC->lattice=vp;
  #else
    GC->lattice = (GAUGE_GROUP **) malloc(param->d_volume * sizeof(GAUGE_GROUP *)); 
    if(GC->lattice == NULL)
      {
      fprintf(stderr, "Problems in allocating the lattice!\n");
      return 1;    
      }
  #endif
  for(i=0; i<(param->d_volume); i++)
     {
     #if HAVE_POSIX_MEMALIGN == 1 
       j=posix_memalign(&vp, DOUBLE_ALIGN, 4 * sizeof(GAUGE_GROUP));
       if(j!=0)
         {
         fprintf(stderr, "Problems in allocating the lattice!\n");
         return 1;    
         }
       GC->lattice[i]=vp;
     #else
       GC->lattice[i] = (GAUGE_GROUP *) malloc( 4 * sizeof(GAUGE_GROUP)); 
       if(GC->lattice[i]==NULL)
         {
         fprintf(stderr, "Problems in allocating the lattice!\n");
         return 1;    
         }
     #endif
     }

  /* initialize geometry */
  init_geometry(&(GC->geo), param);

  /* initialize GC */
  if(param->d_inizio==0)
    {
    GAUGE_GROUP aux1, aux2;
    one(&aux1);

    for(i=0; i<(param->d_volume); i++)
       {
       for(j=0; j<4; j++)
          {
          rand_matrix(&aux2);
          times_equal_real(&aux2, 0.001);
          plus_equal(&aux2, &aux1);
          unitarize(&aux2); 
          equal(&(GC->lattice[i][j]), &aux2);
          }
       }
    }
  if(param->d_inizio==1)
    {
    GAUGE_GROUP aux1;
    for(i=0; i<(param->d_volume); i++)
       {
       for(j=0; j<4; j++)
          {
          rand_matrix(&aux1);
          equal(&(GC->lattice[i][j]), &aux1);
          }
       }
    }
  if(param->d_inizio==2)
    {
    FILE *fp;
    int xl, yl, zl, tl, i;
    double bl;

    fp=fopen(param->conf_file, "r"); /* open the configuration file */
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", param->conf_file);
      return 1;
      }
    else
      {
      i=fscanf(fp, "%d %d %d %d %lg\n", &xl, &yl, &zl, &tl, &bl);
      if(i!=5)
        {
        fprintf(stderr, "Error in reading the file %s\n", param->conf_file);
        return 1;
        }

      if(xl!=param->d_latox || yl!=param->d_latoy || zl!=param->d_latoz || tl!=param->d_latot)
        {
        fclose(fp);
        fprintf(stderr, "The configuration in %s is not of the correct size!\n", param->conf_file);
        return 1;
        }      

      for(i=0; i<(param->d_volume); i++)
         {
         for(j=0; j<4; j++)
            {
            read_from_file(fp, &(GC->lattice[i][j]));
            }
         }
      fclose(fp);
      }
    }

  return 0;
  }


void end_gauge_conf(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param)
  {
  int i;

  /* free lattice */
  for(i=0; i<(param->d_volume); i++)
     {
     free(GC->lattice[i]);
     }
  free(GC->lattice);

  /* free geometry */
  end_geometry(&(GC->geo), param);
  }


void save_on_file(Gauge_Conf const *__restrict__ const GC, Const const *__restrict__ const param)
  {
  int i, j;
  FILE *fp;

  fp=fopen(param->conf_file, "w"); /* open the configuration file */
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s\n", param->conf_file);
    }
  else
    {
    fprintf(fp, "%d %d %d %d %.10g\n", param->d_latox, param->d_latoy, param->d_latoz, param->d_latot, param->d_beta);
    for(i=0; i<(param->d_volume); i++)
       {
       for(j=0; j<4; j++)
          {
          print_on_file(fp, &(GC->lattice[i][j]));
          }
       }
    fclose(fp);
    }
  }


/* allocate GC and initialize with GC2 */
void init_gauge_conf_from_gauge_conf(Gauge_Conf *__restrict__ GC, Gauge_Conf const *__restrict__ const GC2, Const const *__restrict__ const param) 
  {
  int i, j;
  #if HAVE_POSIX_MEMALIGN==1
  void *vp;
  #endif

  /* allocate lattice */
  #if HAVE_POSIX_MEMALIGN==1
    j=posix_memalign(&vp, DOUBLE_ALIGN, param->d_volume * sizeof(GAUGE_GROUP *));
    if(j!=0)
      {
      fprintf(stderr, "Problems in allocating the lattice!\n");
      }
    GC->lattice=vp;
  #else
    GC->lattice = (GAUGE_GROUP **) malloc(param->d_volume * sizeof(GAUGE_GROUP *)); 
    if(GC->lattice == NULL)
      {
      fprintf(stderr, "Problems in allocating the lattice!\n");
      }
  #endif

  for(i=0; i<(param->d_volume); i++)
     {
     #if HAVE_POSIX_MEMALIGN==1
       j=posix_memalign(&vp, DOUBLE_ALIGN, 4 * sizeof(GAUGE_GROUP));
       if(j!=0)
         {
         fprintf(stderr, "Problems in allocating the lattice!\n");
         }
       GC->lattice[i]=vp;
     #else
       GC->lattice[i] = (GAUGE_GROUP *) malloc( 4 * sizeof(GAUGE_GROUP)); 
      if(GC->lattice[i] == NULL)
        {
        fprintf(stderr, "Problems in allocating the lattice!\n");
        }
     #endif

     }

  /* initialize geometry */
  init_geometry(&(GC->geo), param);

  /* initialize GC */
  for(i=0; i<(param->d_volume); i++)
     {
     for(j=0; j<4; j++)
        {
        equal(&(GC->lattice[i][j]), &(GC2->lattice[i][j]) );
        }
     }
  }





#endif
