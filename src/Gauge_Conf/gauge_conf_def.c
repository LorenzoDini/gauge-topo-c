#ifndef GAUGE_CONF_DEF_C
#define GAUGE_CONF_DEF_C

#include<openssl/md5.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../Const/const.h"
#include"../Geometry/geometry.h"
#include"gauge_conf.h"
#include"../Func_Point/function_pointers.h"
#include"../Macro/macro.h"
#include"../Su2/su2.h"

int init_gauge_conf(Gauge_Conf * GC, Const const * const param)
  {
  int i, j;

  /* allocate lattice */
  GC->lattice = (GAUGE_GROUP **) malloc(param->d_volume * sizeof(GAUGE_GROUP *)); 
  if(GC->lattice == NULL)
    {
    fprintf(stderr, "Problems in allocating the lattice!\n");
    return 1;    
    }
  for(i=0; i<(param->d_volume); i++)
     {
     GC->lattice[i] = (GAUGE_GROUP *) malloc( 4 * sizeof(GAUGE_GROUP)); 
     if(GC->lattice[i]==NULL)
       {
       fprintf(stderr, "Problems in allocating the lattice!\n");
       return 1;    
       }
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
    char md5sum_new[2*MD5_DIGEST_LENGTH+1]; 
    char md5sum_old[2*MD5_DIGEST_LENGTH+1];
    double bl;

    fp=fopen(param->conf_file, "r"); /* open the configuration file */
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", param->conf_file);
      return 1;
      }
    else
      {
      i=fscanf(fp, "%d %d %d %d %lg %s\n", &xl, &yl, &zl, &tl, &bl, md5sum_old);
      if(i!=6)
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
      fclose(fp);
      }

    fp=fopen(param->conf_file, "rb"); /* open the configuration file in binary*/
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", param->conf_file);
      return 1;
      }
    else
      {
      // read again the header: xl, yl, zl, tl, bl and hash 
      i=0;
      while(i!='\n')
           {
           i=fgetc(fp);
           }

      // read the configuration
      for(i=0; i<(param->d_volume); i++)
         {
         for(j=0; j<4; j++)
            {
            read_from_binary_file(fp, &(GC->lattice[i][j]));
            }
         }
      fclose(fp);

      // compute the new md5sum
      compute_md5sum(md5sum_new, GC, param);
      if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0) 
        {
        fprintf(stderr, "The computed md5sum %s does not match the stored %s\n", md5sum_new, md5sum_old);
        return 1;
        }
      } // closure of if(param->d_inizio==2)
    }

  return 0;
  }


void end_gauge_conf(Gauge_Conf * GC, Const const * const param)
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


void save_on_file(Gauge_Conf const * const GC, Const const * const param)
  {
  int i, j;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum(md5sum, GC, param);

  fp=fopen(param->conf_file, "w"); /* open the configuration file */
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s\n", param->conf_file);
    }
  else
    {
    fprintf(fp, "%d %d %d %d %.10g %s\n", 
                param->d_latox, 
                param->d_latoy, 
                param->d_latoz, 
                param->d_latot, 
                param->d_beta, 
                md5sum);
    }
  fclose(fp);

  fp=fopen(param->conf_file, "ab"); /* open the configuration file in binary mode*/
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s\n", param->conf_file);
    }
  else
    {
    for(i=0; i<(param->d_volume); i++)
       {
       for(j=0; j<4; j++)
          {
          print_on_binary_file(fp, &(GC->lattice[i][j]));
          }
       }
    fclose(fp);
    }
  }


/* allocate GC and initialize with GC2 */
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC, Gauge_Conf const * const GC2, Const const * const param) 
  {
  int i, j;

  /* allocate lattice */
  GC->lattice = (GAUGE_GROUP **) malloc(param->d_volume * sizeof(GAUGE_GROUP *)); 
  if(GC->lattice == NULL)
    {
    fprintf(stderr, "Problems in allocating the lattice!\n");
    }

  for(i=0; i<(param->d_volume); i++)
     {
     GC->lattice[i] = (GAUGE_GROUP *) malloc( 4 * sizeof(GAUGE_GROUP)); 
     if(GC->lattice[i] == NULL)
       {
       fprintf(stderr, "Problems in allocating the lattice!\n");
       }
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



/* compute the md5sum of the configuration and save it in res, char [2*MD5_DIGEST_LENGTH] */
void compute_md5sum(char *res, Gauge_Conf const * const GC, Const const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  int i, j;
  int bytes=sizeof(GAUGE_GROUP);

  MD5_Init(&mdContext);
  for(i=0; i<(param->d_volume); i++)
     {
     for(j=0; j<4; j++)
        {
        MD5_Update(&mdContext, &(GC->lattice[i][j]), bytes);
        }
     }
  MD5_Final(c, &mdContext);
  
  for(i = 0; i < MD5_DIGEST_LENGTH; i++)
     {
     sprintf(&(res[2*i]), "%02x", c[i]);
     }
  }



#endif
