#ifndef GAUGE_TOPO_H
#define GAUGE_TOPO_H

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include"gauge_topo_include.h"

int main (int argc, char **argv)
    {
    int i, count;
    char *in_file;
    FILE *datafilep;
    Const param;
    Gauge_Conf GC;
    double plaqs, plaqt, polyre, polyim, *charge, *meanplaq, charge_nocooling;
    time_t time1, time2;

    if(argc != 2)
      {
      printf("Program %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compiled with:\n");
      printf("  gauge group: %s\n", QUOTEME(GAUGE_GROUP));
      printf("  number of colors: %d\n", NCOLOR);
      return 0;
      }
    else
      {
      in_file=argv[1];
      }

    /* read input file */
    readinput(in_file, &param);    

    /* initialize random generator */
    initrand(param.d_randseed); 

    /* initialize function_pointers */ 
    init_function_pointers();

    /* open data_file */
    if(param.d_inizio==2)
      {
      datafilep=fopen(param.data_file, "r");
      if(datafilep!=NULL)  /* file exists */
        {
        fclose(datafilep);
        datafilep=fopen(param.data_file, "a");
        }
      else
        {
        datafilep=fopen(param.data_file, "w");
        fprintf(datafilep, "%d %d %d %d %.10f\n", param.d_latox, 
                                                  param.d_latoy, 
                                                  param.d_latoz, 
                                                  param.d_latot,  
                                                  param.d_beta);
        }
      }
    else
      {
      datafilep=fopen(param.data_file, "w");
      fprintf(datafilep, "%d %d %d %d %.10f\n", param.d_latox, 
                                                param.d_latoy, 
                                                param.d_latoz, 
                                                param.d_latot,  
                                                param.d_beta);
      }
    fflush(datafilep);

    /* initialize gauge configuration */
    i=init_gauge_conf(&GC, &param);
    if(i!=0)
      {
      return 1;
      }

    charge = (double *) malloc( (param.d_nummeas)*sizeof(double));
    meanplaq = (double *) malloc( (param.d_nummeas)*sizeof(double));

    /* montecarlo */
    time(&time1);
    for(count=0; count < param.d_campione; count++)
       {
       update(&GC, &param);

       if(count % param.d_measevery ==0 && count >= param.d_thermal)
         {
         plaquette(&GC, &param, &plaqs, &plaqt);
         polyakow(&GC, &param, &polyre, &polyim);
         charge_nocooling=topcharge(&GC, &param);

         fprintf(datafilep, "%.12f %.12f %.12f %.12f %.12f ", plaqs, 
                                                              plaqt,  
                                                              polyre, 
                                                              polyim,  
                                                              charge_nocooling);

         topcharge_cooling(&GC, &param, charge, meanplaq); 
         for(i=0; i<param.d_nummeas; i++)
            {
            fprintf(datafilep, "%.12f %.12f ", charge[i], meanplaq[i]); 
            } 
         fprintf(datafilep, "\n"); 

         fflush(datafilep);
         }
       }
    time(&time2);
    /* montecarlo end */

    free(charge);
    free(meanplaq);

    /* close data file */
    fclose(datafilep);

    /* save configuration */
    if(param.d_saveconf!=0)
      {
      save_on_file(&GC, &param);
      }

    /* print simulation details */
    print_parameters(&param, time1, time2);

    /* free gauge configuration */
    end_gauge_conf(&GC, &param);

    return 0;
    }

#endif
