#ifndef CONST_C
#define CONST_C

#include<stdio.h>
#include<string.h>
#include<time.h>

#include"const.h"
#include"../Macro/macro.h"

int readinput(char *in_file, Const *param)
    {
    FILE *input;
    char str[20], temp_str[20];
    char c;
    double temp_d;
    int temp_i, err, end=1;
    
    input=fopen(in_file, "r");  /* open the input file */
    if(input==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", in_file);
      return 1;
      }
    else
      {
      while(end==1)   /* slide the file */
           {
           err=fscanf(input, "%s", str);
           if(err!=1)
             {
             fprintf(stderr, "Error in reading the file %s\n", in_file);
             printf("err=%d", err);
             return 1;
             }

           if(strncmp(str, "latox", 5)==0)
             { 
             err=fscanf(input, "%d%c", &temp_i, &c);
             if(err!=2)
               {
               fprintf(stderr, "Error in reading the file %s\n", in_file);
               return 1;
               }
             param->d_latox=temp_i;
             /* printf("%d\n", temp_i); */
             }
           else if(strncmp(str, "latoy", 5)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_latoy=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "latoz", 5)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_latoz=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "latot", 5)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_latot=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "beta", 4)==0)
                  { 
                  err=fscanf(input, "%lf%c", &temp_d, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_beta=temp_d;
                  /* printf("%g\n", temp_d); */
                  }
           else if(strncmp(str, "campione", 8)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_campione=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "thermal", 7)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_thermal=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "over", 4)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_over=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "measevery", 9)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_measevery=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "nummeas", 7)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_nummeas=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "cooling", 7)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_cooling=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "inizio", 6)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_inizio=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "saveconf", 8)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_saveconf=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else if(strncmp(str, "conf_file", 9)==0)
                  { 
                  err=fscanf(input, "%s%c", temp_str, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  strcpy(param->conf_file, temp_str);
                  /* printf("%s\n", param->conf_file); */
                  }
           else if(strncmp(str, "data_file", 9)==0)
                  { 
                  err=fscanf(input, "%s%c", temp_str, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  strcpy(param->data_file, temp_str);
                  /* printf("%s\n", temp_str); */
                  }
           else if(strncmp(str, "err_file", 8)==0)
                  { 
                  err=fscanf(input, "%s%c", temp_str, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  strcpy(param->err_file, temp_str);
                  /* printf("%s\n", temp_str); */
                  }
           else if(strncmp(str, "randseed", 8)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s\n", in_file);
                    return 1;
                    }
                  param->d_randseed=temp_i;
                  /* printf("%d\n", temp_i); */
                  }
           else
             {
             fprintf(stderr, "Error: unrecognized option %s in the file %s\n", str, in_file);
             return 1;     
             }
   
           /* discard eventual comments */         
           if(c!='\n')
             { 
             do
               {
               c=getc(input);
               }
             while(c!='\n');
             }

           /* check if the read line is tha last one */
           c=getc(input);
           if(c==EOF)
             {
             end=0;
             }
           else
             {
             ungetc(c, input);
             }
           }

      fclose(input);
      
      /* derived constants */
      param->d_volume=(param->d_latox)*(param->d_latoy)*(param->d_latoz)*(param->d_latot);
      param->d_space_vol=(param->d_latox)*(param->d_latoy)*(param->d_latoz);
      param->d_inv_vol=1.0/((double) param->d_volume);
      param->d_inv_space_vol=1.0/((double) param->d_space_vol);

      return 0;
      }
    }

void print_parameters(Const const * const param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    double diff_sec;

    fp=fopen("dati.log", "w");
    fprintf(fp, "+--------------------+\n");
    fprintf(fp, "| Simulation details |\n");
    fprintf(fp, "+--------------------+\n\n");

    #ifdef ONE_FILE_MODE
      fprintf(fp, "compiled in the single file mode\n");
    #endif  
    fprintf(fp, "precision: double\n");
    fprintf(fp, "lattice: %dx%dx%dx%d\n", param->d_latox, param->d_latoy, param->d_latoz, param->d_latot);
    fprintf(fp, "beta: %g\n", param->d_beta);
    fprintf(fp, "cooling: %d   nummeas: %d\n", param->d_cooling, param->d_nummeas);
    fprintf(fp, "randseed: %d\n", param->d_randseed);
    fprintf(fp, "inizio: %d\n", param->d_inizio);
    fprintf(fp, "saveconf: %d\n", param->d_saveconf);
    fprintf(fp, "campione: %d\n", param->d_campione);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "over: %d\n", param->d_over);
    fprintf(fp, "thermal: %d\n\n", param->d_thermal);

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3f seconds\n\n", diff_sec );

    fclose(fp);
    }


#endif

