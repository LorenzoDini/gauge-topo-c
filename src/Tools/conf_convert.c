#ifndef CONF_CHECK_C
#define CONF_CHECK_C

#include<openssl/md5.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../Func_Point/function_pointers.h"
#include"../Macro/macro.h"

int main (int argc, char **argv)
    {
    char *infile, *outfile;
    FILE *fpi, *fpo;

    int xl, yl, zl, tl, volumel;
    double bl; 
    int i, j;

    char md5sum_old[2*MD5_DIGEST_LENGTH+1]; 
    char md5sum_new[2*MD5_DIGEST_LENGTH+1]; 

    GAUGE_GROUP link;

    MD5_CTX mdContext;
    unsigned char c[MD5_DIGEST_LENGTH];
    int bytes=sizeof(GAUGE_GROUP);


    if(argc != 3)
      {
      printf("Usage: %s input_conf converted_conf\n", argv[0]);
      return 0;
      }
    else
      {
      infile=argv[1];
      outfile=argv[2];
      }


    /* initialize function_pointers */ 
    init_function_pointers();


    /* open the input configuration file to read header */
    fpi=fopen(infile, "r"); 
    if(fpi==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", infile);
      return 1;
      }
    else
      {
      /* xl, yl, zl, tl = lattice sizes */
      /* bl = beta */
      i=fscanf(fpi, "%d %d %d %d %lg %s", &xl, &yl, &zl, &tl, &bl, md5sum_old);
      if(i!=6)
        {
        fprintf(stderr, "Error in reading the file %s\n", infile);
        return 1;
        }

      fclose(fpi);
      }


    /* total volume */
    volumel=xl*yl*zl*tl;


    /* open the input configuration file in binary*/
    fpi=fopen(infile, "rb");
    if(fpi==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", infile);
      return 1;
      }
    else
      {
      /* read again the header: xl, yl, zl, tl, bl and hash  */
      i=0;
      while(i!='\n')
           {
           i=fgetc(fpi);
           }

      /* read links & compute md5sum */
      MD5_Init(&mdContext);
      for(i=0; i<volumel; i++)
         {
         for(j=0; j<4; j++)
            {
            read_from_binary_file_swap(fpi, &link);
            MD5_Update(&mdContext, &link, bytes);
            }
         }
      MD5_Final(c, &mdContext);
      for(i = 0; i < MD5_DIGEST_LENGTH; i++)
         {
         sprintf(&(md5sum_new[2*i]), "%02x", c[i]);
         }

      fclose(fpi);
      }


    /* write header of the new configuration */
    fpo=fopen(outfile, "w"); 
    if(fpo==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", outfile);
      }
    else
      {
      fprintf(fpo, "%d %d %d %d %.10g %s\n", xl, yl, zl, tl, bl, md5sum_new);
      fclose(fpo);
      }


    /* open the output configuration file in binary mode*/
    fpo=fopen(outfile, "ab");    
    if(fpo==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", outfile);
      }
    else
      {
      /* open the input configuration file in binary*/
      fpi=fopen(infile, "rb");
      if(fpi==NULL)
        {
        fprintf(stderr, "Error in opening the file %s\n", infile);
        return 1;
        }
      else
        {
        /* read again the header: xl, yl, zl, tl, bl and hash  */
        i=0;
        while(i!='\n')
             {
             i=fgetc(fpi);
             }

        /* read links & write links */
        for(i=0; i<volumel; i++)
           {
           for(j=0; j<4; j++)
              {
              read_from_binary_file_swap(fpi, &link);
              print_on_binary_file(fpo, &link); 
              }
           }

        fclose(fpi);
        fclose(fpo);
        }
      }

    return 0;
    }

#endif
