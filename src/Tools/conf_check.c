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
    char *infile;
    FILE *fp;

    int xl, yl, zl, tl, volumel;
    double bl; 
    int i, j;

    char md5sum_old[2*MD5_DIGEST_LENGTH+1]; 
    char md5sum_new[2*MD5_DIGEST_LENGTH+1]; 

    GAUGE_GROUP link;

    MD5_CTX mdContext;
    unsigned char c[MD5_DIGEST_LENGTH];
    int bytes=sizeof(GAUGE_GROUP);


    if(argc != 2)
      {
      printf("Usage: %s conf_file\n", argv[0]);
      return 0;
      }
    else
      {
      infile=argv[1];
      }

    /* initialize function_pointers */ 
    init_function_pointers();

    /* open the configuration file to read header */
    fp=fopen(infile, "r"); 
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s\n", infile);
      return 1;
      }
    else
      {
      /* xl, yl, zl, tl = lattice sizes */
      /* bl = beta */
      i=fscanf(fp, "%d %d %d %d %lg %s", &xl, &yl, &zl, &tl, &bl, md5sum_old);
      if(i!=6)
        {
        fprintf(stderr, "Error in reading the file %s\n", infile);
        return 1;
        }

      fclose(fp);
      }

    /* total volume */
    volumel=xl*yl*zl*tl;

    /* open the configuration file in binary*/
    fp=fopen(infile, "rb");
    if(fp==NULL)
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
           i=fgetc(fp);
           }

      /* read the configuration & compute md5sum */
      MD5_Init(&mdContext);
      for(i=0; i<volumel; i++)
         {
         for(j=0; j<4; j++)
            {
            read_from_binary_file(fp, &link);
            MD5_Update(&mdContext, &link, bytes);
            }
         }
      MD5_Final(c, &mdContext);
      for(i = 0; i < MD5_DIGEST_LENGTH; i++)
         {
         sprintf(&(md5sum_new[2*i]), "%02x", c[i]);
         }

      fclose(fp);
      }

    /* check md5sum computed and stored */
    if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)==0) 
      {
      printf("The configuration %s is OK!\n", infile);
      }
    else
      {
      printf("The configuration %s is CORRUPTED!\n", infile);
      }

    return 0;
    }

#endif
