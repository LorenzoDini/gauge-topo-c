#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../Macro/macro.h"

#include"../Const/const.h"
#include"../Func_Point/function_pointers.h"
#include"../G2/g2.h"
#include"gauge_conf.h"
#include"../Rng/random.h"
#include"../Su2/su2.h"
#include"../Su2/su2_upd.h"
#include"../SuN/sun.h"
#include"../SuN/sun_upd.h"


/* compute the staple in position r, direction i and save it in M */
void calcstaples(Gauge_Conf const *__restrict__ const GC, int r, int i, GAUGE_GROUP *__restrict__ M) 
   {
   int j, k;
   GAUGE_GROUP aux;

   zero(M); /* M=0 */

   for(j=0;j<4;j++)
      {
      if(j!=i)
        {
/*
       i ^
         |   (1)
     (b) +----->-----+ (c)
         |           |
         |           |
         |           V (2)
         |           |
         |           |
     (a) +-----<-----+-->   j
       r     (3)   (d)
*/
        equal(&aux, &(GC->lattice[nnp(&(GC->geo), r, i)][j]) );           /* 1 */
        times_equal_dag(&aux, &(GC->lattice[nnp(&(GC->geo), r, j)][i]) ); /* 2 */
        times_equal_dag(&aux, &(GC->lattice[r][j]) );                     /* 3 */

        plus_equal(M, &aux);    /* M+=aux */

/*
       i ^
         |   (1)
     (d) |----<------+ (a)
         |           |
         |           |
     (2) V           |         k=c
         |           | 
         |           | (b)
         +------>----+--->j
      (c)     (3)    r           
*/

        k=nnm(&(GC->geo), r, j);
        equal_dag(&aux, &(GC->lattice[nnp(&(GC->geo), k, i)][j]));  /* 1 */
        times_equal_dag(&aux, &(GC->lattice[k][i]));                /* 2 */
        times_equal(&aux, &(GC->lattice[k][j]));                    /* 3 */

        plus_equal(M, &aux);  /* M+=aux */
        }
      }
    }


/* perform an update with heatbath */
void heatbath(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param, int r, int i)
   {
   GAUGE_GROUP staple;

   calcstaples(GC, r, i, &staple);
   single_heatbath(&(GC->lattice[r][i]), &staple, param); 
   }


/* perform an update with overrelaxation */
void overrelaxation(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param, int r, int i)
   {
   GAUGE_GROUP staple;

   calcstaples(GC, r, i, &staple);
   single_overrelaxation(&(GC->lattice[r][i]), &staple, param);  
   }


/* perform a complete update */
void update(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param)
   {
   int r, j, i;
   #ifdef RAND_GAUGE_TRANSF
   GAUGE_GROUP aux1, aux2;
   #endif

   /* heatbath */
   for(i=0; i<4; i++)
      {
      for(r=0; r<(param->d_volume); r++)
         {
         heatbath(GC, param, r, i);
         } 
      }

   /* overrelax */
   for(j=0; j<param->d_over; j++) 
      {
      for(i=0; i<4; i++)
         {
         for(r=0; r<param->d_volume; r++)
            {
            overrelaxation(GC, param, r, i);
            }
         }
      }
   
   /* random gauge transformation */
   #ifdef RAND_GAUGE_TRANSF
   for(r=0; r<(param->d_volume); r++)
      {
      rand_matrix(&aux1);

      for(i=0; i<4; i++)
         {
         times(&aux2, &aux1, &(GC->lattice[r][i]));   
         /* aux2=aux1*lattice[r][i] */
         equal(&(GC->lattice[r][i]), &aux2);          
         /* lattice[r][i]-> rand * lattice[r][i] */ 

         times_dag2(&aux2, &(GC->lattice[nnm(&(GC->geo),r,i)][i]), &aux1);   
         /* aux2=lattice[r][i] * aux1^dag */
         equal(&(GC->lattice[nnm(&(GC->geo),r,i)][i]), &aux2);               
         /* lattice[nnm(r,i)][i] -> lattice[nnm(r,i)][i] * rand^dag */
         } 
      }
   #endif

   /* final unitarization */
   for(r=0; r<(param->d_volume); r++)
      {
      for(i=0; i<4; i++)
         {
         unitarize(&(GC->lattice[r][i]));
         } 
      }
   }


/* perform n cooling steps */
void cooling(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param, int n)
   {
   GAUGE_GROUP staple;
   #ifdef RAND_GAUGE_TRANSF
   GAUGE_GROUP aux1, aux2;
   #endif
   int r, i, k; 

   for(k=0; k<n; k++)
      {
      /* cooling */
      for(i=0; i<4; i++)
         {
         for(r=0; r<(param->d_volume); r++)
            {
            calcstaples(GC, r, i, &staple);
            cool(&(GC->lattice[r][i]), &staple);  
            }
         }

      /* random gauge transformation */
      #ifdef RAND_GAUGE_TRANSF
      for(r=0; r<(param->d_volume); r++)
         {
         rand_matrix(&aux1);

         for(i=0; i<4; i++)
            {
            times(&aux2, &aux1, &(GC->lattice[r][i]));   
            /* aux2=aux1*lattice[r][i] */
            equal(&(GC->lattice[r][i]), &aux2);          
            /* lattice[r][i]-> rand * lattice[r][i] */ 

            times_dag2(&aux2, &(GC->lattice[nnm(&(GC->geo),r,i)][i]), &aux1);   
            /* aux2=lattice[r][i] * aux1^dag */
            equal(&(GC->lattice[nnm(&(GC->geo),r,i)][i]), &aux2);                
            /* lattice[nnm(r,i)][i] -> lattice[nnm(r,i)][i] * rand^dag */
            } 
         }
      #endif
      }

   /* final unitarization */
   for(r=0; r<(param->d_volume); r++)
      {
      for(i=0; i<4; i++)
         {
         unitarize(&(GC->lattice[r][i]));
         } 
      }
   }

/* perform n cooling steps different method */
void cooling2(Gauge_Conf *__restrict__ GC, Const const *__restrict__ const param, int n)
   {
   GAUGE_GROUP staple;
   #ifdef RAND_GAUGE_TRANSF
   GAUGE_GROUP aux1, aux2;
   #endif
   int r, r1, i, j, k; 

   for(k=0; k<n; k++)
      {
      for(j=2; j<6; j++)
         {
         i=j % 4;
         for(r1=(param->d_volume)/2; r1<3*(param->d_volume)/2; r1++)
            {
            r=r1 % param->d_volume;
            calcstaples(GC, r, i, &staple);
            cool(&(GC->lattice[r][i]), &staple);  
            }
         }

      /* random gauge transformation */
      #ifdef RAND_GAUGE_TRANSF
      for(r=0; r<(param->d_volume); r++)
         {
         rand_matrix(&aux1);

         for(i=0; i<4; i++)
            {
            times(&aux2, &aux1, &(GC->lattice[r][i]));   
            /* aux2=aux1*lattice[r][i] */
            equal(&(GC->lattice[r][i]), &aux2);          
            /* lattice[r][i]-> rand * lattice[r][i] */ 

            times_dag2(&aux2, &(GC->lattice[nnm(&(GC->geo),r,i)][i]), &aux1);   
            /* aux2=lattice[r][i] * aux1^dag */
            equal(&(GC->lattice[nnm(&(GC->geo),r,i)][i]), &aux2);               
            /* lattice[nnm(r,i)][i] -> lattice[nnm(r,i)][i] * rand^dag */
            } 
         }
      #endif
      }

   /* final unitarization */
   for(r=0; r<(param->d_volume); r++)
      {
      for(i=0; i<4; i++)
         {
         unitarize(&(GC->lattice[r][i]));
         } 
      }
   }




/* compute the four-leaf clover in position r, in the plane i,j and save it in M */
void quadrifoglio(Gauge_Conf const *__restrict__ const GC, int r, int j, int i, GAUGE_GROUP *__restrict__ M)
   {
   GAUGE_GROUP aux; 
   int k, p;

   zero(M);

/*
                   i ^
                     |
             (14)    |     (3)
         +-----<-----++-----<-----+
         |           ||           |
         |           ||           |
   (15)  V      (13) ^V (4)       ^ (2)
         |           ||           |
         |   (16)    || r   (1)   |
    p    +----->-----++----->-----+------>   j
         +-----<-----++-----<-----+
         |    (9)    ||   (8)     |
         |           ||           |
    (10) V      (12) ^V (5)       ^ (7) 
         |           ||           |
         |           ||           |
         +------>----++----->-----+
              (11)   k      (6)
*/ 
   /* avanti-avanti */
   equal(&aux, &(GC->lattice[r][j]) );                               /* 1 */
   times_equal(&aux, &(GC->lattice[nnp(&(GC->geo), r, j)][i]) );     /* 2 */
   times_equal_dag(&aux, &(GC->lattice[nnp(&(GC->geo), r, i)][j]) ); /* 3 */
   times_equal_dag(&aux, &(GC->lattice[r][i]) );                     /* 4 */
   plus_equal(M, &aux);
 
   k=nnm(&(GC->geo), r, i);

   /* avanti-indietro */
   equal_dag(&aux, &(GC->lattice[k][i]) );                       /* 5 */
   times_equal(&aux, &(GC->lattice[k][j]) );                     /* 6 */
   times_equal(&aux, &(GC->lattice[nnp(&(GC->geo), k, j)][i]) ); /* 7 */
   times_equal_dag(&aux, &(GC->lattice[r][j]) );                 /* 8 */
   plus_equal(M, &aux);

   p=nnm(&(GC->geo), r, j);

   /* indietro-indietro */
   
   equal_dag(&aux, &(GC->lattice[p][j]) );                           /* 9 */
   times_equal_dag(&aux, &(GC->lattice[nnm(&(GC->geo), k, j)][i]) ); /* 10 */
   times_equal(&aux, &(GC->lattice[nnm(&(GC->geo), k, j)][j]) );     /* 11 */
   times_equal(&aux, &(GC->lattice[k][i]) );                         /* 12 */
   plus_equal(M, &aux);

   /* indietro-avanti */
   equal(&aux, &(GC->lattice[r][i]) );                                /* 13 */
   times_equal_dag(&aux, &(GC->lattice[nnp(&(GC->geo), p, i)][j]) );  /* 14 */
   times_equal_dag(&aux, &(GC->lattice[p][i]) );                      /* 15 */
   times_equal(&aux, &(GC->lattice[p][j]) );                          /* 16 */
   plus_equal(M, &aux);
   }

#endif
