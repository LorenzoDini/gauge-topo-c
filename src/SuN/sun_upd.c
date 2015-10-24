#ifndef SUN_UPD_C
#define SUN_UPD_C

#include"../Macro/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>

#include"../GParam/gparam.h"
#include"../Rng/random.h"
#include"../Su2/su2.h"
#include"../Su2/su2_upd.h"
#include"sun.h"
#include"sun_aux.h"

/* Pseudo-heatbath by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
   Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985)) */
void single_heatbath_SuN(SuN *link, SuN const * const staple, GParam const * const param)
    {
    SuN aux, final;
    Su2 u, v, w;
    double xi, p0; 
    double complex temp0, temp1;
    int i, j, k;

    for(i=0; i<Ncolor-1; i++)
       {
       for(j=i+1; j<Ncolor; j++)
          {
          equal_SuN(&aux, staple);     /* aux=staple */
          times_equal_SuN(&aux, link); /*aux=staple*link */
          ennetodue(&aux, i, j, &xi, &u);

          xi*=(param->d_beta)*2.0/((double) Ncolor);

          if(xi>MIN_VALUE)
            {
            randheat_Su2(xi, &p0);

            equal_dag_Su2(&w, &u); /* w=u^{dag} */
            rand_matrix_p0_Su2(p0, &v);
            times_equal_Su2(&w, &v); /* w*=v */
           
            duetoenne(&w, i, j, &final); 

            /* link*=final */
            for(k=0; k<Ncolor; k++)
               {
               temp0=link->comp[k][i]*final.comp[i][i] + link->comp[k][j]*final.comp[j][i];
               temp1=link->comp[k][i]*final.comp[i][j] + link->comp[k][j]*final.comp[j][j];
               link->comp[k][i]=temp0;
               link->comp[k][j]=temp1;
               }
            }
          }
       }
    }


/* Pseudo-overrelaxation by Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) in the implementation by
   Kennedy, Pendleton (Phys. Lett. B 156, p.393 (1985)) */
void single_overrelaxation_SuN(SuN *link, SuN const * const staple, GParam const * const param)
    {
    SuN aux, final;
    Su2 u,v;
    double xi; 
    double complex temp0, temp1;
    int i, j, k;

    for(i=0; i<Ncolor-1; i++)
       {
       for(j=i+1; j<Ncolor; j++)
          { 
          equal_SuN(&aux, staple);     /* aux=staple */
          times_equal_SuN(&aux, link); /*aux=staple*link */
          ennetodue(&aux, i, j, &xi, &u);

          if(xi>MIN_VALUE)
            {
            equal_dag_Su2(&v, &u);       /* v=u^{dag} */
            times_equal_dag_Su2(&v, &u); /* v=(u^{dag})^2 */

            duetoenne(&v, i, j, &final); 

            /*link*=final */
            for(k=0; k<Ncolor; k++)
               {
               temp0=link->comp[k][i]*final.comp[i][i] + link->comp[k][j]*final.comp[j][i];
               temp1=link->comp[k][i]*final.comp[i][j] + link->comp[k][j]*final.comp[j][j];
               link->comp[k][i]=temp0;
               link->comp[k][j]=temp1;
               }
            }
          }
       }
    }


/* cooling */
void cool_SuN(SuN *link, SuN const * const staple)
  {
  SuN prod, auxN;
  Su2 aux2, aux2bis;
  double complex temp[2];
  double xi; /* not used */
  int i, j, k;

  equal_SuN(&prod, staple);         /* prod=staple */
  times_equal_SuN(&prod, link);     /* prod=staple*link */
 
  for(i=0; i<Ncolor-1; i++)
     {
     for(j=i+1; j<Ncolor; j++)
        {
        ennetodue(&prod, i, j, &xi, &aux2);

        equal_dag_Su2(&aux2bis, &aux2);   /* aux2bis = (aux2)^{dag} */

        duetoenne(&aux2bis, i, j, &auxN);

        /* link*=final */
        for(k=0; k<Ncolor; k++)
           {
           temp[0]=link->comp[k][i]*auxN.comp[i][i] + link->comp[k][j]*auxN.comp[j][i];
           temp[1]=link->comp[k][i]*auxN.comp[i][j] + link->comp[k][j]*auxN.comp[j][j];
           link->comp[k][i]=temp[0];
           link->comp[k][j]=temp[1];
           }

        /* prod*=final */
        for(k=0; k<Ncolor; k++)
           {
           temp[0]=prod.comp[k][i]*auxN.comp[i][i] + prod.comp[k][j]*auxN.comp[j][i];
           temp[1]=prod.comp[k][i]*auxN.comp[i][j] + prod.comp[k][j]*auxN.comp[j][j];
           prod.comp[k][i]=temp[0];
           prod.comp[k][j]=temp[1];
           }
        }
     }
  }

#endif
