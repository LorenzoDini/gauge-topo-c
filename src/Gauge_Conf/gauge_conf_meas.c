#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include<math.h>

#include"../Const/const.h"
#include"../Func_Point/function_pointers.h"
#include"gauge_conf.h"
#include"../Macro/macro.h"


/* computation of the plaquette (the trace of) in position r and positive directions i,j  */
double plaquettep(Gauge_Conf const *__restrict__ const GC, long int r, int i, int j)
   {
   GAUGE_GROUP M;
/*
       ^ i
       |   (2)
       +---<---+
       |       |
   (3) V       ^ (1) 
       |       |
       +--->---+---> j
       r   (4) 
*/

   equal(&M, &(GC->lattice[nnp(&(GC->geo), r, j)][i]) );           /* 1 */
   times_equal_dag(&M, &(GC->lattice[nnp(&(GC->geo), r, i)][j]) ); /* 2 */
   times_equal_dag(&M, &(GC->lattice[r][i]) );                     /* 3 */
   times_equal(&M, &(GC->lattice[r][j]) );                         /* 3 */

   return retr(&M);
   }


/* compute the mean plaquettes (spatial, temporal) */
void plaquette(Gauge_Conf const *__restrict__ const GC, 
               Const const *__restrict__ const param, 
               double *__restrict__ plaqs, 
               double *__restrict__ plaqt) 
   {
   int i, j, r;
   double ps=0.0, pt=0.0;

   for(r=0; r<(param->d_volume); r++)
      {
      i=0;
      for(j=1; j<4; j++)
         {
         pt+=plaquettep(GC, r, i, j);
         }
     
      for(i=1; i<4; i++)
         {
         for(j=i+1; j<4; j++)
            {
            ps+=plaquettep(GC, r, i, j);
            }
         }
      }

   ps*=param->d_inv_vol;
   ps/=3.0;

   pt*=param->d_inv_vol;
   pt/=3.0;

   *plaqs=ps;
   *plaqt=pt;
   }


/* compute the mean Polyakov loop (the trace of) */
void polyakow(Gauge_Conf const *__restrict__ const GC, 
              Const const *__restrict__ const param, 
              double *__restrict__ repoly, 
              double *__restrict__ impoly) 
   {
   int t, r;
   double rp, ip;
   GAUGE_GROUP M;

   rp=0.0;
   ip=0.0;

   for(r=0; r<(param->d_space_vol); r++)
      {
      one(&M);
      for(t=0; t<(param->d_latot); t++)  
         {
         times_equal(&M, &(GC->lattice[param->d_latot*r+t][0]) );
         }

      rp+=retr(&M); 
      ip+=imtr(&M);  
      }
   rp*=param->d_inv_space_vol;
   ip*=param->d_inv_space_vol;

   *repoly=rp;
   *impoly=ip;
   }


/* compute the topological charge */
double topcharge(Gauge_Conf const *__restrict__ const GC, Const const *__restrict__ const param)
   {
   GAUGE_GROUP aux1, aux2, aux3;
   double ris, real1, real2, loc_charge; 
   #ifdef G2_GROUP
     const double chnorm=1.0/(128.0*PI*PI);
   #else
     const double chnorm=1.0/(256.0*PI*PI);
   #endif

   int r, i, dir[4][4], sign;

   dir[1][1] = 1;
   dir[1][2] = 1;
   dir[1][3] = 1;

   dir[2][1] = 2;
   dir[2][2] = 3;
   dir[2][3] = 0;

   dir[3][1] = 3;
   dir[3][2] = 2;
   dir[3][3] = 2;

   dir[0][1] = 0;
   dir[0][2] = 0;
   dir[0][3] = 3;

   ris=0.0;
   for(r=0; r<(param->d_volume); r++)
      {

      sign=1;
      loc_charge=0.0;

      for(i=1; i<4; i++)
         {
         quadrifoglio(GC, r, dir[1][i], dir[2][i], &aux1);
         quadrifoglio(GC, r, dir[3][i], dir[0][i], &aux2);

         times_dag2(&aux3, &aux2, &aux1); /* aux3=aux2*(aux1^{dag}) */
         real1=retr(&aux3)*NCOLOR; /* retr = real part of the trace / COLOR */

         times(&aux3, &aux2, &aux1); /* aux3=aux2*aux1 */
         real2=retr(&aux3)*NCOLOR;
        
         loc_charge+=(sign*(real1-real2));
         sign=-sign;
         }
      ris+=(loc_charge*chnorm); 
      }

   return ris;
   }

/* compute the topological charge density at point "r" */
double topchargedens(Gauge_Conf const *__restrict__ const GC, int r) 
   {
   GAUGE_GROUP aux1, aux2, aux3;
   double real1, real2, loc_charge; 
   #ifdef G2_GROUP
     const double chnorm=1.0/(128.0*PI*PI);
   #else
     const double chnorm=1.0/(256.0*PI*PI);
   #endif
   int dir[4][4], sign, i;

   dir[1][1] = 1;
   dir[1][2] = 1;
   dir[1][3] = 1;

   dir[2][1] = 2;
   dir[2][2] = 3;
   dir[2][3] = 0;

   dir[3][1] = 3;
   dir[3][2] = 2;
   dir[3][3] = 2;

   dir[0][1] = 0;
   dir[0][2] = 0;
   dir[0][3] = 3;

   sign=1;
   loc_charge=0.0;

   for(i=1; i<4; i++)
      {
      quadrifoglio(GC, r, dir[1][i], dir[2][i], &aux1);
      quadrifoglio(GC, r, dir[3][i], dir[0][i], &aux2);
 
         times_dag2(&aux3, &aux2, &aux1); /* aux3=aux2*(aux1^{dag}) */
         real1=retr(&aux3)*NCOLOR; /* retr = real part of the trace / NCOLOR */

         times(&aux3, &aux2, &aux1); /* aux3=aux2*aux1 */
         real2=retr(&aux3)*NCOLOR;
        
      loc_charge+=(sign*(real1-real2));
      sign=-sign;
      }
   return (loc_charge*chnorm); 
   }



/* compute Const::d_nummeas values of the topological charge after some cooling */
void topcharge_cooling(Gauge_Conf const *__restrict__ const GC, 
                       Const const *__restrict__ const param, 
                       double *__restrict__ charge, 
                       double *__restrict__ meanplaq) 
   {
   if(param->d_cooling>0)  /* if using cooling */
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter; 

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);  /* helperconf is a copy of the configuration */
  
     for(iter=0; iter<(param->d_nummeas); iter++)
        {
        cooling(&helperconf, param, param->d_cooling);

        ris=topcharge(&helperconf, param);
        charge[iter]=ris;

        plaquette(&helperconf, param, &plaqs, &plaqt);
        meanplaq[iter]=0.5*(plaqs+plaqt);
        }

     end_gauge_conf(&helperconf, param); 
     }
   else   /* no cooling */
     {
     double ris, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, param);
     plaquette(GC, param, &plaqs, &plaqt);

  
     for(iter=0; iter<(param->d_nummeas); iter++)
        {
        charge[iter]=ris;
        meanplaq[iter]=0.5*(plaqs+plaqt);
        }
     } 
   }


/* compute Const::d_nummeas values of the topological charge after some cooling2 */
void topcharge_cooling2(Gauge_Conf const *__restrict__ const GC, 
                        Const const *__restrict__ const param, 
                        double *__restrict__ charge, 
                        double *__restrict__ meanplaq) 
   {
   if(param->d_cooling>0)  /* if using cooling */
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter; 

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);  /* helperconf is a copy of the configuration */
  
     for(iter=0; iter<(param->d_nummeas); iter++)
        {
        cooling2(&helperconf, param, param->d_cooling);

        ris=topcharge(&helperconf, param);
        charge[iter]=ris;

        plaquette(&helperconf, param, &plaqs, &plaqt);
        meanplaq[iter]=0.5*(plaqs+plaqt);
        }

     end_gauge_conf(&helperconf, param); 
     }
   else   /* no cooling */
     {
     double ris, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, param);
     plaquette(GC, param, &plaqs, &plaqt);

  
     for(iter=0; iter<(param->d_nummeas); iter++)
        {
        charge[iter]=ris;
        meanplaq[iter]=0.5*(plaqs+plaqt);
        }
     } 
   }


#endif
