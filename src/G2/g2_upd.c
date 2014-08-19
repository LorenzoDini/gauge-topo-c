#ifndef G2_UPD_C
#define G2_UPD_C

#include<math.h>

#include"../Const/const.h"
#include"g2.h"
#include"g2_check.h"
#include"g2_un.h"
#include"g2_upd.h"
#include"../Macro/macro.h"
#include"../Rng/random.h"
#include"../Su2/su2_upd.h"

void single_heatbath_G2(G2 *__restrict__ link, G2 const *__restrict__ const staple, Const const *__restrict__ const param)
   {
   G2 aux, mult;
   double w[4], w_mod, k[4], temp[4], temp_mod, norm;
   int count, mode;
   FILE *fp;

   for(count=0;count<6;count++)
      {
      mode=count%3;

      equal_G2(&aux, staple);
      times_equal_G2(&aux, link);

      if(mode==0)
        {
        /* update on the first SU(2) subgroup, see pag. 19 of  Greensite, Langfeld, Olejnik, Reinhardt, Tok   Phys Rev D 75, p.034501 (2007) */

        /* action = (param->d_beta)*(w*k) */
        w[0]=aux.comp[3][3]+aux.comp[4][4]+aux.comp[5][5]+aux.comp[6][6];
        w[1]=aux.comp[3][6]-aux.comp[6][3]+aux.comp[4][5]-aux.comp[5][4];
        w[2]=aux.comp[5][3]-aux.comp[3][5]+aux.comp[4][6]-aux.comp[6][4];
        w[3]=aux.comp[3][4]-aux.comp[4][3]+aux.comp[5][6]-aux.comp[6][5];

        w_mod=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
        norm=w_mod*(param->d_beta);

        if(norm>MIN_VALUE)
          {
          /* generate the component parallel to a w (i.e. t[0]) with heatbath, the other are random */
          randheat(norm, &temp_mod);
          temp[0]=temp_mod;
          temp[1]=1.0-2.0*casuale();
          temp[2]=1.0-2.0*casuale();
          temp[3]=1.0-2.0*casuale();

          /* normalize temp */
          temp_mod=sqrt((1.0-temp[0]*temp[0])/(temp[1]*temp[1]+temp[2]*temp[2]+temp[3]*temp[3]));
          temp[1]*=temp_mod;
          temp[2]*=temp_mod;
          temp[3]*=temp_mod;

          /* normalize w */
          w[0]/=w_mod;
          w[1]/=w_mod;
          w[2]/=w_mod;
          w[3]/=w_mod;

          /* rotate temp in the w direction, result is k */
          if(1.0+w[0] > MIN_VALUE || 1.0+w[0]< -MIN_VALUE) 
            {
            norm=1.0+w[0];

            k[0]=w[0]*temp[0]                -w[1]*temp[1]                -w[2]*temp[2]                -w[3]*temp[3];
            k[1]=w[1]*temp[0]+(1.0-w[1]*w[1]/norm)*temp[1]      -w[2]*w[1]/norm*temp[2]      -w[3]*w[1]/norm*temp[3];
            k[2]=w[2]*temp[0]      -w[1]*w[2]/norm*temp[1]+(1.0-w[2]*w[2]/norm)*temp[2]      -w[3]*w[2]/norm*temp[3];
            k[3]=w[3]*temp[0]      -w[1]*w[3]/norm*temp[1]      -w[2]*w[3]/norm*temp[2]+(1.0-w[3]*w[3]/norm)*temp[3];
            }
          else
            {
            norm=1.0-w[0];
  
            k[0]=w[0]*temp[0]                +w[1]*temp[1]                +w[2]*temp[2]                +w[3]*temp[3];
            k[1]=w[1]*temp[0]+(1.0-w[1]*w[1]/norm)*temp[1]      -w[2]*w[1]/norm*temp[2]      -w[3]*w[1]/norm*temp[3];
            k[2]=w[2]*temp[0]      -w[1]*w[2]/norm*temp[1]+(1.0-w[2]*w[2]/norm)*temp[2]      -w[3]*w[2]/norm*temp[3];
            k[3]=w[3]*temp[0]      -w[1]*w[3]/norm*temp[1]      -w[2]*w[3]/norm*temp[2]+(1.0-w[3]*w[3]/norm)*temp[3];
            }

          zero_G2(&mult);
          mult.comp[0][0]=1.0;
          mult.comp[1][1]=1.0;
          mult.comp[2][2]=1.0;

          mult.comp[3][3]=k[0];
          mult.comp[3][4]=-k[3];
          mult.comp[3][5]=k[2];
          mult.comp[3][6]=-k[1];
       
          mult.comp[4][3]=k[3];
          mult.comp[4][4]=k[0];
          mult.comp[4][5]=-k[1];
          mult.comp[4][6]=-k[2];

          mult.comp[5][3]=-k[2];
          mult.comp[5][4]=k[1];
          mult.comp[5][5]=k[0];
          mult.comp[5][6]=-k[3];
  
          mult.comp[6][3]=k[1];
          mult.comp[6][4]=k[2];
          mult.comp[6][5]=k[3];
          mult.comp[6][6]=k[0];

          /* times_equal_G2(link, &mult); */
          times_equal_4_G2(link, &mult, 3, 4, 5, 6);
          }
       else
          {
          fp=fopen(param->err_file, "a");
          fprintf(fp, "Warning:  in heatbath in g2_upd.cc norm = %g < min_value\n",  norm);
          fclose(fp);
          }
       }

      if(mode==1)
        {
        /* update on the second SU(2) subgroup, see pag. 19 of  Greensite, Langfeld, Olejnik, Reinhardt, Tok   Phys Rev D 75, p.034501 (2007) */

        /* action = (param->d_beta)*(w*k) */
        w[0]=aux.comp[1][1]+aux.comp[2][2]+aux.comp[5][5]+aux.comp[6][6];
        w[1]=aux.comp[6][1]-aux.comp[1][6]+aux.comp[5][2]-aux.comp[2][5];
        w[2]=aux.comp[1][5]-aux.comp[5][1]+aux.comp[6][2]-aux.comp[2][6];
        w[3]=aux.comp[1][2]-aux.comp[2][1]+aux.comp[5][6]-aux.comp[6][5];

        w_mod=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
        norm=w_mod*(param->d_beta);

        if(norm>MIN_VALUE)
          {
          /* generate the component parallel to a w (i.e. t[0]) with heatbath, the other are random */
          randheat(norm, &temp_mod);
          temp[0]=temp_mod;
          temp[1]=1.0-2.0*casuale();
          temp[2]=1.0-2.0*casuale();
          temp[3]=1.0-2.0*casuale();

          /* normalize temp */
          temp_mod=sqrt((1.0-temp[0]*temp[0])/(temp[1]*temp[1]+temp[2]*temp[2]+temp[3]*temp[3]));
          temp[1]*=temp_mod;
          temp[2]*=temp_mod;
          temp[3]*=temp_mod;

          /* normalize w */
          w[0]/=w_mod;
          w[1]/=w_mod;
          w[2]/=w_mod;
          w[3]/=w_mod;

          /* rotate temp in the w direction, result is k */
          if(1.0+w[0] > MIN_VALUE || 1.0+w[0]< -MIN_VALUE) 
            {
            norm=1.0+w[0];

            k[0]=w[0]*temp[0]                -w[1]*temp[1]                -w[2]*temp[2]                -w[3]*temp[3];
            k[1]=w[1]*temp[0]+(1.0-w[1]*w[1]/norm)*temp[1]      -w[2]*w[1]/norm*temp[2]      -w[3]*w[1]/norm*temp[3];
            k[2]=w[2]*temp[0]      -w[1]*w[2]/norm*temp[1]+(1.0-w[2]*w[2]/norm)*temp[2]      -w[3]*w[2]/norm*temp[3];
            k[3]=w[3]*temp[0]      -w[1]*w[3]/norm*temp[1]      -w[2]*w[3]/norm*temp[2]+(1.0-w[3]*w[3]/norm)*temp[3];
            }
          else
            {
            norm=1.0-w[0];

            k[0]=w[0]*temp[0]                +w[1]*temp[1]                +w[2]*temp[2]                +w[3]*temp[3];
            k[1]=w[1]*temp[0]+(1.0-w[1]*w[1]/norm)*temp[1]      -w[2]*w[1]/norm*temp[2]      -w[3]*w[1]/norm*temp[3];
            k[2]=w[2]*temp[0]      -w[1]*w[2]/norm*temp[1]+(1.0-w[2]*w[2]/norm)*temp[2]      -w[3]*w[2]/norm*temp[3];
            k[3]=w[3]*temp[0]      -w[1]*w[3]/norm*temp[1]      -w[2]*w[3]/norm*temp[2]+(1.0-w[3]*w[3]/norm)*temp[3];
            }
        
          zero_G2(&mult);
          mult.comp[0][0]=1.0;
          mult.comp[3][3]=1.0;
          mult.comp[4][4]=1.0;

          mult.comp[1][1]=k[0];
          mult.comp[1][2]=-k[3];
          mult.comp[1][5]=-k[2];
          mult.comp[1][6]=k[1];
       
          mult.comp[2][1]=k[3];
          mult.comp[2][2]=k[0];
          mult.comp[2][5]=k[1];
          mult.comp[2][6]=k[2];

          mult.comp[5][1]=k[2];
          mult.comp[5][2]=-k[1];
          mult.comp[5][5]=k[0];
          mult.comp[5][6]=-k[3];
  
          mult.comp[6][1]=-k[1];
          mult.comp[6][2]=-k[2];
          mult.comp[6][5]=k[3];
          mult.comp[6][6]=k[0];

          /* times_equal_G2(link, &mult); */
          times_equal_4_G2(link, &mult, 1, 2, 5, 6);
          }
        else
          {
          fp=fopen(param->err_file, "a");
          fprintf(fp, "Warning:  in heatbath in g2_upd.cc norm = %g < min_value\n",  norm);
          fclose(fp);
          }
        }

      if(mode==2)
        {
        /* update on the third SU(2) subgroup, see pag. 19 of  Greensite, Langfeld, Olejnik, Reinhardt, Tok   Phys Rev D 75, p.034501 (2007) */

        /* action = (param->d_beta)*(w*k) */
        w[0]=aux.comp[1][1]+aux.comp[2][2]+aux.comp[3][3]+aux.comp[4][4];
        w[1]=aux.comp[4][1]-aux.comp[1][4]+aux.comp[2][3]-aux.comp[3][2];
        w[2]=aux.comp[1][3]-aux.comp[3][1]+aux.comp[2][4]-aux.comp[4][2];
        w[3]=aux.comp[3][4]-aux.comp[4][3]+aux.comp[2][1]-aux.comp[1][2];

        w_mod=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
        norm=w_mod*(param->d_beta);

        if(norm>MIN_VALUE)
          {
          /* generate the component parallel to a w (i.e. t[0]) with heatbath, the other are random */
          randheat(norm, &temp_mod);
          temp[0]=temp_mod;
          temp[1]=1.0-2.0*casuale();
          temp[2]=1.0-2.0*casuale();
          temp[3]=1.0-2.0*casuale();

          /* normalize temp */
          temp_mod=sqrt((1.0-temp[0]*temp[0])/(temp[1]*temp[1]+temp[2]*temp[2]+temp[3]*temp[3]));
          temp[1]*=temp_mod;
          temp[2]*=temp_mod;
          temp[3]*=temp_mod;

          /* normalize w */
          w[0]/=w_mod;
          w[1]/=w_mod;
          w[2]/=w_mod;
          w[3]/=w_mod;

          /* rotate temp in the w direction, result is k */
          if(1.0+w[0] > MIN_VALUE || 1.0+w[0]< -MIN_VALUE) 
            {
            norm=1.0+w[0];

            k[0]=w[0]*temp[0]                -w[1]*temp[1]                -w[2]*temp[2]                -w[3]*temp[3];
            k[1]=w[1]*temp[0]+(1.0-w[1]*w[1]/norm)*temp[1]      -w[2]*w[1]/norm*temp[2]      -w[3]*w[1]/norm*temp[3];
            k[2]=w[2]*temp[0]      -w[1]*w[2]/norm*temp[1]+(1.0-w[2]*w[2]/norm)*temp[2]      -w[3]*w[2]/norm*temp[3];
            k[3]=w[3]*temp[0]      -w[1]*w[3]/norm*temp[1]      -w[2]*w[3]/norm*temp[2]+(1.0-w[3]*w[3]/norm)*temp[3];
            }
          else
            {
            norm=1.0-w[0];

            k[0]=w[0]*temp[0]                +w[1]*temp[1]                +w[2]*temp[2]                +w[3]*temp[3];
            k[1]=w[1]*temp[0]+(1.0-w[1]*w[1]/norm)*temp[1]      -w[2]*w[1]/norm*temp[2]      -w[3]*w[1]/norm*temp[3];
            k[2]=w[2]*temp[0]      -w[1]*w[2]/norm*temp[1]+(1.0-w[2]*w[2]/norm)*temp[2]      -w[3]*w[2]/norm*temp[3];
            k[3]=w[3]*temp[0]      -w[1]*w[3]/norm*temp[1]      -w[2]*w[3]/norm*temp[2]+(1.0-w[3]*w[3]/norm)*temp[3];
            }
       
          zero_G2(&mult); 
          mult.comp[0][0]=1.0;
          mult.comp[5][5]=1.0;
          mult.comp[6][6]=1.0;

          mult.comp[1][1]=k[0];
          mult.comp[1][2]=k[3];
          mult.comp[1][3]=-k[2];
          mult.comp[1][4]=k[1];
       
          mult.comp[2][1]=-k[3];
          mult.comp[2][2]=k[0];
          mult.comp[2][3]=-k[1];
          mult.comp[2][4]=-k[2];

          mult.comp[3][1]=k[2];
          mult.comp[3][2]=k[1];
          mult.comp[3][3]=k[0];
          mult.comp[3][4]=-k[3];
  
          mult.comp[4][1]=-k[1];
          mult.comp[4][2]=k[2];
          mult.comp[4][3]=k[3];
          mult.comp[4][4]=k[0];

          /* times_equal_G2(link, &mult); */
          times_equal_4_G2(link, &mult, 1, 2, 3, 4);
          }
        else
          {
          fp=fopen(param->err_file, "a");
          fprintf(fp, "Warning:  in heatbath in g2_upd.cc norm = %g < min_value\n",  norm);
          fclose(fp);
          }
        }
      }
   }





void single_overrelaxation_G2(G2 *__restrict__ link, G2 const *__restrict__ const staple, Const const *__restrict__ const param)
   {
   G2 M, aux;
   double w0, w1;
   double phi;

   equal_G2(&aux, staple);
   times_equal_G2(&aux, link);

   /* D1 */
   w0=aux.comp[3][3]+aux.comp[4][4]+aux.comp[5][5]+aux.comp[6][6];
   w1=aux.comp[3][6]-aux.comp[6][3]+aux.comp[4][5]-aux.comp[5][4];

   /* phi such that w0=sqrt(w0^2+w1^2)cos(phi)    w1=sqrt(w0^2+w1^2)sin(phi) */
   phi=atan(w1/w0);
   if(w0<0) phi+=PI;
   
   D1(&M, 2.0*phi);
   /*
   times_equal_G2(link, &M);
   times_equal_G2(&aux, &M);
   */
   times_equal_4_G2(link, &M, 3, 4, 5, 6);
   times_equal_4_G2(&aux, &M, 3, 4, 5, 6);

   /* D2 */
   w0=aux.comp[3][3]+aux.comp[4][4]+aux.comp[5][5]+aux.comp[6][6];
   w1=aux.comp[5][3]-aux.comp[3][5]+aux.comp[4][6]-aux.comp[6][4];

   /* phi such that w0=sqrt(w0^2+w1^2)cos(phi)    w1=sqrt(w0^2+w1^2)sin(phi) */
   phi=atan(w1/w0);
   if(w0<0) phi+=PI;
   
   D2(&M, 2.0*phi);
   /*
   times_equal_G2(link, &M);
   times_equal_G2(&aux, &M);
   */
   times_equal_4_G2(link, &M, 3, 4, 5, 6);
   times_equal_4_G2(&aux, &M, 3, 4, 5, 6);

   /* D3 */
   w0=aux.comp[3][3]+aux.comp[4][4]+aux.comp[5][5]+aux.comp[6][6];
   w1=aux.comp[3][4]-aux.comp[4][3]+aux.comp[5][6]-aux.comp[6][5];

   /* phi such that w0=sqrt(w0^2+w1^2)cos(phi)    w1=sqrt(w0^2+w1^2)sin(phi) */
   phi=atan(w1/w0);
   if(w0<0) phi+=PI;
   
   D3(&M, 2.0*phi);
   /*
   times_equal_G2(link, &M);
   times_equal_G2(&aux, &M);
   */
   times_equal_4_G2(link, &M, 3, 4, 5, 6);
   times_equal_4_G2(&aux, &M, 3, 4, 5, 6);

   /* D4 */
   w0=aux.comp[1][1]+aux.comp[2][2]+aux.comp[5][5]+aux.comp[6][6];
   w1=aux.comp[6][1]-aux.comp[1][6]+aux.comp[5][2]-aux.comp[2][5];

   /* phi such that w0=sqrt(w0^2+w1^2)cos(phi)    w1=sqrt(w0^2+w1^2)sin(phi) */
   phi=atan(w1/w0);
   if(w0<0) phi+=PI;
   
   D4(&M, 2.0*phi);
   /*
   times_equal_G2(link, &M);
   times_equal_G2(&aux, &M);
   */
   times_equal_4_G2(link, &M, 1, 2, 5, 6);
   times_equal_4_G2(&aux, &M, 1, 2, 5, 6);

   /* D5 */
   w0=aux.comp[1][1]+aux.comp[2][2]+aux.comp[5][5]+aux.comp[6][6];
   w1=aux.comp[1][5]-aux.comp[5][1]+aux.comp[6][2]-aux.comp[2][6];

   /* phi such that w0=sqrt(w0^2+w1^2)cos(phi)    w1=sqrt(w0^2+w1^2)sin(phi) */
   phi=atan(w1/w0);
   if(w0<0) phi+=PI;
   
   D5(&M, 2.0*phi);
   /*
   times_equal_G2(link, &M);
   times_equal_G2(&aux, &M);
   */
   times_equal_4_G2(link, &M, 1, 2, 5, 6);
   times_equal_4_G2(&aux, &M, 1, 2, 5, 6);


   /* D6 */
   w0=aux.comp[1][1]+aux.comp[2][2]+aux.comp[3][3]+aux.comp[4][4];
   w1=aux.comp[2][3]-aux.comp[3][2]+aux.comp[4][1]-aux.comp[1][4];

   /* phi such that w0=sqrt(w0^2+w1^2)cos(phi)    w1=sqrt(w0^2+w1^2)sin(phi) */
   phi=atan(w1/w0);
   if(w0<0) phi+=PI;
   
   D6(&M, 2.0*phi);
   /*
   times_equal_G2(link, &M);
   times_equal_G2(&aux, &M);
   */
   times_equal_4_G2(link, &M, 1, 2, 3, 4);
   times_equal_4_G2(&aux, &M, 1, 2, 3, 4);

   /* D7 */
   w0=aux.comp[1][1]+aux.comp[2][2]+aux.comp[3][3]+aux.comp[4][4];
   w1=aux.comp[1][3]-aux.comp[3][1]+aux.comp[2][4]-aux.comp[4][2];

   /* phi such that w0=sqrt(w0^2+w1^2)cos(phi)    w1=sqrt(w0^2+w1^2)sin(phi) */
   phi=atan(w1/w0);
   if(w0<0) phi+=PI;
   
   D7(&M, 2.0*phi);
   /*
   times_equal_G2(link, &M);
   */
   times_equal_4_G2(link, &M, 1, 2, 3, 4);
   }


void cool_G2(G2 *__restrict__ link, G2 const *__restrict__ const staple)
   {
   G2 helper, helper2;
   double energy1, energy2;
   double beta_aux=60.0;

   equal_G2(&helper, link); 
   times_G2(&helper2, &helper, staple); /*helper2=link*staple */
   energy1=-retr_G2(&helper2);

   heatbath_aux_G2(&helper, staple, beta_aux);

   times_G2(&helper2, &helper, staple); /*helper2=helper*staple */
   energy2=-retr_G2(&helper2);
   if(energy2<energy1)
     {
     equal_G2(link, &helper);
     }
   }

#endif
