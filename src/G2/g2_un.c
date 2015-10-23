#ifndef G2_UN_C
#define G2_UN_C

#include"../Macro/macro.h"

#include<math.h>

#include"g2.h"
#include"g2_check.h"
#include"g2_un.h"
#include"../Rng/random.h"
#include"../Su2/su2_upd.h"

void heatbath_aux_G2(G2 * link, G2 const * const staple, double beta_aux)
   {
   G2 aux, mult;
   double w[4], w_mod, k[4], temp[4], temp_mod, norm;
   int count, mode;

   for(count=0;count<6;count++)
      {
      mode=count%3;

      equal_G2(&aux, staple);
      times_equal_G2(&aux, link);

      if(mode==0)
        {
        /* action = beta_aux*(w*k) */
        w[0]=aux.comp[3][3]+aux.comp[4][4]+aux.comp[5][5]+aux.comp[6][6];
        w[1]=aux.comp[3][6]-aux.comp[6][3]+aux.comp[4][5]-aux.comp[5][4];
        w[2]=aux.comp[5][3]-aux.comp[3][5]+aux.comp[4][6]-aux.comp[6][4];
        w[3]=aux.comp[3][4]-aux.comp[4][3]+aux.comp[5][6]-aux.comp[6][5];

        w_mod=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
        norm=w_mod*beta_aux;

        /* generate the component parallel to a w (i.e. t[0]) with heatbath, the other are random */
        randheat_Su2(norm, &temp_mod);
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

        times_equal_G2(link, &mult);
        }

      if(mode==1)
        {
        /* action = beta_aux*(w*k) */
        w[0]=aux.comp[1][1]+aux.comp[2][2]+aux.comp[5][5]+aux.comp[6][6];
        w[1]=aux.comp[6][1]-aux.comp[1][6]+aux.comp[5][2]-aux.comp[2][5];
        w[2]=aux.comp[1][5]-aux.comp[5][1]+aux.comp[6][2]-aux.comp[2][6];
        w[3]=aux.comp[1][2]-aux.comp[2][1]+aux.comp[5][6]-aux.comp[6][5];

        w_mod=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
        norm=w_mod*beta_aux;

        /* generate the component parallel to a w (i.e. t[0]) with heatbath, the other are random */
        randheat_Su2(norm, &temp_mod);
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

        times_equal_G2(link, &mult);
        }

      if(mode==2)
        {
        /* action = beta_aux*(w*k) */
        w[0]=aux.comp[1][1]+aux.comp[2][2]+aux.comp[3][3]+aux.comp[4][4];
        w[1]=aux.comp[4][1]-aux.comp[1][4]+aux.comp[2][3]-aux.comp[3][2];
        w[2]=aux.comp[1][3]-aux.comp[3][1]+aux.comp[2][4]-aux.comp[4][2];
        w[3]=aux.comp[3][4]-aux.comp[4][3]+aux.comp[2][1]-aux.comp[1][2];

        w_mod=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
        norm=w_mod*beta_aux;

        /* generate the component parallel to a w (i.e. t[0]) with heatbath, the other are random */
        randheat_Su2(norm, &temp_mod);
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

        times_equal_G2(link, &mult);
        }
      }
   }


void unitarize_G2(G2 * A)
   {
   double nn, cc;
   double beta_aux=1.0e+12;
   int n;
   G2 force, guess, random, helper;
   
   one_G2(&guess);
   while(!check_G2(A, &nn, &cc))
        {
        for(n=0; n<20; n++)
           {
           equal_G2(&force, A);

           rand_matrix_G2(&random);

           times_dag2_G2(&helper, &force, &random);
           times_G2(&force, &random, &helper);    /* force->random*force*(random^{dag}) */

           times_dag2_G2(&helper, &guess, &random);
           times_G2(&guess, &random, &helper);    /* guess->random*guess*(random^{dag}) */

           heatbath_aux_G2(&guess, &force, beta_aux);

           times_G2(&helper, &guess, &random);
           times_dag1_G2(&guess, &random, &helper);    /* guess->random^{dag}*guess*random */
           }
        equal_dag_G2(&force, &guess);
        
        equal_G2(A, &force);
        }
   }


#endif
