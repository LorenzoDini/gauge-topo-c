#ifndef SON_UPD_C
#define SON_UPD_C

#include<math.h>
#include<stdio.h>

#include"../Const/const.h"
#include"../Macro/macro.h"
#include"../Rng/random.h"
#include"son.h"

/* link and staple are here 2-dim real vectors, this is the 
   function for heatbath in U(1) gauge theory 
   see Moriarty Phys. Rev. D 25, pag. 2185 (1982) */
void randheat_U1(double *link, double const * const staple)
    {
    double alpha, theta, theta_staple, q, q_max, r1, r2, x;
    
    alpha=sqrt(staple[0]*staple[0]+staple[1]*staple[1]);
    q_max=exp(0.21051366235301868432776943515*alpha);    /* see the following computation */ 
    theta_staple=atan2(staple[1], staple[0] );

    do
      {
      r1=casuale();
      r2=casuale();
      x = -1.0 + log(1.0 + (exp(2.0*alpha) - 1.0)*r1)/alpha;
      q=exp(alpha*(cos(PI*0.5*(1.0-x)) - x));
      } 
    while(q<=q_max*r2);

    theta = (1.0 - x)*HALF_PI;
    r1=casuale();
    if(r1>0.5) 
      {
      theta=-theta;
      }
    theta+=theta_staple;
     
    link[0]=cos(theta);
    link[1]=sin(theta);
    }

/*
WorkingPrecision->1000;
Q[x_]:=Exp[Cos[Pi/2*(1-x)]-x]
FindRoot[D[Q[x],x]==0, {x,1}, WorkingPrecision->100]
Out[52]= {x -> 0.5606641805798867176366776048997096707812104519411362714885751166519976969907076829844764496162691569}
N[Log[Q[x]/.%52], 100]
Out[53]= 0.210513662353018684327769435155832317434879346989632455087165428289411464536813734686515674410131331
*/


/* Pseudo-heatbath a la Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) */ 
void single_heatbath_SoN(SoN *__restrict__ link, 
                         SoN const *__restrict__ const staple, 
                         Const const *__restrict__ const param)
    {
    SoN aux;
    double link2[2], staple2[2], temp[2], norm;
    int i, j, k;
    FILE *fp; 

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          equal_SoN(&aux, staple);     /* aux=staple */
          times_equal_SoN(&aux, link); /*aux=staple*link */  

          /* update on the SO(2) subgroup of raw, columns = i,j */
 
          staple2[0]=param->d_beta*(aux.comp[i][i]+aux.comp[j][j])/NCOLOR;
          staple2[1]=param->d_beta*(aux.comp[j][i]-aux.comp[i][j])/NCOLOR;

          /* action = Tr(staple2*link2) + const */

          norm=sqrt(staple2[0]*staple2[0]+staple2[1]*staple2[1]);

          if(norm>MIN_VALUE)
            {
            randheat_U1(link2, staple2);

            /* link*=final */
            for(k=0; k<NCOLOR; k++)
               {
               temp[0]=link->comp[k][i]*link2[0] - link->comp[k][j]*link2[1];
               temp[1]=link->comp[k][i]*link2[1] + link->comp[k][j]*link2[0];
               link->comp[k][i]=temp[0];
               link->comp[k][j]=temp[1];
               }
            }
          else
            {
            fp=fopen(param->err_file, "a");
            fprintf(fp, "Warning: nell'heatbath in son_upd.c norm=%g < min_value\n", norm); 
            fclose(fp);
            }
          }
       }
    }


/* Pseudo-overrelaxation a la Cabibbo-Marinari (Phys. Lett. B 119, p.387 (1982)) */ 
void single_overrelaxation_SoN(SoN *__restrict__ link, 
                               SoN const *__restrict__ const staple, 
                               Const const *__restrict__ const param)
    {
    SoN aux;
    double link2[2], staple2[2], temp[2], norm, theta_staple;
    int i, j, k;
    FILE *fp; 

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          equal_SoN(&aux, staple);     /* aux=staple */
          times_equal_SoN(&aux, link); /*aux=staple*link */  

          /* update on the SO(2) subgroup of raw, columns = i,j */
 
          staple2[0]=param->d_beta*(aux.comp[i][i]+aux.comp[j][j])/NCOLOR;
          staple2[1]=param->d_beta*(aux.comp[j][i]-aux.comp[i][j])/NCOLOR;

          /* action = Tr(staple2*link2) + const */

          norm=sqrt(staple2[0]*staple2[0]+staple2[1]*staple2[1]);

          if(norm>MIN_VALUE)
            {
            theta_staple=atan2(staple2[1], staple2[0]);

            link2[0]=cos(2.0*theta_staple);
            link2[1]=sin(2.0*theta_staple);

            /* link*=final */
            for(k=0; k<NCOLOR; k++)
               {
               temp[0]=link->comp[k][i]*link2[0] - link->comp[k][j]*link2[1];
               temp[1]=link->comp[k][i]*link2[1] + link->comp[k][j]*link2[0];
               link->comp[k][i]=temp[0];
               link->comp[k][j]=temp[1];
               }
            }
          else
            {
            fp=fopen(param->err_file, "a");
            fprintf(fp, "Warning: nell'overrelaxation in son_upd.c norm=%g < min_value\n", norm); 
            fclose(fp);
            }
          }
       }
    }


/* cooling */
void cool_SoN(SoN *__restrict__ link, 
              SoN const *__restrict__ const staple)
  {
    SoN aux;
    double link2[2], staple2[2], temp[2], norm;
    int i, j, k;

    for(i=0; i<NCOLOR-1; i++)
       {
       for(j=i+1; j<NCOLOR; j++)
          {
          equal_SoN(&aux, staple);     /* aux=staple */
          times_equal_SoN(&aux, link); /*aux=staple*link */  

          /* update on the SO(2) subgroup of raw, columns = i,j */
 
          staple2[0]=(aux.comp[i][i]+aux.comp[j][j])/NCOLOR;
          staple2[1]=(aux.comp[j][i]-aux.comp[i][j])/NCOLOR;

          /* action = Tr(staple2*link2) + const */

          norm=sqrt(staple2[0]*staple2[0]+staple2[1]*staple2[1]);

          if(norm>MIN_VALUE)
            {

            link2[0]=staple2[0]/norm;
            link2[1]=staple2[1]/norm;

            /* link*=final */
            for(k=0; k<NCOLOR; k++)
               {
               temp[0]=link->comp[k][i]*link2[0] - link->comp[k][j]*link2[1];
               temp[1]=link->comp[k][i]*link2[1] + link->comp[k][j]*link2[0];
               link->comp[k][i]=temp[0];
               link->comp[k][j]=temp[1];
               }
            }
          }
       }
    }

#endif
