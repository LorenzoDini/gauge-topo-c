#ifndef SU2_UPD_C
#define SU2_UPD_C

#include"../Macro/macro.h"

#include<math.h>
#include<stdio.h>

#include"su2.h"
#include"su2_upd.h"
#include"../Const/const.h"
#include"../Rng/random.h"

/* random number generator for heabath (see Kennedy, Pendleton Phys. Lett. B 156, 393 (1985))
   given "k" return "a" in [-1,1] with P(a) = sqrt(1-a*a)*exp(k*a) */
void randheat_Su2(double k, double *__restrict__ out)
    {
    double r, r1, r2, c, r3;

    if(k>1.6847)
      {
      /* kennedy pendleton */
      while(1==1)
        {
        r=casuale();
        r1=casuale();

        r=-log(r);  
        r/=k;
        r1=-log(r1); 
        r1/=k;
 
        r2=casuale();
        c=cos(r2*PI2);
        c*=c;

        r*=c;  
        r1+=r;  

        r3=casuale();
        r=1-r1*0.5-r3*r3;
        if(r>0)
          {
          *out=1-r1;
          break;
          }
        }
      }
    else
      {
      /* creutz */
      r1=exp(-2.0*k);
      while(1==1)
        {
        do
        r=casuale();
        while(r<r1);
        r=1.0+log(r)/k;
 
        r2=casuale();
        if(r*r<r2*(2.0-r2))
          {
          *out=r;
          break;
          }
        }
      }
    }


/* heatbath */
void single_heatbath_Su2(Su2 *__restrict__ link, Su2 const *__restrict__ const staple, Const const *__restrict__ const param)
    {
    double p, p0;
    FILE *fp;  
    Su2 matrix1, matrix2;

    equal_Su2(&matrix1, staple);                    /* matrix1=staple */
    times_equal_real_Su2(&matrix1, param->d_beta);  /* matrix1*=d_beta */

    p=sqrtdet_Su2(&matrix1);
    if(p>MIN_VALUE)
      {
      times_equal_real_Su2(&matrix1, 1.0/p);  /* matrix1 *= 1.0/p */
      equal_dag_Su2(&matrix2, &matrix1);      /* matrix2 = matrix1^{dag} */

      randheat_Su2(p, &p0);
 
      rand_matrix_p0_Su2(p0, link);
      times_equal_Su2(link, &matrix2);        /* link*=matrix2 */
      }
    else
      {
      fp=fopen(param->err_file, "a");
      fprintf(fp, "Warning: in the heatbath in Su2_upd.cc p=%g < min_value\n",  p);
      fclose(fp);
      }
    }


/* overrelaxation */
void single_overrelaxation_Su2(Su2 *__restrict__ link, Su2 const *__restrict__ const staple, Const const *__restrict__ const param)
    {
    double p;
    Su2 matrix1, matrix2;
    FILE *fp; 

    equal_Su2(&matrix1, staple); /* matrix1=staple */

    p=sqrtdet_Su2(&matrix1);   
    if(p>MIN_VALUE)

      {
      times_equal_real_Su2(&matrix1, 1.0/p);     /* matrix1*=(1/p) */
      equal_dag_Su2(&matrix2, &matrix1);         /* matrix2=matrix1^{dag} */

      times_dag1_Su2(&matrix1, &matrix2, link);  /* matrix1=matrix2*link^{dag} */
      times_Su2(link, &matrix1, &matrix2);       /* link=matrix1*matrix2 */
      }
    else
      {
      fp=fopen(param->err_file, "a");
      fprintf(fp, "Warning: in the overrelaxation in Su2_upd.cc p=%g < min_value\n",  p);
      fclose(fp);
      }
    }


/* cooling */
void cool_Su2(Su2 *__restrict__ link, Su2 const *__restrict__ const staple)
    {
    Su2 matrix1;
    
    equal_Su2(&matrix1, staple);   /* matrix1=staple */
    unitarize_Su2(&matrix1);

    equal_dag_Su2(link, &matrix1);        /* link=matrix1^{dag} */
    }


#endif
