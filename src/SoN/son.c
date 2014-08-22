#ifndef SON_C
#define SON_C

#include<math.h>
#include<stdio.h>

#include"../Macro/macro.h"
#include"../Rng/random.h"
#include"son.h"

/* A=1 */
void one_SoN(SoN *__restrict__ A)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        A->comp[i][j]=0.0;
        A->comp[j][i]=0.0;
        }
     A->comp[i][i]=1.0;
     }
  }

/* A=0 */
void zero_SoN(SoN *__restrict__ A)       
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  }

/* A=B */
void equal_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=B->comp[i][j];
        }
     }
  }

/* A=B^{dag} */
void equal_dag_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=(B->comp[j][i]);
        }
     }
  }


/* A+=B */
void plus_equal_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B) 
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]+=B->comp[i][j];
        }
     }
  }


/* A+=B^{dag} */
void plus_equal_dag_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]+=(B->comp[j][i]);
        }
     }
  }

/* A-=B */
void minus_equal_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]-=B->comp[i][j];
        }
     }
  }

/* A-=B^{dag} */
void minus_equal_dag_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]-=(B->comp[j][i]);
        }
     }
  }

/* A=b*B+c*C */
void lin_comb_SoN(SoN *__restrict__ A, 
                  double b, 
                  SoN const *__restrict__ const B, 
                  double c, 
                  SoN const *__restrict__ const C)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=b*(B->comp[i][j])+c*(C->comp[i][j]);
        }
     }
  }

/* A=b*B^{dag}+c*C */
void lin_comb_dag1_SoN(SoN *__restrict__ A, 
                       double b, 
                       SoN const *__restrict__ const B, 
                       double c, 
                       SoN const *__restrict__ const C)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=b*(B->comp[j][i])+c*(C->comp[i][j]);
        }
     }
  }

/* A=b*B+c*C^{dag} */
void lin_comb_dag2_SoN(SoN *__restrict__ A, 
                       double b, 
                       SoN const *__restrict__ const B, 
                       double c, 
                       SoN const *__restrict__ const C)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=b*(B->comp[i][j])+c*(C->comp[j][i]);
        }
     }
  }

/* A=b*B^{dag}+c*C^{dag} */
void lin_comb_dag12_SoN(SoN *__restrict__ A, 
                        double b, 
                        SoN const *__restrict__ const B, 
                        double c, 
                        SoN const *__restrict__ const C)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=b*(B->comp[j][i])+c*(C->comp[j][i]);
        }
     }
  }

/* A*=r */
void times_equal_real_SoN(SoN *__restrict__ A, double r)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]*=r;
        }
     }
  }

/* A*=B */
void times_equal_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B)
   {
   int i, j, k;
   double sum, aux[NCOLOR];

   for(i=0; i<NCOLOR; i++)
      {
      for(j=0; j<NCOLOR; j++)
         {
         aux[j]=A->comp[i][j];
         }
     
      for(j=0; j<NCOLOR; j++)
         {
         sum=0.0;
         for(k=0; k<NCOLOR; k++)
            {
            sum+=aux[k]*(B->comp[k][j]);
            }
         A->comp[i][j]=sum;
         }
      }
   }
 
/* A*=B^{dag} */
void times_equal_dag_SoN(SoN *__restrict__ A, SoN const *__restrict__ const B)
   {
   int i, j, k;
   double sum, aux[NCOLOR];

   for(i=0; i<NCOLOR; i++)
      {
      for(j=0; j<NCOLOR; j++)
         {
         aux[j]=A->comp[i][j];
         }
     
      for(j=0; j<NCOLOR; j++)
         {
         sum=0.0;
         for(k=0; k<NCOLOR; k++)
            {
            sum+=aux[k]*(B->comp[j][k]);
            }
         A->comp[i][j]=sum;
         }
      }
   }

/* A=B*C */
void times_SoN(SoN *__restrict__ A, 
               SoN const *__restrict__ const B, 
               SoN const *__restrict__ const C)
  {
  int i, j, k;
  double sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[i][k])*(C->comp[k][j]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B^{dag}*C */
void times_dag1_SoN(SoN *__restrict__ A, 
                    SoN const *__restrict__ const B, 
                    SoN const *__restrict__ const C) 
  {
  int i, j, k;
  double sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[k][i])*(C->comp[k][j]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B*C^{dag} */
void times_dag2_SoN(SoN *__restrict__ A, 
                    SoN const *__restrict__ const B, 
                    SoN const *__restrict__ const C) 
  {
  int i, j, k;
  double sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[i][k])*(C->comp[j][k]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B^{dag}*C^{dag} */
void times_dag12_SoN(SoN *__restrict__ A, 
                     SoN const *__restrict__ const B, 
                     SoN const *__restrict__ const C)
  {
  int i, j, k;
  double sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[k][i])*(C->comp[j][k]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* SO(N) random matrix
   generated a la Cabibbo Marinari with N(N-1)/2 SO(2) random matrices */
void rand_matrix_SoN(SoN *__restrict__ A)
  {
  int i, j, k;
  double p0, p1, p, temp[2];

  one_SoN(A);

  for(i=0; i<NCOLOR-1; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        /* SO(2) random components */
        p=2.0;
        while(p>1.0)
             { 
             p0=1.0-2.0*casuale();
             p1=1.0-2.0*casuale();
             p=sqrt(p0*p0+p1*p1);
             }

        p0/=p;
        p1/=p;

        for(k=0; k<NCOLOR; k++)
           {
           temp[0]=A->comp[k][i]*p0 - A->comp[k][j]*p1;
           temp[1]=A->comp[k][i]*p1 + A->comp[k][j]*p0;
           A->comp[k][i]=temp[0];
           A->comp[k][j]=temp[1];
           }
        }
     }
  }

/* l2 norm of the matrix */
double norm_SoN(SoN const *__restrict__ const A)
  {
  int i, j;
  double aux, ris;

  ris=0.0;
  for(i=0; i<NCOLOR; i++)
     { 
     for(j=0; j<NCOLOR; j++)
        {
        aux=A->comp[i][j];
        ris+=aux*aux;
        }
     }
  return sqrt(ris);
  }

/* real part of the trace /N */
double retr_SoN(SoN const *__restrict__ const A)
  {
  int i;
  double ris, tr;

  tr=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[i][i];
     }
  ris=tr/(double)NCOLOR;
  return ris;
  }

/* imaginary part of the trace /N */
double imtr_SoN(SoN const *__restrict__ const A)
  {
  return 0.0;
  }

/* LU decomposition with partial pivoting
   from Numerical Recipes in C, pag 46 */
void LU_SoN(SoN const *__restrict__ const A, SoN *__restrict__ ris, int *sign) 
  {
  int i, imax, j, k;
  double  big, temp, sum, dum, vv[NCOLOR];

  imax=0;
  equal_SoN(ris, A);

  (*sign)=1;
  for(i=0; i<NCOLOR; i++)
     {
     big=0.0;
     for(j=0; j<NCOLOR; j++) 
        {
        temp = fabs(ris->comp[i][j]);
        if( temp>big ) big=temp;
        }
     vv[i]=1.0/big;
     }

  for(j=0; j<NCOLOR; j++)
     {
     for(i=0; i<j; i++)
        {
        sum=ris->comp[i][j];
        for(k=0; k<i; k++)
           { 
           sum-=(ris->comp[i][k])*(ris->comp[k][j]);
           }
        ris->comp[i][j]=sum;
        }
     
     big=0.0;
     for(i=j; i<NCOLOR; i++)
        {
        sum=ris->comp[i][j];
        for(k=0; k<j; k++)
           {
           sum-=(ris->comp[i][k])*(ris->comp[k][j]);
           }
        ris->comp[i][j]=sum;
       
        temp = vv[i]*fabs(sum);
        if(temp >= big)
          {
          big=temp;
          imax=i;
          }
        }

     if(j!= imax)
       {
       for(k=0; k<NCOLOR; k++)
          {
          dum=ris->comp[imax][k];
          ris->comp[imax][k]=ris->comp[j][k];
          ris->comp[j][k]=dum;
          }
       (*sign)*=(-1);
       vv[imax]=vv[j];
       }
 
     if(j!= NCOLOR-1)
       {
       dum=1.0/(ris->comp[j][j]);
       for(i=j+1; i<NCOLOR; i++)
          { 
          (ris->comp[i][j])*=dum;
          }
       }
     }
  }

/* determinant */
double det_SoN(SoN const *__restrict__ const A)
  {
  int i;
  double ris;
  SoN lu; 

  LU_SoN(A, &lu, &i); 

  if(i>0)
    {
    ris=1.0;
    }
  else
   {
   ris=-1.0;
   }

  for(i=0; i<NCOLOR; i++)
     {
     ris*=(lu.comp[i][i]);
     }

  return ris;
  }

/* gives 0 if the matrix is in SU(N) and 1 otherwise */
int scheck_SoN(SoN const *__restrict__ const A)
  {
  int i, j, k, ris;
  double aux;

  ris=0;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           aux+=(A->comp[i][k])*(A->comp[j][k]);
           }
        if(i==j) aux-=1.0;
        if(fabs(aux)>MIN_VALUE) ris=1;
        }
     }

  if(ris==0)
    {
    if(det_SoN(A)<-0.5)
      {
      ris=1;
      }
    }

  return ris;
  }


/* sunitarize, first part: Gram–Schmidt process and check for det=+1 */
void unitarize_aux_SoN(SoN *__restrict__ A)
  {
  int i, j, k;
  double norm, c[NCOLOR];

  for(i=0; i<NCOLOR; i++) /* loop on the lines */
     {
     /* scalar product with previous lines */
     for(j=0; j<i; j++) 
        {
        c[j]=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           c[j]+=(A->comp[i][k])*(A->comp[j][k]);
           }
        }
     
     /* orthogonalize with respect to previous lines */
     for(j=0; j<i; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           A->comp[i][k]-=c[j]*(A->comp[j][k]);
           }
        }
     
     /* normalizze the line */
     norm=0.0;
     for(k=0; k<NCOLOR; k++)
        {
        norm+=(A->comp[i][k])*(A->comp[i][k]);
        }
     norm=1.0/sqrt(norm);
     for(k=0; k<NCOLOR; k++)
        {
        (A->comp[i][k])*=norm;
        }
     }

  norm=det_SoN(A);  /* norm =+1 o -1 */

  if(norm<=-0.5)
    {
    for(i=0;i<NCOLOR;i++) 
       {
       (A->comp[NCOLOR-1][i])*=(-1.0);
       }
    } 
  }


/* sunitarize */
void unitarize_SoN(SoN *__restrict__ A)
  {
  while(scheck_SoN(A)!=0) unitarize_aux_SoN(A);
  }

void print_on_screen_SoN(SoN const *__restrict__ const A)
  {
  int i, j;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        printf("%.16f ", A->comp[i][j]);
        }
     }
  printf("\n");
  }

void print_on_file_SoN(FILE *fp, SoN const *__restrict__ const A)
  {
  int i, j;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        fprintf(fp, "%.16f ", A->comp[i][j]);
        }
     }
  fprintf(fp, "\n");
  }

void read_from_file_SoN(FILE *fp, SoN *__restrict__ A)
  {
  int i, j, err;
  double data;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fscanf(fp, "%lg ", &data);
        if(err!=1)
          {
          fprintf(stderr, "Problems reading SoN matrix from file");
          }
        A->comp[i][j]=data; 
        }
     }
  }
  

#endif
