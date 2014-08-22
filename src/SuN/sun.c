#ifndef SUN_C
#define SUN_C

#include<complex.h>
#include<math.h>
#include<stdio.h>

#include"../Macro/macro.h"
#include"../Rng/random.h"
#include"sun.h"

/* A=1 */
void one_SuN(SuN *__restrict__ A)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        A->comp[i][j]=0.0+0.0*I;
        A->comp[j][i]=0.0+0.0*I;
        }
     A->comp[i][i]=1.0+0.0*I;
     }
  }

/* A=0 */
void zero_SuN(SuN *__restrict__ A)       
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=0.0+0.0*I;
        }
     }
  }

/* A=B */
void equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B)
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
void equal_dag_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=conj(B->comp[j][i]);
        }
     }
  }


/* A+=B */
void plus_equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B) 
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
void plus_equal_dag_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]+=conj(B->comp[j][i]);
        }
     }
  }

/* A-=B */
void minus_equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B)
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
void minus_equal_dag_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]-=conj(B->comp[j][i]);
        }
     }
  }

/* A=b*B+c*C */
void lin_comb_SuN(SuN *__restrict__ A, 
                  double b, 
                  SuN const *__restrict__ const B, 
                  double c, 
                  SuN const *__restrict__ const C)
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
void lin_comb_dag1_SuN(SuN *__restrict__ A, 
                       double b, 
                       SuN const *__restrict__ const B, 
                       double c, 
                       SuN const *__restrict__ const C)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=b*conj(B->comp[j][i])+c*(C->comp[i][j]);
        }
     }
  }

/* A=b*B+c*C^{dag} */
void lin_comb_dag2_SuN(SuN *__restrict__ A, 
                       double b, 
                       SuN const *__restrict__ const B, 
                       double c, 
                       SuN const *__restrict__ const C)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=b*(B->comp[i][j])+c*conj(C->comp[j][i]);
        }
     }
  }

/* A=b*B^{dag}+c*C^{dag} */
void lin_comb_dag12_SuN(SuN *__restrict__ A, 
                        double b, 
                        SuN const *__restrict__ const B, 
                        double c, 
                        SuN const *__restrict__ const C)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]=b*conj(B->comp[j][i])+c*conj(C->comp[j][i]);
        }
     }
  }

/* A*=r */
void times_equal_real_SuN(SuN *__restrict__ A, double r)
  {
  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[i][j]*=(r+0.0*I);
        }
     }
  }

/* A*=B */
void times_equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B)
   {
   #if NCOLOR != 3
   int i, j, k;
   double complex sum, aux[NCOLOR];

   for(i=0; i<NCOLOR; i++)
      {
      for(j=0; j<NCOLOR; j++)
         {
         aux[j]=A->comp[i][j];
         }
     
      for(j=0; j<NCOLOR; j++)
         {
         sum=0.0+0.0*I;
         for(k=0; k<NCOLOR; k++)
            {
            sum+=aux[k]*(B->comp[k][j]);
            }
         A->comp[i][j]=sum;
         }
      }
   #else
   int i, j;
   register double complex sum, aux0, aux1, aux2;

   for(i=0; i<NCOLOR; i++)
      {
      aux0=A->comp[i][0];
      aux1=A->comp[i][1];
      aux2=A->comp[i][2];

      for(j=0; j<NCOLOR; j++)
         {
         sum=aux0*(B->comp[0][j]);
         sum+=aux1*(B->comp[1][j]);
         sum+=aux2*(B->comp[2][j]);

         A->comp[i][j]=sum;
         }
      }
   #endif
   }
 
/* A*=B^{dag} */
void times_equal_dag_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B)
   {
   #if NCOLOR != 3
   int i, j, k;
   double complex sum, aux[NCOLOR];

   for(i=0; i<NCOLOR; i++)
      {
      for(j=0; j<NCOLOR; j++)
         {
         aux[j]=A->comp[i][j];
         }
     
      for(j=0; j<NCOLOR; j++)
         {
         sum=0.0+0.0*I;
         for(k=0; k<NCOLOR; k++)
            {
            sum+=aux[k]*conj(B->comp[j][k]);
            }
         A->comp[i][j]=sum;
         }
      }
   #else
   int i, j;
   register double complex sum, aux0, aux1, aux2;

   for(i=0; i<NCOLOR; i++)
      {
      aux0=A->comp[i][0];
      aux1=A->comp[i][1];
      aux2=A->comp[i][2];

      for(j=0; j<NCOLOR; j++)
         {
         sum=aux0*conj(B->comp[j][0]);
         sum+=aux1*conj(B->comp[j][1]);
         sum+=aux2*conj(B->comp[j][2]);

         A->comp[i][j]=sum;
         }
      }
   #endif
   }

/* A=B*C */
void times_SuN(SuN *__restrict__ A, 
               SuN const *__restrict__ const B, 
               SuN const *__restrict__ const C)
  {
  int i, j, k;
  double complex sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[i][k])*(C->comp[k][j]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B^{dag}*C */
void times_dag1_SuN(SuN *__restrict__ A, 
                    SuN const *__restrict__ const B, 
                    SuN const *__restrict__ const C) 
  {
  int i, j, k;
  double complex sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=conj(B->comp[k][i])*(C->comp[k][j]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B*C^{dag} */
void times_dag2_SuN(SuN *__restrict__ A, 
                    SuN const *__restrict__ const B, 
                    SuN const *__restrict__ const C) 
  {
  int i, j, k;
  double complex sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[i][k])*conj(C->comp[j][k]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B^{dag}*C^{dag} */
void times_dag12_SuN(SuN *__restrict__ A, 
                     SuN const *__restrict__ const B, 
                     SuN const *__restrict__ const C)
  {
  int i, j, k;
  double complex sum;
  
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=conj(B->comp[k][i])*conj(C->comp[j][k]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* SU(N) random matrix
   generated a la Cabibbo Marinari with N(N-1)/2 SU(2) random matrices */
void rand_matrix_SuN(SuN *__restrict__ A)
  {
  int i, j, k;
  double p0, p1, p2, p3, p;
  double complex aux[2][2], temp[2];

  one_SuN(A);

  for(i=0; i<NCOLOR-1; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        /* SU(2) random components */
        p=2.0;
        while(p>1.0)
             { 
             p0=1.0-2.0*casuale();
             p1=1.0-2.0*casuale();
             p2=1.0-2.0*casuale();
             p3=1.0-2.0*casuale();
             p=sqrt(p0*p0+p1*p1+p2*p2+p3*p3);
             }

        p0/=p;
        p1/=p;
        p2/=p;
        p3/=p;

        aux[0][0]= p0 + p3*I;
        aux[0][1]= p2 + p1*I;
        aux[1][0]=-p2 + p1*I;
        aux[1][1]= p0 - p3*I;

        for(k=0; k<NCOLOR; k++)
           {
           temp[0]=A->comp[k][i]*aux[0][0] + A->comp[k][j]*aux[1][0];
           temp[1]=A->comp[k][i]*aux[0][1] + A->comp[k][j]*aux[1][1];
           A->comp[k][i]=temp[0];
           A->comp[k][j]=temp[1];
           }
        }
     }
  }

/* l2 norm of the matrix */
double norm_SuN(SuN const *__restrict__ const A)
  {
  int i, j;
  double aux, ris;

  ris=0.0;
  for(i=0; i<NCOLOR; i++)
     { 
     for(j=0; j<NCOLOR; j++)
        {
        aux=cabs(A->comp[i][j]);
        ris+=aux*aux;
        }
     }
  return sqrt(ris);
  }

/* real part of the trace /N */
double retr_SuN(SuN const *__restrict__ const A)
  {
  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[i][i];
     }
  ris=creal(tr)/(double)NCOLOR;
  return ris;
  }

/* imaginary part of the trace /N */
double imtr_SuN(SuN const *__restrict__ const A)
  {
  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[i][i];
     }
  ris=cimag(tr)/(double)NCOLOR;
  return ris;
  }

/* LU decomposition with partial pivoting
   from Numerical Recipes in C, pag 46 */
void LU_SuN(SuN const *__restrict__ const A, SuN *__restrict__ ris, int *sign) 
  {
  int i, imax, j, k;
  double  big, temp;
  double complex sum, dum;
  double vv[NCOLOR];

  imax=0;
  equal_SuN(ris, A);

  (*sign)=1;
  for(i=0; i<NCOLOR; i++)
     {
     big=0.0;
     for(j=0; j<NCOLOR; j++) 
        {
        temp = cabs(ris->comp[i][j]);
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
       
        temp = vv[i]*cabs(sum);
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
       dum=(1.0+0.0*I)/(ris->comp[j][j]);
       for(i=j+1; i<NCOLOR; i++)
          { 
          (ris->comp[i][j])*=dum;
          }
       }
     }
  }

/* determinant */
complex double det_SuN(SuN const *__restrict__ const A)
  {
  #if NCOLOR==3
    complex double ris=0.0+0.0*I;
 
    ris+=(A->comp[0][0])*(A->comp[1][1])*(A->comp[2][2]);
    ris+=(A->comp[1][0])*(A->comp[2][1])*(A->comp[0][2]);
    ris+=(A->comp[2][0])*(A->comp[0][1])*(A->comp[1][2]);
    ris-=(A->comp[2][0])*(A->comp[1][1])*(A->comp[0][2]);
    ris-=(A->comp[1][0])*(A->comp[0][1])*(A->comp[2][2]);
    ris-=(A->comp[0][0])*(A->comp[2][1])*(A->comp[1][2]);

    return ris;
  #else
    int i;
    double complex ris;
    SuN lu; 

    LU_SuN(A, &lu, &i); 

    if(i>0)
      {
      ris=1.0+0.0*I;
      }
    else
     {
     ris=-1.0+0.0*I;
     }

    for(i=0; i<NCOLOR; i++)
       {
       ris*=(lu.comp[i][i]);
       }

    return ris;
  #endif
  }

/* gives 0 if the matrix is in SU(N) and 1 otherwise */
int scheck_SuN(SuN const *__restrict__ const A)
  {
  int i, j, k, ris;
  double complex aux;

  ris=0;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           aux+=(A->comp[i][k])*conj(A->comp[j][k]);
           }
        if(i==j) aux-=(1.0+0.0*I);
        if(cabs(aux)>MIN_VALUE) ris=1;
        }
     }

  if(ris==0)
    {
    if(creal(det_SuN(A))<-0.5)
      {
      ris=1;
      }
    }

  return ris;
  }


/* sunitarize, first part: Gram–Schmidt process and check for det=+1 */
void unitarize_aux_SuN(SuN *__restrict__ A)
  {
  int i, j, k;
  double complex c[NCOLOR];
  double norm;

  for(i=0; i<NCOLOR; i++) /* loop on the lines */
     {
     /* scalar product with previous lines */
     for(j=0; j<i; j++) 
        {
        c[j]=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           c[j]+=(A->comp[i][k])*conj(A->comp[j][k]);
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
        norm+=(creal(A->comp[i][k])*creal(A->comp[i][k])
              +cimag(A->comp[i][k])*cimag(A->comp[i][k]));
        }
     norm=1.0/sqrt(norm);
     for(k=0; k<NCOLOR; k++)
        {
        (A->comp[i][k])*=norm;
        }
     }

  norm=creal(det_SuN(A));  /* norm =+1 o -1 */

  if(norm<=-0.5)
    {
    for(i=0;i<NCOLOR;i++) 
       {
       (A->comp[NCOLOR-1][i])*=(-1.0);
       }
    } 
  }


/* sunitarize */
void unitarize_SuN(SuN *__restrict__ A)
  {
  while(scheck_SuN(A)!=0) unitarize_aux_SuN(A);
  }

void print_on_screen_SuN(SuN const *__restrict__ const A)
  {
  int i, j;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        printf("%.16f %.16f ", creal(A->comp[i][j]), cimag(A->comp[i][j]));
        }
     }
  printf("\n");
  }

void print_on_file_SuN(FILE *fp, SuN const *__restrict__ const A)
  {
  int i, j;
 
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        fprintf(fp, "%.16f %.16f ", creal(A->comp[i][j]), cimag(A->comp[i][j]));
        }
     }
  fprintf(fp, "\n");
  }

void read_from_file_SuN(FILE *fp, SuN *__restrict__ A)
  {
  int i, j, err;
  double re, im;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        err=fscanf(fp, "%lg %lg", &re, &im);
        if(err!=2)
          {
          fprintf(stderr, "Problems reading SuN matrix from file");
          }
        A->comp[i][j]=re+im*I; 
        }
     }
  }
  

#endif
