#ifndef SUN_C
#define SUN_C

#include"../Macro/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<string.h>

#include"../Endian/endianness.h"
#include"../Rng/random.h"
#include"sun.h"

/* A=1 */
void one_SuN(SuN *A)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=i+1; j<Ncolor; j++)
        {
        A->comp[i][j]=0.0+0.0*I;
        A->comp[j][i]=0.0+0.0*I;
        }
     A->comp[i][i]=1.0+0.0*I;
     }
  }

/* A=0 */
void zero_SuN(SuN *A)       
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]=0.0+0.0*I;
        }
     }
  }

/* A=B */
void equal_SuN(SuN *A, SuN const * const B)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]=B->comp[i][j];
        }
     }
  }

/* A=B^{dag} */
void equal_dag_SuN(SuN *A, SuN const * const B)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]=conj(B->comp[j][i]);
        }
     }
  }


/* A+=B */
void plus_equal_SuN(SuN *A, SuN const * const B) 
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]+=B->comp[i][j];
        }
     }
  }


/* A+=B^{dag} */
void plus_equal_dag_SuN(SuN *A, SuN const * const B)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]+=conj(B->comp[j][i]);
        }
     }
  }

/* A-=B */
void minus_equal_SuN(SuN *A, SuN const * const B)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]-=B->comp[i][j];
        }
     }
  }

/* A-=B^{dag} */
void minus_equal_dag_SuN(SuN *A, SuN const * const B)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]-=conj(B->comp[j][i]);
        }
     }
  }

/* A=b*B+c*C */
void lin_comb_SuN(SuN *A, double b, SuN const * const B, double c, SuN const * const C)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]=b*(B->comp[i][j])+c*(C->comp[i][j]);
        }
     }
  }

/* A=b*B^{dag}+c*C */
void lin_comb_dag1_SuN(SuN *A, double b, SuN const * const B, double c, SuN const * const C)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]=b*conj(B->comp[j][i])+c*(C->comp[i][j]);
        }
     }
  }

/* A=b*B+c*C^{dag} */
void lin_comb_dag2_SuN(SuN *A, double b, SuN const * const B, double c, SuN const * const C)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]=b*(B->comp[i][j])+c*conj(C->comp[j][i]);
        }
     }
  }

/* A=b*B^{dag}+c*C^{dag} */
void lin_comb_dag12_SuN(SuN *A, double b, SuN const * const B, double c, SuN const * const C)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]=b*conj(B->comp[j][i])+c*conj(C->comp[j][i]);
        }
     }
  }

/* A*=r */
void times_equal_real_SuN(SuN *A, double r)
  {
  int i, j;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        A->comp[i][j]*=(r+0.0*I);
        }
     }
  }

/* A*=B */
void times_equal_SuN(SuN *__restrict__ A, SuN const *__restrict__ const B)
   {
   #if Ncolor != 3
     int i, j, k;
     double complex sum, aux[Ncolor];

     for(i=0; i<Ncolor; i++)
        {
        for(j=0; j<Ncolor; j++)
           {
           aux[j]=A->comp[i][j];
           }
     
        for(j=0; j<Ncolor; j++)
           {
           sum=0.0+0.0*I;
           for(k=0; k<Ncolor; k++)
              {
              sum+=aux[k]*(B->comp[k][j]);
              }
           A->comp[i][j]=sum;
           }
        }
   #else
     int i, j;
     register double complex sum, aux0, aux1, aux2;

     for(i=0; i<3; i++)
        {
        aux0=A->comp[i][0];
        aux1=A->comp[i][1];
        aux2=A->comp[i][2];

        for(j=0; j<3; j++)
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
   #if Ncolor != 3
     int i, j, k;
     double complex sum, aux[Ncolor];

     for(i=0; i<Ncolor; i++)
        {
        for(j=0; j<Ncolor; j++)
           {
           aux[j]=A->comp[i][j];
           }
     
        for(j=0; j<Ncolor; j++)
           {
           sum=0.0+0.0*I;
           for(k=0; k<Ncolor; k++)
              {
              sum+=aux[k]*conj(B->comp[j][k]);
              }
           A->comp[i][j]=sum;
           }
        }
   #else
     int i, j;
     register double complex sum, aux0, aux1, aux2;

     for(i=0; i<3; i++)
        {
        aux0=A->comp[i][0];
        aux1=A->comp[i][1];
        aux2=A->comp[i][2];

        for(j=0; j<3; j++)
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
void times_SuN(SuN * A, SuN const * const B, SuN const * const C)
  {
  #if Ncolor != 3
    int i, j, k;
    double complex sum;
  
    for(i=0; i<Ncolor; i++)
       {
       for(j=0; j<Ncolor; j++)
          {
          sum=0.0+0.0*I;
          for(k=0; k<Ncolor; k++)
             {
             sum+=(B->comp[i][k])*(C->comp[k][j]);
             }
          A->comp[i][j]=sum;
          }
       } 
  #else
    int i, j;
    register double complex aux0, aux1, aux2, sum;

    for(i=0; i<3; i++)
       {
       aux0=B->comp[i][0];
       aux1=B->comp[i][1];
       aux2=B->comp[i][2];
       for(j=0; j<3; j++)
          {
          sum =aux0*(C->comp[0][j]);
          sum+=aux1*(C->comp[1][j]);
          sum+=aux2*(C->comp[2][j]);

          A->comp[i][j]=sum;
          }
       } 
  #endif
  }

/* A=B^{dag}*C */
void times_dag1_SuN(SuN *A, SuN const * const B, SuN const * const C) 
  {
  #if Ncolor != 3
    int i, j, k;
    double complex sum;
  
    for(i=0; i<Ncolor; i++)
       {
       for(j=0; j<Ncolor; j++)
          {
          sum=0.0+0.0*I;
          for(k=0; k<Ncolor; k++)
             {
             sum+=conj(B->comp[k][i])*(C->comp[k][j]);
             }
          A->comp[i][j]=sum;
          }
       } 
  #else
    int i, j;
    register double complex aux0, aux1, aux2, sum;

    for(i=0; i<3; i++)
       {
       aux0=conj(B->comp[0][i]);
       aux1=conj(B->comp[1][i]);
       aux2=conj(B->comp[2][i]);
       for(j=0; j<3; j++)
          {
          sum =aux0*(C->comp[0][j]);
          sum+=aux1*(C->comp[1][j]);
          sum+=aux2*(C->comp[2][j]);

          A->comp[i][j]=sum;
          }
       } 
  #endif
  }

/* A=B*C^{dag} */
void times_dag2_SuN(SuN *A, SuN const * const B, SuN const * const C) 
  {
  #if Ncolor != 3
    int i, j, k;
    double complex sum;
  
    for(i=0; i<Ncolor; i++)
       {
       for(j=0; j<Ncolor; j++)
          {
          sum=0.0+0.0*I;
          for(k=0; k<Ncolor; k++)
             {
             sum+=(B->comp[i][k])*conj(C->comp[j][k]);
             }
          A->comp[i][j]=sum;
          }
       } 
  #else
    int i, j;
    register double complex aux0, aux1, aux2, sum;

    for(i=0; i<3; i++)
       {
       aux0=B->comp[i][0];
       aux1=B->comp[i][1];
       aux2=B->comp[i][2];
       for(j=0; j<3; j++)
          {
          sum =aux0*conj(C->comp[j][0]);
          sum+=aux1*conj(C->comp[j][1]);
          sum+=aux2*conj(C->comp[j][2]);

          A->comp[i][j]=sum;
          }
       } 
  #endif
  }

/* A=B^{dag}*C^{dag} */
void times_dag12_SuN(SuN *A, SuN const * const B, SuN const * const C)
  {
  #if Ncolor != 3
  int i, j, k;
  double complex sum;
  
  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<Ncolor; k++)
           {
           sum+=conj(B->comp[k][i])*conj(C->comp[j][k]);
           }
        A->comp[i][j]=sum;
        }
     } 
  #else
    int i, j;
    register double complex aux0, aux1, aux2, sum;

    for(i=0; i<3; i++)
       {
       aux0=conj(B->comp[0][i]);
       aux1=conj(B->comp[1][i]);
       aux2=conj(B->comp[2][i]);
       for(j=0; j<3; j++)
          {
          sum =aux0*conj(C->comp[j][0]);
          sum+=aux1*conj(C->comp[j][1]);
          sum+=aux2*conj(C->comp[j][2]);

          A->comp[i][j]=sum;
          }
       } 
  #endif
  }

/* SU(N) random matrix
   generated a la Cabibbo Marinari with N(N-1)/2 SU(2) random matrices */
void rand_matrix_SuN(SuN *A)
  {
  int i, j, k;
  double p0, p1, p2, p3, p;
  double complex aux00, aux01, aux10, aux11, temp0, temp1;

  one_SuN(A);

  for(i=0; i<Ncolor-1; i++)
     {
     for(j=i+1; j<Ncolor; j++)
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

        aux00= p0 + p3*I;
        aux01= p2 + p1*I;
        aux10=-p2 + p1*I;
        aux11= p0 - p3*I;

        for(k=0; k<Ncolor; k++)
           {
           temp0=A->comp[k][i]*aux00 + A->comp[k][j]*aux10;
           temp1=A->comp[k][i]*aux01 + A->comp[k][j]*aux11;
           A->comp[k][i]=temp0;
           A->comp[k][j]=temp1;
           }
        }
     }
  }

/* l2 norm of the matrix */
double norm_SuN(SuN const * const A)
  {
  int i, j;
  double aux, ris;

  ris=0.0;
  for(i=0; i<Ncolor; i++)
     { 
     for(j=0; j<Ncolor; j++)
        {
        aux=cabs(A->comp[i][j]);
        ris+=aux*aux;
        }
     }
  return sqrt(ris);
  }

/* real part of the trace /N */
double retr_SuN(SuN const * const A)
  {
  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<Ncolor; i++)
     {
     tr+=A->comp[i][i];
     }
  ris=creal(tr)/(double)Ncolor;
  return ris;
  }

/* imaginary part of the trace /N */
double imtr_SuN(SuN const * const A)
  {
  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<Ncolor; i++)
     {
     tr+=A->comp[i][i];
     }
  ris=cimag(tr)/(double)Ncolor;
  return ris;
  }

/* LU decomposition with partial pivoting
   from Numerical Recipes in C, pag 46 */
void LU_SuN(SuN const * const A, SuN *ris, int *sign) 
  {
  int i, imax, j, k;
  double  big, temp;
  double complex sum, dum;
  double vv[Ncolor];

  imax=0;
  equal_SuN(ris, A);

  (*sign)=1;
  for(i=0; i<Ncolor; i++)
     {
     big=0.0;
     for(j=0; j<Ncolor; j++) 
        {
        temp = cabs(ris->comp[i][j]);
        if( temp>big ) big=temp;
        }
     vv[i]=1.0/big;
     }

  for(j=0; j<Ncolor; j++)
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
     for(i=j; i<Ncolor; i++)
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
       for(k=0; k<Ncolor; k++)
          {
          dum=ris->comp[imax][k];
          ris->comp[imax][k]=ris->comp[j][k];
          ris->comp[j][k]=dum;
          }
       (*sign)*=(-1);
       vv[imax]=vv[j];
       }
 
     if(j!= Ncolor-1)
       {
       dum=(1.0+0.0*I)/(ris->comp[j][j]);
       for(i=j+1; i<Ncolor; i++)
          { 
          (ris->comp[i][j])*=dum;
          }
       }
     }
  }

/* determinant */
complex double det_SuN(SuN const * const A)
  {
  #if Ncolor==3
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

    for(i=0; i<Ncolor; i++)
       {
       ris*=(lu.comp[i][i]);
       }

    return ris;
  #endif
  }

/* gives 0 if the matrix is in SU(N) and 1 otherwise */
int scheck_SuN(SuN const * const A)
  {
  int i, j, k, ris;
  double complex aux;

  ris=0;
 
  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        aux=0.0+0.0*I;
        for(k=0; k<Ncolor; k++)
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
void unitarize_aux_SuN(SuN *A)
  {
  int i, j, k;
  double complex c[Ncolor];
  double norm;

  for(i=0; i<Ncolor; i++) /* loop on the lines */
     {
     /* scalar product with previous lines */
     for(j=0; j<i; j++) 
        {
        c[j]=0.0+0.0*I;
        for(k=0; k<Ncolor; k++)
           {
           c[j]+=(A->comp[i][k])*conj(A->comp[j][k]);
           }
        }
     
     /* orthogonalize with respect to previous lines */
     for(j=0; j<i; j++)
        {
        for(k=0; k<Ncolor; k++)
           {
           A->comp[i][k]-=c[j]*(A->comp[j][k]);
           }
        }
     
     /* normalizze the line */
     norm=0.0;
     for(k=0; k<Ncolor; k++)
        {
        norm+=(creal(A->comp[i][k])*creal(A->comp[i][k])
              +cimag(A->comp[i][k])*cimag(A->comp[i][k]));
        }
     norm=1.0/sqrt(norm);
     for(k=0; k<Ncolor; k++)
        {
        (A->comp[i][k])*=norm;
        }
     }

  norm=creal(det_SuN(A));  /* norm =+1 o -1 */

  if(norm<=-0.5)
    {
    for(i=0;i<Ncolor;i++) 
       {
       (A->comp[Ncolor-1][i])*=(-1.0);
       }
    } 
  }


/* sunitarize */
void unitarize_SuN(SuN *A)
  {
  while(scheck_SuN(A)!=0) unitarize_aux_SuN(A);
  }

void print_on_screen_SuN(SuN const * const A)
  {
  int i, j;
 
  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        printf("%.16f %.16f ", creal(A->comp[i][j]), cimag(A->comp[i][j]));
        }
     }
  printf("\n");
  }

void print_on_file_SuN(FILE *fp, SuN const * const A)
  {
  int i, j;
 
  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        fprintf(fp, "%.16f %.16f ", creal(A->comp[i][j]), cimag(A->comp[i][j]));
        }
     }
  fprintf(fp, "\n");
  }

void print_on_binary_file_SuN(FILE *fp, SuN const * const A)
  {
  int i, j;
  double re, im;
 
  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        re=creal(A->comp[i][j]);
        im=cimag(A->comp[i][j]);

        fwrite(&re, 1, sizeof(double), fp);  
        fwrite(&im, 1, sizeof(double), fp);  
        }
     }
  }

void read_from_file_SuN(FILE *fp, SuN *A)
  {
  int i, j, err;
  double re, im;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
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

void read_from_binary_file_SuN(FILE *fp, SuN *A)
  {
  int i, j, err;
  double re, im;
  double aux[2];
 
  err=0;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        err+=fread(&re, sizeof(double), 1, fp);
        err+=fread(&im, sizeof(double), 1, fp);
        aux[0]=re;
        aux[1]=im;

        memcpy((void *)&(A->comp[i][j]), (void*)aux, sizeof(aux));
        //equivalent to A->comp[i][j]=re+im*I;
        }
     }

  if(err!=2*Ncolor*Ncolor)
    {
    fprintf(stderr, "Problems reading SuN matrix from file");
    }
  }
 

/* read from binary file and change endianness */ 
void read_from_binary_file_swap_SuN(FILE *fp, SuN *A)
  {
  int i, j, err;
  double re, im;
  double aux[2];
  
  err=0;

  for(i=0; i<Ncolor; i++)
     {
     for(j=0; j<Ncolor; j++)
        {
        err+=fread(&re, sizeof(double), 1, fp);
        err+=fread(&im, sizeof(double), 1, fp);
        SwapBytesDouble(&re);
        SwapBytesDouble(&im);
        aux[0]=re;
        aux[1]=im;

        memcpy((void *)&(A->comp[i][j]), (void*)aux, sizeof(aux));
        // equivalent to A->comp[i][j]=re+im*I;
        }
     }

  if(err!=2*Ncolor*Ncolor)
    {
    fprintf(stderr, "Problems reading SuN matrix from file");
    }
  }
 
 
#endif
