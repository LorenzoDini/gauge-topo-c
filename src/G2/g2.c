#ifndef G2_C
#define G2_C

/* for definitions see
   Greensite, Langfeld, Olejnik, Reinhardt, Tok   Phys Rev D 75, p.034501 (2007) */

#include<math.h>
#include<stdio.h>

#include"g2.h"
#include"../Macro/macro.h"
#include"../Rng/random.h"

/* A=1 */
void one_G2(G2 *__restrict__ A)
   {
  int i, j;

   for(i=0; i<7; i++)
      {
      for(j=i+1; j<7; j++)
         {
         A->comp[i][j]=0.0;
         A->comp[j][i]=0.0;
         }
      A->comp[i][i]=1.0;
      }
   }

/* A=0 */
void zero_G2(G2 *__restrict__ A)
   {
   int i, j;
   
   for(i=0; i<7; i++)
      {
      for(j=0; j<7; j++)
         {
         A->comp[i][j]=0.0; 
         }
      }
   }

/* A=B */
void equal_G2(G2 *__restrict__ A, G2 const *__restrict__ const B)
  { 
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]=B->comp[i][j];
        }
     }
  }

/* A=B^{dag} */
void equal_dag_G2(G2 *__restrict__ A, G2 const *__restrict__ const B)
  { 
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]=B->comp[j][i];
        }
     }
  }

/* A+=B */
void plus_equal_G2(G2 *__restrict__ A, G2 const *__restrict__ const B)
  { 
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]+=B->comp[i][j];
        }
     }
  }

/* A+=B^{dag} */
void plus_equal_dag_G2(G2 *__restrict__ A, G2 const *__restrict__ const B)
  { 
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]+=B->comp[j][i];
        }
     }
  }

/* A-=B */
void minus_equal_G2(G2 *__restrict__ A, G2 const *__restrict__ const B)
  { 
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]-=B->comp[i][j];
        }
     }
  }

/* A-=B^{dag} */
void minus_equal_dag_G2(G2 *__restrict__ A, G2 const *__restrict__ const B)
  { 
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]-=B->comp[j][i];
        }
     }
  }

/* A=b*B+c*C */
void lin_comb_G2(G2 *__restrict__ A, double b, G2 const *__restrict__ const B, double c, G2 const *__restrict__ const C)
  {
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]=b*(B->comp[i][j])+c*(C->comp[i][j]);
        }
     }
  }

/* A=b*B^{dag}+c*C */
void lin_comb_dag1_G2(G2 *__restrict__ A, double b, G2 const *__restrict__ const B, double c, G2 const *__restrict__ const C)  
  {
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]=b*(B->comp[j][i])+c*(C->comp[i][j]);
        }
     }
  }

/* A=b*B+c*C^{dag} */
void lin_comb_dag2_G2(G2 *__restrict__ A, double b, G2 const *__restrict__ const B, double c, G2 const *__restrict__ const C) 
  {
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]=b*(B->comp[i][j])+c*(C->comp[j][i]);
        }
     }
  }

/* A=b*B^{dag}+c*C^{dag} */
void lin_comb_dag12_G2(G2 *__restrict__ A, double b, G2 const *__restrict__ const B, double c, G2 const *__restrict__ const C)
  {
  int i, j;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        A->comp[i][j]=b*(B->comp[j][i])+c*(C->comp[j][i]);
        }
     }
  }

/* A*=r */
void times_equal_real_G2(G2 *__restrict__ A, double r)
  {
  int i, j;

  for(i=0; i<7; i++)
     {  
     for(j=0; j<7; j++)
        {
        (A->comp[i][j])*=r;
        }
     }
  }

/* A*=B */
/*
void times_equal_G2(G2 * A, G2 const * const B)
  {
  int i, j, k; 
  double sum, aux[7];

  for(i=0; i< 7; i++)
     {
     for(j=0; j<7; j++)
        {
        aux[j]=A->comp[i][j];
        }

     for(j=0; j<7; j++)
        {
        sum=0.0;
        for(k=0; k<7; k++)
           { 
           sum+=aux[k]*(B->comp[k][j]);
           } 

        A->comp[i][j]=sum;
        }
     }
  }
*/
void times_equal_G2(G2 *__restrict__ A, G2 const *__restrict__ const B)
  {
  int i, j; 
  register double sum, aux0, aux1, aux2, aux3, aux4, aux5, aux6;

  for(i=0; i<7; i++)
     {
     aux0=A->comp[i][0];
     aux1=A->comp[i][1];
     aux2=A->comp[i][2];
     aux3=A->comp[i][3];
     aux4=A->comp[i][4];
     aux5=A->comp[i][5];
     aux6=A->comp[i][6];

     for(j=0; j<7; j++)
        {
        sum=aux0*(B->comp[0][j]);
        sum+=aux1*(B->comp[1][j]);
        sum+=aux2*(B->comp[2][j]);
        sum+=aux3*(B->comp[3][j]);
        sum+=aux4*(B->comp[4][j]);
        sum+=aux5*(B->comp[5][j]);
        sum+=aux6*(B->comp[6][j]);

        A->comp[i][j]=sum;
        }
     }
  }

/* A*=B^{dag} */
/*
void times_equal_dag_G2(G2 *A, G2 const * const B)
  {
  int i, j, k;
  double sum, aux[7];

  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        aux[j]=A->comp[i][j];
        }
     
     for(j=0; j<7; j++)
        {
        sum=0.0;
        for(k=0; k<7; k++)
           {
           sum+=aux[k]*(B->comp[j][k]);
           }
        A->comp[i][j]=sum;
        }
     }
  }
*/
void times_equal_dag_G2(G2 *__restrict__ A, G2 const *__restrict__ const B)
  {
  int i, j; 
  register double sum, aux0, aux1, aux2, aux3, aux4, aux5, aux6;

  for(i=0; i<7; i++)
     {
     aux0=A->comp[i][0];
     aux1=A->comp[i][1];
     aux2=A->comp[i][2];
     aux3=A->comp[i][3];
     aux4=A->comp[i][4];
     aux5=A->comp[i][5];
     aux6=A->comp[i][6];

     for(j=0; j<7; j++)
        {
        sum=aux0*(B->comp[j][0]);
        sum+=aux1*(B->comp[j][1]);
        sum+=aux2*(B->comp[j][2]);
        sum+=aux3*(B->comp[j][3]);
        sum+=aux4*(B->comp[j][4]);
        sum+=aux5*(B->comp[j][5]);
        sum+=aux6*(B->comp[j][6]);

        A->comp[i][j]=sum;
        }
     }

  }

/* A=B*C */
void times_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, G2 const *__restrict__ const C)    
  {
  int i, j, k;
  double sum;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        sum=0.0;
        for(k=0; k<7; k++)
           {
           sum+=(B->comp[i][k])*(C->comp[k][j]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B^{dag}*C */
void times_dag1_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, G2 const *__restrict__ const C)
  {
  int i, j, k;
  double sum;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        sum=0.0;
        for(k=0; k<7; k++)
           {
           sum+=(B->comp[k][i])*(C->comp[k][j]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B*C^{dag} */
void times_dag2_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, G2 const *__restrict__ const C) 
  {
  int i, j, k;
  double sum;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        sum=0.0;
        for(k=0; k<7; k++)
           {
           sum+=(B->comp[i][k])*(C->comp[j][k]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

/* A=B^{dag}*C^{dag} */
void times_dag12_G2(G2 *__restrict__ A, G2 const *__restrict__ const B, G2 const *__restrict__ const C)
  {
  int i, j, k;
  double sum;
  
  for(i=0; i<7; i++)
     {
     for(j=0; j<7; j++)
        {
        sum=0.0;
        for(k=0; k<7; k++)
           {
           sum+=(B->comp[k][j])*(C->comp[j][k]);
           }
        A->comp[i][j]=sum;
        }
     } 
  }

                               /* ROTATION MATRICES */

/* D1 */
void D1(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c;

  s=sin(x);
  c=cos(x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[0][0]=1.0;
  A->comp[1][1]=1.0;
  A->comp[2][2]=1.0;

  A->comp[3][3]=c;
  A->comp[3][6]=-s;

  A->comp[4][4]=c;
  A->comp[4][5]=-s;

  A->comp[5][4]=s;
  A->comp[5][5]=c;

  A->comp[6][3]=s;
  A->comp[6][6]=c;
  }

/* D2 */
void D2(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c;

  s=sin(x);
  c=cos(x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[0][0]=1.0;
  A->comp[1][1]=1.0;
  A->comp[2][2]=1.0;

  A->comp[3][3]=c;
  A->comp[3][5]=s;

  A->comp[4][4]=c;
  A->comp[4][6]=-s;

  A->comp[5][3]=-s;
  A->comp[5][5]=c;

  A->comp[6][4]=s;
  A->comp[6][6]=c;
  }

/* D3 */
void D3(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c;

  s=sin(x);
  c=cos(x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[0][0]=1.0;
  A->comp[1][1]=1.0;
  A->comp[2][2]=1.0;

  A->comp[3][3]=c;
  A->comp[3][4]=-s;

  A->comp[4][3]=s;
  A->comp[4][4]=c;

  A->comp[5][5]=c;
  A->comp[5][6]=-s;

  A->comp[6][5]=s;
  A->comp[6][6]=c;
  }

/* D4 */
void D4(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c;

  s=sin(x);
  c=cos(x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[0][0]=1.0;
  A->comp[3][3]=1.0;
  A->comp[4][4]=1.0;

  A->comp[1][1]=c;
  A->comp[1][6]=s;

  A->comp[2][2]=c;
  A->comp[2][5]=s;

  A->comp[5][2]=-s;
  A->comp[5][5]=c;

  A->comp[6][1]=-s;
  A->comp[6][6]=c;
  }

/* D5 */
void D5(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c;

  s=sin(x);
  c=cos(x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[0][0]=1.0;
  A->comp[3][3]=1.0;
  A->comp[4][4]=1.0;

  A->comp[1][1]=c;
  A->comp[1][5]=-s;

  A->comp[2][2]=c;
  A->comp[2][6]=s;

  A->comp[5][1]=s;
  A->comp[5][5]=c;

  A->comp[6][2]=-s;
  A->comp[6][6]=c;
  }

/* D6 */
void D6(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c;

  s=sin(x);
  c=cos(x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[0][0]=1.0;
  A->comp[5][5]=1.0;
  A->comp[6][6]=1.0;

  A->comp[1][1]=c;
  A->comp[1][4]=s;

  A->comp[2][2]=c;
  A->comp[2][3]=-s;

  A->comp[3][2]=s;
  A->comp[3][3]=c;

  A->comp[4][1]=-s;
  A->comp[4][4]=c;
  }

/* D7 */
void D7(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c;

  s=sin(x);
  c=cos(x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[0][0]=1.0;
  A->comp[5][5]=1.0;
  A->comp[6][6]=1.0;

  A->comp[1][1]=c;
  A->comp[1][3]=-s;

  A->comp[2][2]=c;
  A->comp[2][4]=-s;

  A->comp[3][1]=s;
  A->comp[3][3]=c;

  A->comp[4][2]=s;
  A->comp[4][4]=c;
  }

/* D8 */
void D8(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c, s2, c2;

  s=sin(x);
  c=cos(x);

  s2=sin(2.0*x);
  c2=cos(2.0*x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[0][0]=1.0;

  A->comp[1][1]=c2;
  A->comp[1][2]=-s2;

  A->comp[2][1]=s2;
  A->comp[2][2]=c2;

  A->comp[3][3]=c;
  A->comp[3][4]=s;

  A->comp[4][3]=-s;
  A->comp[4][4]=c;

  A->comp[5][5]=c;
  A->comp[5][6]=-s;

  A->comp[6][5]=s;
  A->comp[6][6]=c;
  }

/* D9 */
void D9(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c, s2, c2;

  s=sin(x);
  c=cos(x);

  s2=sin(2.0*x);
  c2=cos(2.0*x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[2][2]=1.0;

  A->comp[0][0]=c2;
  A->comp[0][1]=-s2;

  A->comp[1][0]=s2;
  A->comp[1][1]=c2;

  A->comp[3][3]=c;
  A->comp[3][6]=s;

  A->comp[4][4]=c;
  A->comp[4][5]=-s;

  A->comp[5][4]=s;
  A->comp[5][5]=c;

  A->comp[6][3]=-s;
  A->comp[6][6]=c;
  }

/* D10 */
void D10(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c, s2, c2;

  s=sin(x);
  c=cos(x);

  s2=sin(2.0*x);
  c2=cos(2.0*x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[1][1]=1.0;

  A->comp[0][0]=c2;
  A->comp[0][2]=-s2;

  A->comp[2][0]=s2;
  A->comp[2][2]=c2;

  A->comp[3][3]=c;
  A->comp[3][5]=-s;

  A->comp[4][4]=c;
  A->comp[4][6]=-s;

  A->comp[5][3]=s;
  A->comp[5][5]=c;

  A->comp[6][4]=s;
  A->comp[6][6]=c;
  }

/* D11 */
void D11(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c, s2, c2;

  s=sin(x);
  c=cos(x);

  s2=sin(2.0*x);
  c2=cos(2.0*x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[4][4]=1.0;

  A->comp[0][0]=c2;
  A->comp[0][3]=-s2;

  A->comp[3][0]=s2;
  A->comp[3][3]=c2;

  A->comp[1][1]=c;
  A->comp[1][6]=-s;

  A->comp[2][2]=c;
  A->comp[2][5]=s;

  A->comp[5][2]=-s;
  A->comp[5][5]=c;

  A->comp[6][1]=s;
  A->comp[6][6]=c;
  }

/* D12 */
void D12(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c, s2, c2;

  s=sin(x);
  c=cos(x);

  s2=sin(2.0*x);
  c2=cos(2.0*x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[3][3]=1.0;

  A->comp[0][0]=c2;
  A->comp[0][4]=-s2;

  A->comp[4][0]=s2;
  A->comp[4][4]=c2;

  A->comp[1][1]=c;
  A->comp[1][5]=s;

  A->comp[2][2]=c;
  A->comp[2][6]=s;

  A->comp[5][1]=-s;
  A->comp[5][5]=c;

  A->comp[6][2]=-s;
  A->comp[6][6]=c;
  }

/* D13 */
void D13(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c, s2, c2;

  s=sin(x);
  c=cos(x);

  s2=sin(2.0*x);
  c2=cos(2.0*x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[6][6]=1.0;

  A->comp[0][0]=c2;
  A->comp[0][5]=-s2;

  A->comp[5][0]=s2;
  A->comp[5][5]=c2;

  A->comp[1][1]=c;
  A->comp[1][4]=-s;

  A->comp[2][2]=c;
  A->comp[2][3]=-s;

  A->comp[3][2]=s;
  A->comp[3][3]=c;

  A->comp[4][1]=s;
  A->comp[4][4]=c;
  }

/* D14 */
void D14(G2 *__restrict__ A, double x)
  {
  int i, j;
  double s, c, s2, c2;

  s=sin(x);
  c=cos(x);

  s2=sin(2.0*x);
  c2=cos(2.0*x);

  for(i=0;i<7;i++)
     { 
     for(j=0;j<7; j++)
        {
        A->comp[i][j]=0.0;
        }
     }
  A->comp[5][5]=1.0;

  A->comp[0][0]=c2;
  A->comp[0][6]=-s2;

  A->comp[6][0]=s2;
  A->comp[6][6]=c2;

  A->comp[1][1]=c;
  A->comp[1][3]=s;

  A->comp[2][2]=c;
  A->comp[2][4]=-s;

  A->comp[3][1]=-s;
  A->comp[3][3]=c;

  A->comp[4][2]=s;
  A->comp[4][4]=c;
  }
                                        /* END OF ROTATION MATRICES */

/* random matrix */
void rand_matrix_G2(G2 *__restrict__ A)
  {
  double aux, alpha8;
  G2 M, N;

  /* D8(alha1) */
  aux=PI2*casuale();
  D8(&M, aux);
 
  /* D9(alpha2) */
  aux=HALF_PI*casuale();
  D9(&N, aux);
  times_equal_G2(&M, &N);

  /* D8(alpha3) */
  aux=PI*casuale();
  D8(&N, aux);
  times_equal_G2(&M, &N);

  /* D3(alhpa4) */
  aux=PI2*casuale();
  D3(&N, aux);
  times_equal_G2(&M, &N);
   
  /* D2(alpha5) */
  aux=HALF_PI*casuale();
  D2(&N, aux);
  times_equal_G2(&M, &N);

  /* D3(alpha6) */
  aux=PI*casuale();
  D3(&N, aux);
  times_equal_G2(&M, &N);

  alpha8=HALF_PI*casuale();

  /* D11(alpha7) */
  aux=2.0/3.0*alpha8*casuale();
  D11(&N, aux);
  times_equal_G2(&M, &N);

  /* D5(alpha8) */
  D5(&N, alpha8);
  times_equal_G2(&M, &N);

  /* D3(alpha9) */
  aux=PI2*casuale();
  D3(&N, aux);
  times_equal_G2(&M, &N);

  /* D2(alpha10) */
  aux=HALF_PI*casuale();
  D10(&N, aux);
  times_equal_G2(&M, &N);

  /* D3(alpha11) */
  aux=PI*casuale();
  D3(&N, aux);
  times_equal_G2(&M, &N);

  /* D8(alpha12) */
  aux=PI2*casuale();
  D8(&N, aux);
  times_equal_G2(&M, &N);

  /* D9(alpha13) */
  aux=HALF_PI*casuale();
  D9(&N, aux);
  times_equal_G2(&M, &N);

  /* D8(alpha14) */
  aux=PI*casuale();
  D8(&N, aux);
  times_equal_G2(&M, &N);

  equal_G2(A, &M);
  }

/* 1/7 of the real part of the trace */
double retr_G2(G2 const *__restrict__ const A)
   {
   return (A->comp[0][0]+A->comp[1][1]+A->comp[2][2]+A->comp[3][3]+A->comp[4][4]+A->comp[5][5]+A->comp[6][6])/7.0;
   }

/* 1/7 of the imaginary part of the trace */
double imtr_G2(G2 const *__restrict__ const A)
   {
   return 0.0;
   }

/* L^2 norm */
double norm_G2(G2 const *__restrict__ const A)
   {
   int i, j;
   double aux=0.0;
  
   for(i=0; i<7; i++)
      {
      for(j=0; j<7; j++)
         {
         aux+=A->comp[i][j]*A->comp[i][j];
         }
      }

   return sqrt(aux);
   }


void print_on_screen_G2(G2 const *__restrict__ const A) 
   {
   int i;
   
   for(i=0; i<7; i++)
      {
      printf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", A->comp[i][0],
                                                            A->comp[i][1],
                                                            A->comp[i][2],
                                                            A->comp[i][3],
                                                            A->comp[i][4],
                                                            A->comp[i][5],
                                                            A->comp[i][6]);
      }
   }

void print_on_file_G2(FILE *fp, G2 const *__restrict__ const A)
   {
   int i;
   
   for(i=0; i<7; i++)
      {
      fprintf(fp, "%.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", A->comp[i][0],
                                                                 A->comp[i][1],
                                                                 A->comp[i][2],
                                                                 A->comp[i][3],
                                                                 A->comp[i][4],
                                                                 A->comp[i][5],
                                                                 A->comp[i][6]);
      }
   }

void read_from_file_G2(FILE *fp, G2 *__restrict__ A)
   {
   int i, err;
   
   for(i=0; i<7; i++)
      {
      err=fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg\n", &(A->comp[i][0]),
                                                      &(A->comp[i][1]),
                                                      &(A->comp[i][2]),
                                                      &(A->comp[i][3]),
                                                      &(A->comp[i][4]),
                                                      &(A->comp[i][5]),
                                                      &(A->comp[i][6]));
      if(err!=7)
        {
        fprintf(stderr, "Problems reading G2 matrix from file");
        }
      }
   }



#endif
