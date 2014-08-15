#ifndef SU2_C
#define SU2_C

#include<math.h>
#include<stdio.h>

#include"../Rng/random.h"

#include"su2.h"

void init_Su2(Su2 *__restrict__ A, double vec[4])
  {
  A->comp[0]=vec[0];
  A->comp[1]=vec[1];
  A->comp[2]=vec[2];
  A->comp[3]=vec[3];
  }
  

/* A=1 */
void one_Su2(Su2 *__restrict__ A)
  {
  A->comp[0]=1.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


/* A=0 */
void zero_Su2(Su2 *__restrict__ A)
  {
  A->comp[0]=0.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


/* A=B */
void equal_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B)
  {
  A->comp[0]=B->comp[0];
  A->comp[1]=B->comp[1];
  A->comp[2]=B->comp[2];
  A->comp[3]=B->comp[3];
  }


/* A=B^{dag} */
void equal_dag_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B)
  {
  A->comp[0]=B->comp[0];
  A->comp[1]=-B->comp[1];
  A->comp[2]=-B->comp[2];
  A->comp[3]=-B->comp[3];
  }


/* A+=B */
void plus_equal_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B)
  {
  A->comp[0]+=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


/* A+=B^{dag} */
void plus_equal_dag_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B)
  {
  A->comp[0]+=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


/* A-=B */
void minus_equal_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B)
  {
  A->comp[0]-=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


/* A-=B^{dag} */
void minus_equal_dag_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B)
  {
  A->comp[0]-=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


/* A=b*B+c*C */
void lin_comb_Su2(Su2 *__restrict__ A, double b, Su2 const *__restrict__ const B, double c, Su2 const *__restrict__ const C)
  {
  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] + c*C->comp[1];
  A->comp[2]= b*B->comp[2] + c*C->comp[2];
  A->comp[3]= b*B->comp[3] + c*C->comp[3];
  }


/* A=b*B^{dag}+c*C */
void lin_comb_dag1_Su2(Su2 *__restrict__ A, double b, Su2 const *__restrict__ const B, double c, Su2 const *__restrict__ const C)
  {
  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] + c*C->comp[1];
  A->comp[2]= -b*B->comp[2] + c*C->comp[2];
  A->comp[3]= -b*B->comp[3] + c*C->comp[3];
  }


/* A=b*B+c*C^{dag} */
void lin_comb_dag2_Su2(Su2 *__restrict__ A, double b, Su2 const *__restrict__ const B, double c, Su2 const *__restrict__ const C)
  {
  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] - c*C->comp[1];
  A->comp[2]= b*B->comp[2] - c*C->comp[2];
  A->comp[3]= b*B->comp[3] - c*C->comp[3];
  }


/* A=b*B^{dag}+c*C^{dag} */
void lin_comb_dag12_Su2(Su2 *__restrict__ A, double b, Su2 const *__restrict__ const B, double c, Su2 const *__restrict__ const C)
  {
  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] - c*C->comp[1];
  A->comp[2]= -b*B->comp[2] - c*C->comp[2];
  A->comp[3]= -b*B->comp[3] - c*C->comp[3];
  }


/* A*=r */
void times_equal_real_Su2(Su2 *__restrict__ A, double r)
  {
  A->comp[0]*=r;
  A->comp[1]*=r;
  A->comp[2]*=r;
  A->comp[3]*=r;
  }


/* A*=B       */
void times_equal_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B)
  {
  double a0, a1, a2, a3;

  a0=A->comp[0];
  a1=A->comp[1];
  a2=A->comp[2];
  a3=A->comp[3];
 
  A->comp[0]= a0*B->comp[0] - a1*B->comp[1] - a2*B->comp[2] - a3*B->comp[3];
  A->comp[1]= a0*B->comp[1] + a1*B->comp[0] - a2*B->comp[3] + a3*B->comp[2];
  A->comp[2]= a0*B->comp[2] + a2*B->comp[0] + a1*B->comp[3] - a3*B->comp[1];
  A->comp[3]= a0*B->comp[3] + a3*B->comp[0] - a1*B->comp[2] + a2*B->comp[1];
  }


/* A*=B^{dag} */
void times_equal_dag_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B)
  {
  double a0, a1, a2, a3;

  a0=A->comp[0];
  a1=A->comp[1];
  a2=A->comp[2];
  a3=A->comp[3];
 
  A->comp[0]=  a0*B->comp[0] + a1*B->comp[1] + a2*B->comp[2] + a3*B->comp[3];
  A->comp[1]= -a0*B->comp[1] + a1*B->comp[0] + a2*B->comp[3] - a3*B->comp[2];
  A->comp[2]= -a0*B->comp[2] + a2*B->comp[0] - a1*B->comp[3] + a3*B->comp[1];
  A->comp[3]= -a0*B->comp[3] + a3*B->comp[0] + a1*B->comp[2] - a2*B->comp[1];
  }


/* A=B*C             */
void times_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B, Su2 const *__restrict__ const C)
  {
  A->comp[0]= B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


/* A=B^{dag}*C       */
void times_dag1_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B, Su2 const *__restrict__ const C)
  {
  A->comp[0]= B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


/* A=B*C^{dag}       */
void times_dag2_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B, Su2 const *__restrict__ const C)
  {
  A->comp[0]=  B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


/* A=B^{dag}*C^{dag}       */
void times_dag12_Su2(Su2 *__restrict__ A, Su2 const *__restrict__ const B, Su2 const *__restrict__ const C)
  {
  A->comp[0]=  B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


/* random SU(2) matrix */
void rand_matrix_Su2(Su2 *__restrict__ A)
  {
  double p0, p1, p2, p3, p;

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

  A->comp[0]=p0;
  A->comp[1]=p1;
  A->comp[2]=p2;
  A->comp[3]=p3;
  }


/* random SU(2) matrix with p0 given (used in the update) */
void rand_matrix_p0_Su2(double p0, Su2 *__restrict__ A)
  {
  register double p1, p2, p3, p;

  p=2.0;
  while(p>1.0)
       { 
       p1=1.0-2.0*casuale();
       p2=1.0-2.0*casuale();
       p3=1.0-2.0*casuale();
       p=p1*p1+p2*p2+p3*p3;
       }

  p/=(1.0-p0*p0);
  p=sqrt(p);

  p1/=p;
  p2/=p;
  p3/=p;

  A->comp[0]=p0;
  A->comp[1]=p1;
  A->comp[2]=p2;
  A->comp[3]=p3;
  }


/* real part of the trace /2 */
double retr_Su2(Su2 const *__restrict__ const A)
   {
   return A->comp[0];
   }


/* imaginary part of the trace /2 */
double imtr_Su2(Su2 const *__restrict__ const A)
   {
   return 0.0;
   }

/* sqrt of the determinant */
double sqrtdet_Su2(Su2 const *__restrict__ const A)
   {
   return sqrt(A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1] + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3]);
   }


/* unitarize the matrix */
void unitarize_Su2(Su2 *__restrict__ A)
   {
   double p;

   p=A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1] + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3];
   p=1.0/sqrt(p);

   A->comp[0]*=p;
   A->comp[1]*=p;
   A->comp[2]*=p;
   A->comp[3]*=p;
   }


/* print on screen */
void print_on_screen_Su2(Su2 const *__restrict__ const A)
  {
  printf("%.16g %.16g %.16g %.16g\n", A->comp[0], A->comp[1], A->comp[2], A->comp[3]);
  }


/* print on file */
void print_on_file_Su2(FILE *fp, Su2 const *__restrict__ const A)
  {
  fprintf(fp, "%.16g %.16g %.16g %.16g\n", A->comp[0], A->comp[1], A->comp[2], A->comp[3]);
  }


/* read form file */
void read_from_file_Su2(FILE *fp, Su2 *__restrict__ A)
  {
  int err=fscanf(fp, "%lg %lg %lg %lg", &(A->comp[0]), &(A->comp[1]), &(A->comp[2]), &(A->comp[3]));

  if(err!=4)
    {
    fprintf(stderr, "Problems reading Su2 matrix from file");
    }
  }


#endif
