#ifndef SUN_AUX_C
#define SUN_AUX_C

#include"../Macro/macro.h"

#include<complex.h>
#include<math.h>

#include"sun.h"
#include"../Su2/su2.h"

/* given the matrix N*N "in" extracts the i, j lines and column and
   gives "xi" real number and "u" in SU(2) 
   4 xi^2 = redet2[s-s^(dag)+1*tr(s^(dag))]
   u = [s-s^(dag)+1*tr(s^(dag))]/2/xi 
   (see Kennedy, Pendleton Phys. Lett. B 156, 393 (1985))  */
void ennetodue(SuN const *__restrict__ const in, 
               int i,  
               int j, 
               double *xi, 
               Su2 *__restrict__ u)
   {
   double s[2][2][2], auxr[2][2], auxi[2][2];
   double p;

   s[0][0][0]=creal(in->comp[i][i]);
   s[0][0][1]=cimag(in->comp[i][i]);

   s[0][1][0]=creal(in->comp[i][j]);
   s[0][1][1]=cimag(in->comp[i][j]);
   
   s[1][0][0]=creal(in->comp[j][i]);
   s[1][0][1]=cimag(in->comp[j][i]);

   s[1][1][0]=creal(in->comp[j][j]);
   s[1][1][1]=cimag(in->comp[j][j]);

   auxr[0][0]=s[0][0][0]+s[1][1][0];
   auxi[0][0]=s[0][0][1]-s[1][1][1];

   auxr[0][1]=s[0][1][0]-s[1][0][0];
   auxi[0][1]=s[0][1][1]+s[1][0][1];

   auxr[1][0]=s[1][0][0]-s[0][1][0];
   auxi[1][0]=s[1][0][1]+s[0][1][1];

   auxr[1][1]=s[0][0][0]+s[1][1][0];
   auxi[1][1]=s[1][1][1]+s[1][1][1]-s[0][0][1]-s[1][1][1];

   p=auxr[0][0]*auxr[1][1]-auxi[0][0]*auxi[1][1]-auxr[0][1]*auxr[1][0]+auxi[0][1]*auxi[1][0];
   p=sqrt(p);
  
   (*xi)=p/2.0;

   if(*xi>MIN_VALUE)
     {
     auxr[0][0]/=p;
     auxi[0][1]/=p;
     auxr[0][1]/=p;
     auxi[0][0]/=p;
     }
   u->comp[0]=auxr[0][0];
   u->comp[1]=auxi[0][1];
   u->comp[2]=auxr[0][1];
   u->comp[3]=auxi[0][0];
   }


/* given a 2*2 matrix extend to N*N with 1 on the diagonal */
void duetoenne(Su2 const *__restrict__ const in, 
               int i, 
               int j, 
               SuN *__restrict__ out)
   {
   one_SuN(out);

   out->comp[i][i]= in->comp[0] + (in->comp[3])*I;
   out->comp[i][j]= in->comp[2] + (in->comp[1])*I;
   out->comp[j][i]=-in->comp[2] + (in->comp[1])*I;
   out->comp[j][j]= in->comp[0] - (in->comp[3])*I;
   }

#endif
