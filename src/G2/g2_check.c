#ifndef G2_CHECK_C
#define G2_CHECK_C

#include"../Macro/macro.h"

#include<math.h>

#include"g2.h"
#include"../GParam/gparam.h"

/* verify the cubic condition needed for an SO(7) matrix to be in G2
   return 1 if ok, 0 else */
int cc_check_G2(G2 const * const A, double *cc)
  {
  int i;
  double  test[35];
  double test_result;

  /*1,2,3 */
  test[0] = -1.0  -(A->comp[0][2]*A->comp[1][1]*A->comp[2][0]) + A->comp[0][1]*A->comp[1][2]*A->comp[2][0] + A->comp[0][2]*A->comp[1][0]*A->comp[2][1] - 
            A->comp[0][0]*A->comp[1][2]*A->comp[2][1] - A->comp[0][1]*A->comp[1][0]*A->comp[2][2] + A->comp[0][0]*A->comp[1][1]*A->comp[2][2] - 
            A->comp[0][2]*A->comp[3][1]*A->comp[4][0] + A->comp[0][1]*A->comp[3][2]*A->comp[4][0] + A->comp[0][2]*A->comp[3][0]*A->comp[4][1] - 
            A->comp[0][0]*A->comp[3][2]*A->comp[4][1] - A->comp[0][1]*A->comp[3][0]*A->comp[4][2] + A->comp[0][0]*A->comp[3][1]*A->comp[4][2] - 
            A->comp[1][2]*A->comp[3][1]*A->comp[5][0] + A->comp[1][1]*A->comp[3][2]*A->comp[5][0] + A->comp[2][2]*A->comp[4][1]*A->comp[5][0] - 
            A->comp[2][1]*A->comp[4][2]*A->comp[5][0] + A->comp[1][2]*A->comp[3][0]*A->comp[5][1] - A->comp[1][0]*A->comp[3][2]*A->comp[5][1] - 
            A->comp[2][2]*A->comp[4][0]*A->comp[5][1] + A->comp[2][0]*A->comp[4][2]*A->comp[5][1] - A->comp[1][1]*A->comp[3][0]*A->comp[5][2] +
            A->comp[1][0]*A->comp[3][1]*A->comp[5][2] + A->comp[2][1]*A->comp[4][0]*A->comp[5][2] - A->comp[2][0]*A->comp[4][1]*A->comp[5][2] - 
            A->comp[2][2]*A->comp[3][1]*A->comp[6][0] + A->comp[2][1]*A->comp[3][2]*A->comp[6][0] - A->comp[1][2]*A->comp[4][1]*A->comp[6][0] + 
            A->comp[1][1]*A->comp[4][2]*A->comp[6][0] + A->comp[0][2]*A->comp[5][1]*A->comp[6][0] - A->comp[0][1]*A->comp[5][2]*A->comp[6][0] + 
            A->comp[2][2]*A->comp[3][0]*A->comp[6][1] - A->comp[2][0]*A->comp[3][2]*A->comp[6][1] + A->comp[1][2]*A->comp[4][0]*A->comp[6][1] - 
            A->comp[1][0]*A->comp[4][2]*A->comp[6][1] - A->comp[0][2]*A->comp[5][0]*A->comp[6][1] + A->comp[0][0]*A->comp[5][2]*A->comp[6][1] - 
            A->comp[2][1]*A->comp[3][0]*A->comp[6][2] + A->comp[2][0]*A->comp[3][1]*A->comp[6][2] - A->comp[1][1]*A->comp[4][0]*A->comp[6][2] + 
            A->comp[1][0]*A->comp[4][1]*A->comp[6][2] + A->comp[0][1]*A->comp[5][0]*A->comp[6][2] - A->comp[0][0]*A->comp[5][1]*A->comp[6][2];

  /* 1,2,4 */
  test[1] = -(A->comp[0][3]*A->comp[1][1]*A->comp[2][0]) + A->comp[0][1]*A->comp[1][3]*A->comp[2][0] + A->comp[0][3]*A->comp[1][0]*A->comp[2][1] - 
            A->comp[0][0]*A->comp[1][3]*A->comp[2][1] - A->comp[0][1]*A->comp[1][0]*A->comp[2][3] + A->comp[0][0]*A->comp[1][1]*A->comp[2][3] - 
            A->comp[0][3]*A->comp[3][1]*A->comp[4][0] + A->comp[0][1]*A->comp[3][3]*A->comp[4][0] + A->comp[0][3]*A->comp[3][0]*A->comp[4][1] - 
            A->comp[0][0]*A->comp[3][3]*A->comp[4][1] - A->comp[0][1]*A->comp[3][0]*A->comp[4][3] + A->comp[0][0]*A->comp[3][1]*A->comp[4][3] - 
            A->comp[1][3]*A->comp[3][1]*A->comp[5][0] + A->comp[1][1]*A->comp[3][3]*A->comp[5][0] + A->comp[2][3]*A->comp[4][1]*A->comp[5][0] - 
            A->comp[2][1]*A->comp[4][3]*A->comp[5][0] + A->comp[1][3]*A->comp[3][0]*A->comp[5][1] - A->comp[1][0]*A->comp[3][3]*A->comp[5][1] - 
            A->comp[2][3]*A->comp[4][0]*A->comp[5][1] + A->comp[2][0]*A->comp[4][3]*A->comp[5][1] - A->comp[1][1]*A->comp[3][0]*A->comp[5][3] +
            A->comp[1][0]*A->comp[3][1]*A->comp[5][3] + A->comp[2][1]*A->comp[4][0]*A->comp[5][3] - A->comp[2][0]*A->comp[4][1]*A->comp[5][3] - 
            A->comp[2][3]*A->comp[3][1]*A->comp[6][0] + A->comp[2][1]*A->comp[3][3]*A->comp[6][0] - A->comp[1][3]*A->comp[4][1]*A->comp[6][0] + 
            A->comp[1][1]*A->comp[4][3]*A->comp[6][0] + A->comp[0][3]*A->comp[5][1]*A->comp[6][0] - A->comp[0][1]*A->comp[5][3]*A->comp[6][0] + 
            A->comp[2][3]*A->comp[3][0]*A->comp[6][1] - A->comp[2][0]*A->comp[3][3]*A->comp[6][1] + A->comp[1][3]*A->comp[4][0]*A->comp[6][1] - 
            A->comp[1][0]*A->comp[4][3]*A->comp[6][1] - A->comp[0][3]*A->comp[5][0]*A->comp[6][1] + A->comp[0][0]*A->comp[5][3]*A->comp[6][1] - 
            A->comp[2][1]*A->comp[3][0]*A->comp[6][3] + A->comp[2][0]*A->comp[3][1]*A->comp[6][3] - A->comp[1][1]*A->comp[4][0]*A->comp[6][3] + 
            A->comp[1][0]*A->comp[4][1]*A->comp[6][3] + A->comp[0][1]*A->comp[5][0]*A->comp[6][3] - A->comp[0][0]*A->comp[5][1]*A->comp[6][3];

  /* 1,2,5 */
  test[2] = -(A->comp[0][4]*A->comp[1][1]*A->comp[2][0]) + A->comp[0][1]*A->comp[1][4]*A->comp[2][0] + A->comp[0][4]*A->comp[1][0]*A->comp[2][1] - 
           A->comp[0][0]*A->comp[1][4]*A->comp[2][1] - A->comp[0][1]*A->comp[1][0]*A->comp[2][4] + A->comp[0][0]*A->comp[1][1]*A->comp[2][4] - 
           A->comp[0][4]*A->comp[3][1]*A->comp[4][0] + A->comp[0][1]*A->comp[3][4]*A->comp[4][0] + A->comp[0][4]*A->comp[3][0]*A->comp[4][1] - 
           A->comp[0][0]*A->comp[3][4]*A->comp[4][1] - A->comp[0][1]*A->comp[3][0]*A->comp[4][4] + A->comp[0][0]*A->comp[3][1]*A->comp[4][4] - 
           A->comp[1][4]*A->comp[3][1]*A->comp[5][0] + A->comp[1][1]*A->comp[3][4]*A->comp[5][0] + A->comp[2][4]*A->comp[4][1]*A->comp[5][0] - 
           A->comp[2][1]*A->comp[4][4]*A->comp[5][0] + A->comp[1][4]*A->comp[3][0]*A->comp[5][1] - A->comp[1][0]*A->comp[3][4]*A->comp[5][1] - 
           A->comp[2][4]*A->comp[4][0]*A->comp[5][1] + A->comp[2][0]*A->comp[4][4]*A->comp[5][1] - A->comp[1][1]*A->comp[3][0]*A->comp[5][4] +
           A->comp[1][0]*A->comp[3][1]*A->comp[5][4] + A->comp[2][1]*A->comp[4][0]*A->comp[5][4] - A->comp[2][0]*A->comp[4][1]*A->comp[5][4] - 
           A->comp[2][4]*A->comp[3][1]*A->comp[6][0] + A->comp[2][1]*A->comp[3][4]*A->comp[6][0] - A->comp[1][4]*A->comp[4][1]*A->comp[6][0] + 
           A->comp[1][1]*A->comp[4][4]*A->comp[6][0] + A->comp[0][4]*A->comp[5][1]*A->comp[6][0] - A->comp[0][1]*A->comp[5][4]*A->comp[6][0] + 
           A->comp[2][4]*A->comp[3][0]*A->comp[6][1] - A->comp[2][0]*A->comp[3][4]*A->comp[6][1] + A->comp[1][4]*A->comp[4][0]*A->comp[6][1] - 
           A->comp[1][0]*A->comp[4][4]*A->comp[6][1] - A->comp[0][4]*A->comp[5][0]*A->comp[6][1] + A->comp[0][0]*A->comp[5][4]*A->comp[6][1] - 
           A->comp[2][1]*A->comp[3][0]*A->comp[6][4] + A->comp[2][0]*A->comp[3][1]*A->comp[6][4] - A->comp[1][1]*A->comp[4][0]*A->comp[6][4] + 
           A->comp[1][0]*A->comp[4][1]*A->comp[6][4] + A->comp[0][1]*A->comp[5][0]*A->comp[6][4] - A->comp[0][0]*A->comp[5][1]*A->comp[6][4];

  /* 1,2,6 */
  test[3] = -(A->comp[0][5]*A->comp[1][1]*A->comp[2][0]) + A->comp[0][1]*A->comp[1][5]*A->comp[2][0] + A->comp[0][5]*A->comp[1][0]*A->comp[2][1] - 
           A->comp[0][0]*A->comp[1][5]*A->comp[2][1] - A->comp[0][1]*A->comp[1][0]*A->comp[2][5] + A->comp[0][0]*A->comp[1][1]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][1]*A->comp[4][0] + A->comp[0][1]*A->comp[3][5]*A->comp[4][0] + A->comp[0][5]*A->comp[3][0]*A->comp[4][1] - 
           A->comp[0][0]*A->comp[3][5]*A->comp[4][1] - A->comp[0][1]*A->comp[3][0]*A->comp[4][5] + A->comp[0][0]*A->comp[3][1]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][1]*A->comp[5][0] + A->comp[1][1]*A->comp[3][5]*A->comp[5][0] + A->comp[2][5]*A->comp[4][1]*A->comp[5][0] - 
           A->comp[2][1]*A->comp[4][5]*A->comp[5][0] + A->comp[1][5]*A->comp[3][0]*A->comp[5][1] - A->comp[1][0]*A->comp[3][5]*A->comp[5][1] - 
           A->comp[2][5]*A->comp[4][0]*A->comp[5][1] + A->comp[2][0]*A->comp[4][5]*A->comp[5][1] - A->comp[1][1]*A->comp[3][0]*A->comp[5][5] + 
           A->comp[1][0]*A->comp[3][1]*A->comp[5][5] + A->comp[2][1]*A->comp[4][0]*A->comp[5][5] - A->comp[2][0]*A->comp[4][1]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][1]*A->comp[6][0] + A->comp[2][1]*A->comp[3][5]*A->comp[6][0] - A->comp[1][5]*A->comp[4][1]*A->comp[6][0] + 
           A->comp[1][1]*A->comp[4][5]*A->comp[6][0] + A->comp[0][5]*A->comp[5][1]*A->comp[6][0] - A->comp[0][1]*A->comp[5][5]*A->comp[6][0] + 
           A->comp[2][5]*A->comp[3][0]*A->comp[6][1] - A->comp[2][0]*A->comp[3][5]*A->comp[6][1] + A->comp[1][5]*A->comp[4][0]*A->comp[6][1] - 
           A->comp[1][0]*A->comp[4][5]*A->comp[6][1] - A->comp[0][5]*A->comp[5][0]*A->comp[6][1] + A->comp[0][0]*A->comp[5][5]*A->comp[6][1] - 
           A->comp[2][1]*A->comp[3][0]*A->comp[6][5] + A->comp[2][0]*A->comp[3][1]*A->comp[6][5] - A->comp[1][1]*A->comp[4][0]*A->comp[6][5] + 
           A->comp[1][0]*A->comp[4][1]*A->comp[6][5] + A->comp[0][1]*A->comp[5][0]*A->comp[6][5] - A->comp[0][0]*A->comp[5][1]*A->comp[6][5];

  /* 1,2,7 */
  test[4] = -(A->comp[0][6]*A->comp[1][1]*A->comp[2][0]) + A->comp[0][1]*A->comp[1][6]*A->comp[2][0] + A->comp[0][6]*A->comp[1][0]*A->comp[2][1] - 
           A->comp[0][0]*A->comp[1][6]*A->comp[2][1] - A->comp[0][1]*A->comp[1][0]*A->comp[2][6] + A->comp[0][0]*A->comp[1][1]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][1]*A->comp[4][0] + A->comp[0][1]*A->comp[3][6]*A->comp[4][0] + A->comp[0][6]*A->comp[3][0]*A->comp[4][1] - 
           A->comp[0][0]*A->comp[3][6]*A->comp[4][1] - A->comp[0][1]*A->comp[3][0]*A->comp[4][6] + A->comp[0][0]*A->comp[3][1]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][1]*A->comp[5][0] + A->comp[1][1]*A->comp[3][6]*A->comp[5][0] + A->comp[2][6]*A->comp[4][1]*A->comp[5][0] - 
           A->comp[2][1]*A->comp[4][6]*A->comp[5][0] + A->comp[1][6]*A->comp[3][0]*A->comp[5][1] - A->comp[1][0]*A->comp[3][6]*A->comp[5][1] - 
           A->comp[2][6]*A->comp[4][0]*A->comp[5][1] + A->comp[2][0]*A->comp[4][6]*A->comp[5][1] - A->comp[1][1]*A->comp[3][0]*A->comp[5][6] +
           A->comp[1][0]*A->comp[3][1]*A->comp[5][6] + A->comp[2][1]*A->comp[4][0]*A->comp[5][6] - A->comp[2][0]*A->comp[4][1]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][1]*A->comp[6][0] + A->comp[2][1]*A->comp[3][6]*A->comp[6][0] - A->comp[1][6]*A->comp[4][1]*A->comp[6][0] + 
           A->comp[1][1]*A->comp[4][6]*A->comp[6][0] + A->comp[0][6]*A->comp[5][1]*A->comp[6][0] - A->comp[0][1]*A->comp[5][6]*A->comp[6][0] + 
           A->comp[2][6]*A->comp[3][0]*A->comp[6][1] - A->comp[2][0]*A->comp[3][6]*A->comp[6][1] + A->comp[1][6]*A->comp[4][0]*A->comp[6][1] - 
           A->comp[1][0]*A->comp[4][6]*A->comp[6][1] - A->comp[0][6]*A->comp[5][0]*A->comp[6][1] + A->comp[0][0]*A->comp[5][6]*A->comp[6][1] - 
           A->comp[2][1]*A->comp[3][0]*A->comp[6][6] + A->comp[2][0]*A->comp[3][1]*A->comp[6][6] - A->comp[1][1]*A->comp[4][0]*A->comp[6][6] + 
           A->comp[1][0]*A->comp[4][1]*A->comp[6][6] + A->comp[0][1]*A->comp[5][0]*A->comp[6][6] - A->comp[0][0]*A->comp[5][1]*A->comp[6][6];
  
  /* 1,3,4 */
  test[5] = -(A->comp[0][3]*A->comp[1][2]*A->comp[2][0]) + A->comp[0][2]*A->comp[1][3]*A->comp[2][0] + A->comp[0][3]*A->comp[1][0]*A->comp[2][2] - 
           A->comp[0][0]*A->comp[1][3]*A->comp[2][2] - A->comp[0][2]*A->comp[1][0]*A->comp[2][3] + A->comp[0][0]*A->comp[1][2]*A->comp[2][3] - 
           A->comp[0][3]*A->comp[3][2]*A->comp[4][0] + A->comp[0][2]*A->comp[3][3]*A->comp[4][0] + A->comp[0][3]*A->comp[3][0]*A->comp[4][2] - 
           A->comp[0][0]*A->comp[3][3]*A->comp[4][2] - A->comp[0][2]*A->comp[3][0]*A->comp[4][3] + A->comp[0][0]*A->comp[3][2]*A->comp[4][3] - 
           A->comp[1][3]*A->comp[3][2]*A->comp[5][0] + A->comp[1][2]*A->comp[3][3]*A->comp[5][0] + A->comp[2][3]*A->comp[4][2]*A->comp[5][0] - 
           A->comp[2][2]*A->comp[4][3]*A->comp[5][0] + A->comp[1][3]*A->comp[3][0]*A->comp[5][2] - A->comp[1][0]*A->comp[3][3]*A->comp[5][2] - 
           A->comp[2][3]*A->comp[4][0]*A->comp[5][2] + A->comp[2][0]*A->comp[4][3]*A->comp[5][2] - A->comp[1][2]*A->comp[3][0]*A->comp[5][3] +
           A->comp[1][0]*A->comp[3][2]*A->comp[5][3] + A->comp[2][2]*A->comp[4][0]*A->comp[5][3] - A->comp[2][0]*A->comp[4][2]*A->comp[5][3] - 
           A->comp[2][3]*A->comp[3][2]*A->comp[6][0] + A->comp[2][2]*A->comp[3][3]*A->comp[6][0] - A->comp[1][3]*A->comp[4][2]*A->comp[6][0] + 
           A->comp[1][2]*A->comp[4][3]*A->comp[6][0] + A->comp[0][3]*A->comp[5][2]*A->comp[6][0] - A->comp[0][2]*A->comp[5][3]*A->comp[6][0] + 
           A->comp[2][3]*A->comp[3][0]*A->comp[6][2] - A->comp[2][0]*A->comp[3][3]*A->comp[6][2] + A->comp[1][3]*A->comp[4][0]*A->comp[6][2] - 
           A->comp[1][0]*A->comp[4][3]*A->comp[6][2] - A->comp[0][3]*A->comp[5][0]*A->comp[6][2] + A->comp[0][0]*A->comp[5][3]*A->comp[6][2] - 
           A->comp[2][2]*A->comp[3][0]*A->comp[6][3] + A->comp[2][0]*A->comp[3][2]*A->comp[6][3] - A->comp[1][2]*A->comp[4][0]*A->comp[6][3] + 
           A->comp[1][0]*A->comp[4][2]*A->comp[6][3] + A->comp[0][2]*A->comp[5][0]*A->comp[6][3] - A->comp[0][0]*A->comp[5][2]*A->comp[6][3];

  /* 1,3,5 */
  test[6] = -(A->comp[0][4]*A->comp[1][2]*A->comp[2][0]) + A->comp[0][2]*A->comp[1][4]*A->comp[2][0] + A->comp[0][4]*A->comp[1][0]*A->comp[2][2] - 
           A->comp[0][0]*A->comp[1][4]*A->comp[2][2] - A->comp[0][2]*A->comp[1][0]*A->comp[2][4] + A->comp[0][0]*A->comp[1][2]*A->comp[2][4] - 
           A->comp[0][4]*A->comp[3][2]*A->comp[4][0] + A->comp[0][2]*A->comp[3][4]*A->comp[4][0] + A->comp[0][4]*A->comp[3][0]*A->comp[4][2] - 
           A->comp[0][0]*A->comp[3][4]*A->comp[4][2] - A->comp[0][2]*A->comp[3][0]*A->comp[4][4] + A->comp[0][0]*A->comp[3][2]*A->comp[4][4] - 
           A->comp[1][4]*A->comp[3][2]*A->comp[5][0] + A->comp[1][2]*A->comp[3][4]*A->comp[5][0] + A->comp[2][4]*A->comp[4][2]*A->comp[5][0] - 
           A->comp[2][2]*A->comp[4][4]*A->comp[5][0] + A->comp[1][4]*A->comp[3][0]*A->comp[5][2] - A->comp[1][0]*A->comp[3][4]*A->comp[5][2] - 
           A->comp[2][4]*A->comp[4][0]*A->comp[5][2] + A->comp[2][0]*A->comp[4][4]*A->comp[5][2] - A->comp[1][2]*A->comp[3][0]*A->comp[5][4] +
           A->comp[1][0]*A->comp[3][2]*A->comp[5][4] + A->comp[2][2]*A->comp[4][0]*A->comp[5][4] - A->comp[2][0]*A->comp[4][2]*A->comp[5][4] - 
           A->comp[2][4]*A->comp[3][2]*A->comp[6][0] + A->comp[2][2]*A->comp[3][4]*A->comp[6][0] - A->comp[1][4]*A->comp[4][2]*A->comp[6][0] + 
           A->comp[1][2]*A->comp[4][4]*A->comp[6][0] + A->comp[0][4]*A->comp[5][2]*A->comp[6][0] - A->comp[0][2]*A->comp[5][4]*A->comp[6][0] + 
           A->comp[2][4]*A->comp[3][0]*A->comp[6][2] - A->comp[2][0]*A->comp[3][4]*A->comp[6][2] + A->comp[1][4]*A->comp[4][0]*A->comp[6][2] - 
           A->comp[1][0]*A->comp[4][4]*A->comp[6][2] - A->comp[0][4]*A->comp[5][0]*A->comp[6][2] + A->comp[0][0]*A->comp[5][4]*A->comp[6][2] - 
           A->comp[2][2]*A->comp[3][0]*A->comp[6][4] + A->comp[2][0]*A->comp[3][2]*A->comp[6][4] - A->comp[1][2]*A->comp[4][0]*A->comp[6][4] + 
           A->comp[1][0]*A->comp[4][2]*A->comp[6][4] + A->comp[0][2]*A->comp[5][0]*A->comp[6][4] - A->comp[0][0]*A->comp[5][2]*A->comp[6][4];

  /* 1,3,6 */
  test[7] = -(A->comp[0][5]*A->comp[1][2]*A->comp[2][0]) + A->comp[0][2]*A->comp[1][5]*A->comp[2][0] + A->comp[0][5]*A->comp[1][0]*A->comp[2][2] - 
           A->comp[0][0]*A->comp[1][5]*A->comp[2][2] - A->comp[0][2]*A->comp[1][0]*A->comp[2][5] + A->comp[0][0]*A->comp[1][2]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][2]*A->comp[4][0] + A->comp[0][2]*A->comp[3][5]*A->comp[4][0] + A->comp[0][5]*A->comp[3][0]*A->comp[4][2] - 
           A->comp[0][0]*A->comp[3][5]*A->comp[4][2] - A->comp[0][2]*A->comp[3][0]*A->comp[4][5] + A->comp[0][0]*A->comp[3][2]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][2]*A->comp[5][0] + A->comp[1][2]*A->comp[3][5]*A->comp[5][0] + A->comp[2][5]*A->comp[4][2]*A->comp[5][0] - 
           A->comp[2][2]*A->comp[4][5]*A->comp[5][0] + A->comp[1][5]*A->comp[3][0]*A->comp[5][2] - A->comp[1][0]*A->comp[3][5]*A->comp[5][2] - 
           A->comp[2][5]*A->comp[4][0]*A->comp[5][2] + A->comp[2][0]*A->comp[4][5]*A->comp[5][2] - A->comp[1][2]*A->comp[3][0]*A->comp[5][5] +
           A->comp[1][0]*A->comp[3][2]*A->comp[5][5] + A->comp[2][2]*A->comp[4][0]*A->comp[5][5] - A->comp[2][0]*A->comp[4][2]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][2]*A->comp[6][0] + A->comp[2][2]*A->comp[3][5]*A->comp[6][0] - A->comp[1][5]*A->comp[4][2]*A->comp[6][0] + 
           A->comp[1][2]*A->comp[4][5]*A->comp[6][0] + A->comp[0][5]*A->comp[5][2]*A->comp[6][0] - A->comp[0][2]*A->comp[5][5]*A->comp[6][0] + 
           A->comp[2][5]*A->comp[3][0]*A->comp[6][2] - A->comp[2][0]*A->comp[3][5]*A->comp[6][2] + A->comp[1][5]*A->comp[4][0]*A->comp[6][2] - 
           A->comp[1][0]*A->comp[4][5]*A->comp[6][2] - A->comp[0][5]*A->comp[5][0]*A->comp[6][2] + A->comp[0][0]*A->comp[5][5]*A->comp[6][2] - 
           A->comp[2][2]*A->comp[3][0]*A->comp[6][5] + A->comp[2][0]*A->comp[3][2]*A->comp[6][5] - A->comp[1][2]*A->comp[4][0]*A->comp[6][5] + 
           A->comp[1][0]*A->comp[4][2]*A->comp[6][5] + A->comp[0][2]*A->comp[5][0]*A->comp[6][5] - A->comp[0][0]*A->comp[5][2]*A->comp[6][5];

  /* 1,3,7 */
  test[8] = -(A->comp[0][6]*A->comp[1][2]*A->comp[2][0]) + A->comp[0][2]*A->comp[1][6]*A->comp[2][0] + A->comp[0][6]*A->comp[1][0]*A->comp[2][2] - 
           A->comp[0][0]*A->comp[1][6]*A->comp[2][2] - A->comp[0][2]*A->comp[1][0]*A->comp[2][6] + A->comp[0][0]*A->comp[1][2]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][2]*A->comp[4][0] + A->comp[0][2]*A->comp[3][6]*A->comp[4][0] + A->comp[0][6]*A->comp[3][0]*A->comp[4][2] - 
           A->comp[0][0]*A->comp[3][6]*A->comp[4][2] - A->comp[0][2]*A->comp[3][0]*A->comp[4][6] + A->comp[0][0]*A->comp[3][2]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][2]*A->comp[5][0] + A->comp[1][2]*A->comp[3][6]*A->comp[5][0] + A->comp[2][6]*A->comp[4][2]*A->comp[5][0] - 
           A->comp[2][2]*A->comp[4][6]*A->comp[5][0] + A->comp[1][6]*A->comp[3][0]*A->comp[5][2] - A->comp[1][0]*A->comp[3][6]*A->comp[5][2] - 
           A->comp[2][6]*A->comp[4][0]*A->comp[5][2] + A->comp[2][0]*A->comp[4][6]*A->comp[5][2] - A->comp[1][2]*A->comp[3][0]*A->comp[5][6] +
           A->comp[1][0]*A->comp[3][2]*A->comp[5][6] + A->comp[2][2]*A->comp[4][0]*A->comp[5][6] - A->comp[2][0]*A->comp[4][2]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][2]*A->comp[6][0] + A->comp[2][2]*A->comp[3][6]*A->comp[6][0] - A->comp[1][6]*A->comp[4][2]*A->comp[6][0] + 
           A->comp[1][2]*A->comp[4][6]*A->comp[6][0] + A->comp[0][6]*A->comp[5][2]*A->comp[6][0] - A->comp[0][2]*A->comp[5][6]*A->comp[6][0] + 
           A->comp[2][6]*A->comp[3][0]*A->comp[6][2] - A->comp[2][0]*A->comp[3][6]*A->comp[6][2] + A->comp[1][6]*A->comp[4][0]*A->comp[6][2] -  
           A->comp[1][0]*A->comp[4][6]*A->comp[6][2] - A->comp[0][6]*A->comp[5][0]*A->comp[6][2] + A->comp[0][0]*A->comp[5][6]*A->comp[6][2] - 
           A->comp[2][2]*A->comp[3][0]*A->comp[6][6] + A->comp[2][0]*A->comp[3][2]*A->comp[6][6] - A->comp[1][2]*A->comp[4][0]*A->comp[6][6] + 
           A->comp[1][0]*A->comp[4][2]*A->comp[6][6] + A->comp[0][2]*A->comp[5][0]*A->comp[6][6] - A->comp[0][0]*A->comp[5][2]*A->comp[6][6];

  /* 1.4.5 */
  test[9] = -1.0 -(A->comp[0][4]*A->comp[1][3]*A->comp[2][0]) + A->comp[0][3]*A->comp[1][4]*A->comp[2][0] + A->comp[0][4]*A->comp[1][0]*A->comp[2][3] - 
           A->comp[0][0]*A->comp[1][4]*A->comp[2][3] - A->comp[0][3]*A->comp[1][0]*A->comp[2][4] + A->comp[0][0]*A->comp[1][3]*A->comp[2][4] - 
           A->comp[0][4]*A->comp[3][3]*A->comp[4][0] + A->comp[0][3]*A->comp[3][4]*A->comp[4][0] + A->comp[0][4]*A->comp[3][0]*A->comp[4][3] - 
           A->comp[0][0]*A->comp[3][4]*A->comp[4][3] - A->comp[0][3]*A->comp[3][0]*A->comp[4][4] + A->comp[0][0]*A->comp[3][3]*A->comp[4][4] - 
           A->comp[1][4]*A->comp[3][3]*A->comp[5][0] + A->comp[1][3]*A->comp[3][4]*A->comp[5][0] + A->comp[2][4]*A->comp[4][3]*A->comp[5][0] - 
           A->comp[2][3]*A->comp[4][4]*A->comp[5][0] + A->comp[1][4]*A->comp[3][0]*A->comp[5][3] - A->comp[1][0]*A->comp[3][4]*A->comp[5][3] - 
           A->comp[2][4]*A->comp[4][0]*A->comp[5][3] + A->comp[2][0]*A->comp[4][4]*A->comp[5][3] - A->comp[1][3]*A->comp[3][0]*A->comp[5][4] +
           A->comp[1][0]*A->comp[3][3]*A->comp[5][4] + A->comp[2][3]*A->comp[4][0]*A->comp[5][4] - A->comp[2][0]*A->comp[4][3]*A->comp[5][4] - 
           A->comp[2][4]*A->comp[3][3]*A->comp[6][0] + A->comp[2][3]*A->comp[3][4]*A->comp[6][0] - A->comp[1][4]*A->comp[4][3]*A->comp[6][0] + 
           A->comp[1][3]*A->comp[4][4]*A->comp[6][0] + A->comp[0][4]*A->comp[5][3]*A->comp[6][0] - A->comp[0][3]*A->comp[5][4]*A->comp[6][0] + 
           A->comp[2][4]*A->comp[3][0]*A->comp[6][3] - A->comp[2][0]*A->comp[3][4]*A->comp[6][3] + A->comp[1][4]*A->comp[4][0]*A->comp[6][3] - 
           A->comp[1][0]*A->comp[4][4]*A->comp[6][3] - A->comp[0][4]*A->comp[5][0]*A->comp[6][3] + A->comp[0][0]*A->comp[5][4]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][0]*A->comp[6][4] + A->comp[2][0]*A->comp[3][3]*A->comp[6][4] - A->comp[1][3]*A->comp[4][0]*A->comp[6][4] + 
           A->comp[1][0]*A->comp[4][3]*A->comp[6][4] + A->comp[0][3]*A->comp[5][0]*A->comp[6][4] - A->comp[0][0]*A->comp[5][3]*A->comp[6][4];
  
  /* 1,4,6 */
  test[10] = -(A->comp[0][5]*A->comp[1][3]*A->comp[2][0]) + A->comp[0][3]*A->comp[1][5]*A->comp[2][0] + A->comp[0][5]*A->comp[1][0]*A->comp[2][3] - 
           A->comp[0][0]*A->comp[1][5]*A->comp[2][3] - A->comp[0][3]*A->comp[1][0]*A->comp[2][5] + A->comp[0][0]*A->comp[1][3]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][3]*A->comp[4][0] + A->comp[0][3]*A->comp[3][5]*A->comp[4][0] + A->comp[0][5]*A->comp[3][0]*A->comp[4][3] - 
           A->comp[0][0]*A->comp[3][5]*A->comp[4][3] - A->comp[0][3]*A->comp[3][0]*A->comp[4][5] + A->comp[0][0]*A->comp[3][3]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][3]*A->comp[5][0] + A->comp[1][3]*A->comp[3][5]*A->comp[5][0] + A->comp[2][5]*A->comp[4][3]*A->comp[5][0] - 
           A->comp[2][3]*A->comp[4][5]*A->comp[5][0] + A->comp[1][5]*A->comp[3][0]*A->comp[5][3] - A->comp[1][0]*A->comp[3][5]*A->comp[5][3] - 
           A->comp[2][5]*A->comp[4][0]*A->comp[5][3] + A->comp[2][0]*A->comp[4][5]*A->comp[5][3] - A->comp[1][3]*A->comp[3][0]*A->comp[5][5] +
           A->comp[1][0]*A->comp[3][3]*A->comp[5][5] + A->comp[2][3]*A->comp[4][0]*A->comp[5][5] - A->comp[2][0]*A->comp[4][3]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][3]*A->comp[6][0] + A->comp[2][3]*A->comp[3][5]*A->comp[6][0] - A->comp[1][5]*A->comp[4][3]*A->comp[6][0] + 
           A->comp[1][3]*A->comp[4][5]*A->comp[6][0] + A->comp[0][5]*A->comp[5][3]*A->comp[6][0] - A->comp[0][3]*A->comp[5][5]*A->comp[6][0] + 
           A->comp[2][5]*A->comp[3][0]*A->comp[6][3] - A->comp[2][0]*A->comp[3][5]*A->comp[6][3] + A->comp[1][5]*A->comp[4][0]*A->comp[6][3] - 
           A->comp[1][0]*A->comp[4][5]*A->comp[6][3] - A->comp[0][5]*A->comp[5][0]*A->comp[6][3] + A->comp[0][0]*A->comp[5][5]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][0]*A->comp[6][5] + A->comp[2][0]*A->comp[3][3]*A->comp[6][5] - A->comp[1][3]*A->comp[4][0]*A->comp[6][5] + 
           A->comp[1][0]*A->comp[4][3]*A->comp[6][5] + A->comp[0][3]*A->comp[5][0]*A->comp[6][5] - A->comp[0][0]*A->comp[5][3]*A->comp[6][5];

  /* 1,4,7 */
  test[11] = -(A->comp[0][6]*A->comp[1][3]*A->comp[2][0]) + A->comp[0][3]*A->comp[1][6]*A->comp[2][0] + A->comp[0][6]*A->comp[1][0]*A->comp[2][3] - 
           A->comp[0][0]*A->comp[1][6]*A->comp[2][3] - A->comp[0][3]*A->comp[1][0]*A->comp[2][6] + A->comp[0][0]*A->comp[1][3]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][3]*A->comp[4][0] + A->comp[0][3]*A->comp[3][6]*A->comp[4][0] + A->comp[0][6]*A->comp[3][0]*A->comp[4][3] - 
           A->comp[0][0]*A->comp[3][6]*A->comp[4][3] - A->comp[0][3]*A->comp[3][0]*A->comp[4][6] + A->comp[0][0]*A->comp[3][3]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][3]*A->comp[5][0] + A->comp[1][3]*A->comp[3][6]*A->comp[5][0] + A->comp[2][6]*A->comp[4][3]*A->comp[5][0] - 
           A->comp[2][3]*A->comp[4][6]*A->comp[5][0] + A->comp[1][6]*A->comp[3][0]*A->comp[5][3] - A->comp[1][0]*A->comp[3][6]*A->comp[5][3] - 
           A->comp[2][6]*A->comp[4][0]*A->comp[5][3] + A->comp[2][0]*A->comp[4][6]*A->comp[5][3] - A->comp[1][3]*A->comp[3][0]*A->comp[5][6] +
           A->comp[1][0]*A->comp[3][3]*A->comp[5][6] + A->comp[2][3]*A->comp[4][0]*A->comp[5][6] - A->comp[2][0]*A->comp[4][3]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][3]*A->comp[6][0] + A->comp[2][3]*A->comp[3][6]*A->comp[6][0] - A->comp[1][6]*A->comp[4][3]*A->comp[6][0] + 
           A->comp[1][3]*A->comp[4][6]*A->comp[6][0] + A->comp[0][6]*A->comp[5][3]*A->comp[6][0] - A->comp[0][3]*A->comp[5][6]*A->comp[6][0] + 
           A->comp[2][6]*A->comp[3][0]*A->comp[6][3] - A->comp[2][0]*A->comp[3][6]*A->comp[6][3] + A->comp[1][6]*A->comp[4][0]*A->comp[6][3] - 
           A->comp[1][0]*A->comp[4][6]*A->comp[6][3] - A->comp[0][6]*A->comp[5][0]*A->comp[6][3] + A->comp[0][0]*A->comp[5][6]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][0]*A->comp[6][6] + A->comp[2][0]*A->comp[3][3]*A->comp[6][6] - A->comp[1][3]*A->comp[4][0]*A->comp[6][6] + 
           A->comp[1][0]*A->comp[4][3]*A->comp[6][6] + A->comp[0][3]*A->comp[5][0]*A->comp[6][6] - A->comp[0][0]*A->comp[5][3]*A->comp[6][6];
  
  /* 1,5,6 */
  test[12] = -(A->comp[0][5]*A->comp[1][4]*A->comp[2][0]) + A->comp[0][4]*A->comp[1][5]*A->comp[2][0] + A->comp[0][5]*A->comp[1][0]*A->comp[2][4] - 
           A->comp[0][0]*A->comp[1][5]*A->comp[2][4] - A->comp[0][4]*A->comp[1][0]*A->comp[2][5] + A->comp[0][0]*A->comp[1][4]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][4]*A->comp[4][0] + A->comp[0][4]*A->comp[3][5]*A->comp[4][0] + A->comp[0][5]*A->comp[3][0]*A->comp[4][4] - 
           A->comp[0][0]*A->comp[3][5]*A->comp[4][4] - A->comp[0][4]*A->comp[3][0]*A->comp[4][5] + A->comp[0][0]*A->comp[3][4]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][4]*A->comp[5][0] + A->comp[1][4]*A->comp[3][5]*A->comp[5][0] + A->comp[2][5]*A->comp[4][4]*A->comp[5][0] - 
           A->comp[2][4]*A->comp[4][5]*A->comp[5][0] + A->comp[1][5]*A->comp[3][0]*A->comp[5][4] - A->comp[1][0]*A->comp[3][5]*A->comp[5][4] - 
           A->comp[2][5]*A->comp[4][0]*A->comp[5][4] + A->comp[2][0]*A->comp[4][5]*A->comp[5][4] - A->comp[1][4]*A->comp[3][0]*A->comp[5][5] +
           A->comp[1][0]*A->comp[3][4]*A->comp[5][5] + A->comp[2][4]*A->comp[4][0]*A->comp[5][5] - A->comp[2][0]*A->comp[4][4]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][4]*A->comp[6][0] + A->comp[2][4]*A->comp[3][5]*A->comp[6][0] - A->comp[1][5]*A->comp[4][4]*A->comp[6][0] + 
           A->comp[1][4]*A->comp[4][5]*A->comp[6][0] + A->comp[0][5]*A->comp[5][4]*A->comp[6][0] - A->comp[0][4]*A->comp[5][5]*A->comp[6][0] + 
           A->comp[2][5]*A->comp[3][0]*A->comp[6][4] - A->comp[2][0]*A->comp[3][5]*A->comp[6][4] + A->comp[1][5]*A->comp[4][0]*A->comp[6][4] - 
           A->comp[1][0]*A->comp[4][5]*A->comp[6][4] - A->comp[0][5]*A->comp[5][0]*A->comp[6][4] + A->comp[0][0]*A->comp[5][5]*A->comp[6][4] - 
           A->comp[2][4]*A->comp[3][0]*A->comp[6][5] + A->comp[2][0]*A->comp[3][4]*A->comp[6][5] - A->comp[1][4]*A->comp[4][0]*A->comp[6][5] + 
           A->comp[1][0]*A->comp[4][4]*A->comp[6][5] + A->comp[0][4]*A->comp[5][0]*A->comp[6][5] - A->comp[0][0]*A->comp[5][4]*A->comp[6][5];

  /* 1,5,7 */
  test[13] = -(A->comp[0][6]*A->comp[1][4]*A->comp[2][0]) + A->comp[0][4]*A->comp[1][6]*A->comp[2][0] + A->comp[0][6]*A->comp[1][0]*A->comp[2][4] - 
           A->comp[0][0]*A->comp[1][6]*A->comp[2][4] - A->comp[0][4]*A->comp[1][0]*A->comp[2][6] + A->comp[0][0]*A->comp[1][4]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][4]*A->comp[4][0] + A->comp[0][4]*A->comp[3][6]*A->comp[4][0] + A->comp[0][6]*A->comp[3][0]*A->comp[4][4] - 
           A->comp[0][0]*A->comp[3][6]*A->comp[4][4] - A->comp[0][4]*A->comp[3][0]*A->comp[4][6] + A->comp[0][0]*A->comp[3][4]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][4]*A->comp[5][0] + A->comp[1][4]*A->comp[3][6]*A->comp[5][0] + A->comp[2][6]*A->comp[4][4]*A->comp[5][0] - 
           A->comp[2][4]*A->comp[4][6]*A->comp[5][0] + A->comp[1][6]*A->comp[3][0]*A->comp[5][4] - A->comp[1][0]*A->comp[3][6]*A->comp[5][4] - 
           A->comp[2][6]*A->comp[4][0]*A->comp[5][4] + A->comp[2][0]*A->comp[4][6]*A->comp[5][4] - A->comp[1][4]*A->comp[3][0]*A->comp[5][6] +
           A->comp[1][0]*A->comp[3][4]*A->comp[5][6] + A->comp[2][4]*A->comp[4][0]*A->comp[5][6] - A->comp[2][0]*A->comp[4][4]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][4]*A->comp[6][0] + A->comp[2][4]*A->comp[3][6]*A->comp[6][0] - A->comp[1][6]*A->comp[4][4]*A->comp[6][0] + 
           A->comp[1][4]*A->comp[4][6]*A->comp[6][0] + A->comp[0][6]*A->comp[5][4]*A->comp[6][0] - A->comp[0][4]*A->comp[5][6]*A->comp[6][0] + 
           A->comp[2][6]*A->comp[3][0]*A->comp[6][4] - A->comp[2][0]*A->comp[3][6]*A->comp[6][4] + A->comp[1][6]*A->comp[4][0]*A->comp[6][4] - 
           A->comp[1][0]*A->comp[4][6]*A->comp[6][4] - A->comp[0][6]*A->comp[5][0]*A->comp[6][4] + A->comp[0][0]*A->comp[5][6]*A->comp[6][4] - 
           A->comp[2][4]*A->comp[3][0]*A->comp[6][6] + A->comp[2][0]*A->comp[3][4]*A->comp[6][6] - A->comp[1][4]*A->comp[4][0]*A->comp[6][6] + 
           A->comp[1][0]*A->comp[4][4]*A->comp[6][6] + A->comp[0][4]*A->comp[5][0]*A->comp[6][6] - A->comp[0][0]*A->comp[5][4]*A->comp[6][6];

  /* 1,6,7 */
  test[14] = 1.0 -(A->comp[0][6]*A->comp[1][5]*A->comp[2][0]) + A->comp[0][5]*A->comp[1][6]*A->comp[2][0] + A->comp[0][6]*A->comp[1][0]*A->comp[2][5] - 
           A->comp[0][0]*A->comp[1][6]*A->comp[2][5] - A->comp[0][5]*A->comp[1][0]*A->comp[2][6] + A->comp[0][0]*A->comp[1][5]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][5]*A->comp[4][0] + A->comp[0][5]*A->comp[3][6]*A->comp[4][0] + A->comp[0][6]*A->comp[3][0]*A->comp[4][5] - 
           A->comp[0][0]*A->comp[3][6]*A->comp[4][5] - A->comp[0][5]*A->comp[3][0]*A->comp[4][6] + A->comp[0][0]*A->comp[3][5]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][5]*A->comp[5][0] + A->comp[1][5]*A->comp[3][6]*A->comp[5][0] + A->comp[2][6]*A->comp[4][5]*A->comp[5][0] - 
           A->comp[2][5]*A->comp[4][6]*A->comp[5][0] + A->comp[1][6]*A->comp[3][0]*A->comp[5][5] - A->comp[1][0]*A->comp[3][6]*A->comp[5][5] - 
           A->comp[2][6]*A->comp[4][0]*A->comp[5][5] + A->comp[2][0]*A->comp[4][6]*A->comp[5][5] - A->comp[1][5]*A->comp[3][0]*A->comp[5][6] +
           A->comp[1][0]*A->comp[3][5]*A->comp[5][6] + A->comp[2][5]*A->comp[4][0]*A->comp[5][6] - A->comp[2][0]*A->comp[4][5]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][5]*A->comp[6][0] + A->comp[2][5]*A->comp[3][6]*A->comp[6][0] - A->comp[1][6]*A->comp[4][5]*A->comp[6][0] + 
           A->comp[1][5]*A->comp[4][6]*A->comp[6][0] + A->comp[0][6]*A->comp[5][5]*A->comp[6][0] - A->comp[0][5]*A->comp[5][6]*A->comp[6][0] + 
           A->comp[2][6]*A->comp[3][0]*A->comp[6][5] - A->comp[2][0]*A->comp[3][6]*A->comp[6][5] + A->comp[1][6]*A->comp[4][0]*A->comp[6][5] - 
           A->comp[1][0]*A->comp[4][6]*A->comp[6][5] - A->comp[0][6]*A->comp[5][0]*A->comp[6][5] + A->comp[0][0]*A->comp[5][6]*A->comp[6][5] - 
           A->comp[2][5]*A->comp[3][0]*A->comp[6][6] + A->comp[2][0]*A->comp[3][5]*A->comp[6][6] - A->comp[1][5]*A->comp[4][0]*A->comp[6][6] + 
           A->comp[1][0]*A->comp[4][5]*A->comp[6][6] + A->comp[0][5]*A->comp[5][0]*A->comp[6][6] - A->comp[0][0]*A->comp[5][5]*A->comp[6][6];

  /* 2,3,4 */
  test[15] = -(A->comp[0][3]*A->comp[1][2]*A->comp[2][1]) + A->comp[0][2]*A->comp[1][3]*A->comp[2][1] + A->comp[0][3]*A->comp[1][1]*A->comp[2][2] - 
           A->comp[0][1]*A->comp[1][3]*A->comp[2][2] - A->comp[0][2]*A->comp[1][1]*A->comp[2][3] + A->comp[0][1]*A->comp[1][2]*A->comp[2][3] - 
           A->comp[0][3]*A->comp[3][2]*A->comp[4][1] + A->comp[0][2]*A->comp[3][3]*A->comp[4][1] + A->comp[0][3]*A->comp[3][1]*A->comp[4][2] - 
           A->comp[0][1]*A->comp[3][3]*A->comp[4][2] - A->comp[0][2]*A->comp[3][1]*A->comp[4][3] + A->comp[0][1]*A->comp[3][2]*A->comp[4][3] - 
           A->comp[1][3]*A->comp[3][2]*A->comp[5][1] + A->comp[1][2]*A->comp[3][3]*A->comp[5][1] + A->comp[2][3]*A->comp[4][2]*A->comp[5][1] - 
           A->comp[2][2]*A->comp[4][3]*A->comp[5][1] + A->comp[1][3]*A->comp[3][1]*A->comp[5][2] - A->comp[1][1]*A->comp[3][3]*A->comp[5][2] - 
           A->comp[2][3]*A->comp[4][1]*A->comp[5][2] + A->comp[2][1]*A->comp[4][3]*A->comp[5][2] - A->comp[1][2]*A->comp[3][1]*A->comp[5][3] +
           A->comp[1][1]*A->comp[3][2]*A->comp[5][3] + A->comp[2][2]*A->comp[4][1]*A->comp[5][3] - A->comp[2][1]*A->comp[4][2]*A->comp[5][3] - 
           A->comp[2][3]*A->comp[3][2]*A->comp[6][1] + A->comp[2][2]*A->comp[3][3]*A->comp[6][1] - A->comp[1][3]*A->comp[4][2]*A->comp[6][1] + 
           A->comp[1][2]*A->comp[4][3]*A->comp[6][1] + A->comp[0][3]*A->comp[5][2]*A->comp[6][1] - A->comp[0][2]*A->comp[5][3]*A->comp[6][1] + 
           A->comp[2][3]*A->comp[3][1]*A->comp[6][2] - A->comp[2][1]*A->comp[3][3]*A->comp[6][2] + A->comp[1][3]*A->comp[4][1]*A->comp[6][2] - 
           A->comp[1][1]*A->comp[4][3]*A->comp[6][2] - A->comp[0][3]*A->comp[5][1]*A->comp[6][2] + A->comp[0][1]*A->comp[5][3]*A->comp[6][2] - 
           A->comp[2][2]*A->comp[3][1]*A->comp[6][3] + A->comp[2][1]*A->comp[3][2]*A->comp[6][3] - A->comp[1][2]*A->comp[4][1]*A->comp[6][3] + 
           A->comp[1][1]*A->comp[4][2]*A->comp[6][3] + A->comp[0][2]*A->comp[5][1]*A->comp[6][3] - A->comp[0][1]*A->comp[5][2]*A->comp[6][3];

  /* 2,3,5 */
  test[16] = -(A->comp[0][4]*A->comp[1][2]*A->comp[2][1]) + A->comp[0][2]*A->comp[1][4]*A->comp[2][1] + A->comp[0][4]*A->comp[1][1]*A->comp[2][2] - 
           A->comp[0][1]*A->comp[1][4]*A->comp[2][2] - A->comp[0][2]*A->comp[1][1]*A->comp[2][4] + A->comp[0][1]*A->comp[1][2]*A->comp[2][4] - 
           A->comp[0][4]*A->comp[3][2]*A->comp[4][1] + A->comp[0][2]*A->comp[3][4]*A->comp[4][1] + A->comp[0][4]*A->comp[3][1]*A->comp[4][2] - 
           A->comp[0][1]*A->comp[3][4]*A->comp[4][2] - A->comp[0][2]*A->comp[3][1]*A->comp[4][4] + A->comp[0][1]*A->comp[3][2]*A->comp[4][4] - 
           A->comp[1][4]*A->comp[3][2]*A->comp[5][1] + A->comp[1][2]*A->comp[3][4]*A->comp[5][1] + A->comp[2][4]*A->comp[4][2]*A->comp[5][1] - 
           A->comp[2][2]*A->comp[4][4]*A->comp[5][1] + A->comp[1][4]*A->comp[3][1]*A->comp[5][2] - A->comp[1][1]*A->comp[3][4]*A->comp[5][2] - 
           A->comp[2][4]*A->comp[4][1]*A->comp[5][2] + A->comp[2][1]*A->comp[4][4]*A->comp[5][2] - A->comp[1][2]*A->comp[3][1]*A->comp[5][4] +
           A->comp[1][1]*A->comp[3][2]*A->comp[5][4] + A->comp[2][2]*A->comp[4][1]*A->comp[5][4] - A->comp[2][1]*A->comp[4][2]*A->comp[5][4] - 
           A->comp[2][4]*A->comp[3][2]*A->comp[6][1] + A->comp[2][2]*A->comp[3][4]*A->comp[6][1] - A->comp[1][4]*A->comp[4][2]*A->comp[6][1] + 
           A->comp[1][2]*A->comp[4][4]*A->comp[6][1] + A->comp[0][4]*A->comp[5][2]*A->comp[6][1] - A->comp[0][2]*A->comp[5][4]*A->comp[6][1] + 
           A->comp[2][4]*A->comp[3][1]*A->comp[6][2] - A->comp[2][1]*A->comp[3][4]*A->comp[6][2] + A->comp[1][4]*A->comp[4][1]*A->comp[6][2] - 
           A->comp[1][1]*A->comp[4][4]*A->comp[6][2] - A->comp[0][4]*A->comp[5][1]*A->comp[6][2] + A->comp[0][1]*A->comp[5][4]*A->comp[6][2] - 
           A->comp[2][2]*A->comp[3][1]*A->comp[6][4] + A->comp[2][1]*A->comp[3][2]*A->comp[6][4] - A->comp[1][2]*A->comp[4][1]*A->comp[6][4] + 
           A->comp[1][1]*A->comp[4][2]*A->comp[6][4] + A->comp[0][2]*A->comp[5][1]*A->comp[6][4] - A->comp[0][1]*A->comp[5][2]*A->comp[6][4];

  /* 2,3,6 */
  test[17] = -(A->comp[0][5]*A->comp[1][2]*A->comp[2][1]) + A->comp[0][2]*A->comp[1][5]*A->comp[2][1] + A->comp[0][5]*A->comp[1][1]*A->comp[2][2] - 
           A->comp[0][1]*A->comp[1][5]*A->comp[2][2] - A->comp[0][2]*A->comp[1][1]*A->comp[2][5] + A->comp[0][1]*A->comp[1][2]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][2]*A->comp[4][1] + A->comp[0][2]*A->comp[3][5]*A->comp[4][1] + A->comp[0][5]*A->comp[3][1]*A->comp[4][2] - 
           A->comp[0][1]*A->comp[3][5]*A->comp[4][2] - A->comp[0][2]*A->comp[3][1]*A->comp[4][5] + A->comp[0][1]*A->comp[3][2]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][2]*A->comp[5][1] + A->comp[1][2]*A->comp[3][5]*A->comp[5][1] + A->comp[2][5]*A->comp[4][2]*A->comp[5][1] - 
           A->comp[2][2]*A->comp[4][5]*A->comp[5][1] + A->comp[1][5]*A->comp[3][1]*A->comp[5][2] - A->comp[1][1]*A->comp[3][5]*A->comp[5][2] - 
           A->comp[2][5]*A->comp[4][1]*A->comp[5][2] + A->comp[2][1]*A->comp[4][5]*A->comp[5][2] - A->comp[1][2]*A->comp[3][1]*A->comp[5][5] +
           A->comp[1][1]*A->comp[3][2]*A->comp[5][5] + A->comp[2][2]*A->comp[4][1]*A->comp[5][5] - A->comp[2][1]*A->comp[4][2]*A->comp[5][5] -
           A->comp[2][5]*A->comp[3][2]*A->comp[6][1] + A->comp[2][2]*A->comp[3][5]*A->comp[6][1] - A->comp[1][5]*A->comp[4][2]*A->comp[6][1] + 
           A->comp[1][2]*A->comp[4][5]*A->comp[6][1] + A->comp[0][5]*A->comp[5][2]*A->comp[6][1] - A->comp[0][2]*A->comp[5][5]*A->comp[6][1] + 
           A->comp[2][5]*A->comp[3][1]*A->comp[6][2] - A->comp[2][1]*A->comp[3][5]*A->comp[6][2] + A->comp[1][5]*A->comp[4][1]*A->comp[6][2] - 
           A->comp[1][1]*A->comp[4][5]*A->comp[6][2] - A->comp[0][5]*A->comp[5][1]*A->comp[6][2] + A->comp[0][1]*A->comp[5][5]*A->comp[6][2] - 
           A->comp[2][2]*A->comp[3][1]*A->comp[6][5] + A->comp[2][1]*A->comp[3][2]*A->comp[6][5] - A->comp[1][2]*A->comp[4][1]*A->comp[6][5] + 
           A->comp[1][1]*A->comp[4][2]*A->comp[6][5] + A->comp[0][2]*A->comp[5][1]*A->comp[6][5] - A->comp[0][1]*A->comp[5][2]*A->comp[6][5];

  /* 2,3,7 */
  test[18] = -(A->comp[0][6]*A->comp[1][2]*A->comp[2][1]) + A->comp[0][2]*A->comp[1][6]*A->comp[2][1] + A->comp[0][6]*A->comp[1][1]*A->comp[2][2] - 
           A->comp[0][1]*A->comp[1][6]*A->comp[2][2] - A->comp[0][2]*A->comp[1][1]*A->comp[2][6] + A->comp[0][1]*A->comp[1][2]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][2]*A->comp[4][1] + A->comp[0][2]*A->comp[3][6]*A->comp[4][1] + A->comp[0][6]*A->comp[3][1]*A->comp[4][2] - 
           A->comp[0][1]*A->comp[3][6]*A->comp[4][2] - A->comp[0][2]*A->comp[3][1]*A->comp[4][6] + A->comp[0][1]*A->comp[3][2]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][2]*A->comp[5][1] + A->comp[1][2]*A->comp[3][6]*A->comp[5][1] + A->comp[2][6]*A->comp[4][2]*A->comp[5][1] - 
           A->comp[2][2]*A->comp[4][6]*A->comp[5][1] + A->comp[1][6]*A->comp[3][1]*A->comp[5][2] - A->comp[1][1]*A->comp[3][6]*A->comp[5][2] - 
           A->comp[2][6]*A->comp[4][1]*A->comp[5][2] + A->comp[2][1]*A->comp[4][6]*A->comp[5][2] - A->comp[1][2]*A->comp[3][1]*A->comp[5][6] +
           A->comp[1][1]*A->comp[3][2]*A->comp[5][6] + A->comp[2][2]*A->comp[4][1]*A->comp[5][6] - A->comp[2][1]*A->comp[4][2]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][2]*A->comp[6][1] + A->comp[2][2]*A->comp[3][6]*A->comp[6][1] - A->comp[1][6]*A->comp[4][2]*A->comp[6][1] + 
           A->comp[1][2]*A->comp[4][6]*A->comp[6][1] + A->comp[0][6]*A->comp[5][2]*A->comp[6][1] - A->comp[0][2]*A->comp[5][6]*A->comp[6][1] + 
           A->comp[2][6]*A->comp[3][1]*A->comp[6][2] - A->comp[2][1]*A->comp[3][6]*A->comp[6][2] + A->comp[1][6]*A->comp[4][1]*A->comp[6][2] - 
           A->comp[1][1]*A->comp[4][6]*A->comp[6][2] - A->comp[0][6]*A->comp[5][1]*A->comp[6][2] + A->comp[0][1]*A->comp[5][6]*A->comp[6][2] - 
           A->comp[2][2]*A->comp[3][1]*A->comp[6][6] + A->comp[2][1]*A->comp[3][2]*A->comp[6][6] - A->comp[1][2]*A->comp[4][1]*A->comp[6][6] + 
           A->comp[1][1]*A->comp[4][2]*A->comp[6][6] + A->comp[0][2]*A->comp[5][1]*A->comp[6][6] - A->comp[0][1]*A->comp[5][2]*A->comp[6][6];

  /* 2,4,5 */
  test[19] = -(A->comp[0][4]*A->comp[1][3]*A->comp[2][1]) + A->comp[0][3]*A->comp[1][4]*A->comp[2][1] + A->comp[0][4]*A->comp[1][1]*A->comp[2][3] - 
           A->comp[0][1]*A->comp[1][4]*A->comp[2][3] - A->comp[0][3]*A->comp[1][1]*A->comp[2][4] + A->comp[0][1]*A->comp[1][3]*A->comp[2][4] - 
           A->comp[0][4]*A->comp[3][3]*A->comp[4][1] + A->comp[0][3]*A->comp[3][4]*A->comp[4][1] + A->comp[0][4]*A->comp[3][1]*A->comp[4][3] - 
           A->comp[0][1]*A->comp[3][4]*A->comp[4][3] - A->comp[0][3]*A->comp[3][1]*A->comp[4][4] + A->comp[0][1]*A->comp[3][3]*A->comp[4][4] - 
           A->comp[1][4]*A->comp[3][3]*A->comp[5][1] + A->comp[1][3]*A->comp[3][4]*A->comp[5][1] + A->comp[2][4]*A->comp[4][3]*A->comp[5][1] - 
           A->comp[2][3]*A->comp[4][4]*A->comp[5][1] + A->comp[1][4]*A->comp[3][1]*A->comp[5][3] - A->comp[1][1]*A->comp[3][4]*A->comp[5][3] - 
           A->comp[2][4]*A->comp[4][1]*A->comp[5][3] + A->comp[2][1]*A->comp[4][4]*A->comp[5][3] - A->comp[1][3]*A->comp[3][1]*A->comp[5][4] +
           A->comp[1][1]*A->comp[3][3]*A->comp[5][4] + A->comp[2][3]*A->comp[4][1]*A->comp[5][4] - A->comp[2][1]*A->comp[4][3]*A->comp[5][4] - 
           A->comp[2][4]*A->comp[3][3]*A->comp[6][1] + A->comp[2][3]*A->comp[3][4]*A->comp[6][1] - A->comp[1][4]*A->comp[4][3]*A->comp[6][1] + 
           A->comp[1][3]*A->comp[4][4]*A->comp[6][1] + A->comp[0][4]*A->comp[5][3]*A->comp[6][1] - A->comp[0][3]*A->comp[5][4]*A->comp[6][1] + 
           A->comp[2][4]*A->comp[3][1]*A->comp[6][3] - A->comp[2][1]*A->comp[3][4]*A->comp[6][3] + A->comp[1][4]*A->comp[4][1]*A->comp[6][3] - 
           A->comp[1][1]*A->comp[4][4]*A->comp[6][3] - A->comp[0][4]*A->comp[5][1]*A->comp[6][3] + A->comp[0][1]*A->comp[5][4]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][1]*A->comp[6][4] + A->comp[2][1]*A->comp[3][3]*A->comp[6][4] - A->comp[1][3]*A->comp[4][1]*A->comp[6][4] + 
           A->comp[1][1]*A->comp[4][3]*A->comp[6][4] + A->comp[0][3]*A->comp[5][1]*A->comp[6][4] - A->comp[0][1]*A->comp[5][3]*A->comp[6][4];

  /* 2,4,6 */
  test[20] = -1.0 -(A->comp[0][5]*A->comp[1][3]*A->comp[2][1]) + A->comp[0][3]*A->comp[1][5]*A->comp[2][1] + A->comp[0][5]*A->comp[1][1]*A->comp[2][3] - 
           A->comp[0][1]*A->comp[1][5]*A->comp[2][3] - A->comp[0][3]*A->comp[1][1]*A->comp[2][5] + A->comp[0][1]*A->comp[1][3]*A->comp[2][5] -  
           A->comp[0][5]*A->comp[3][3]*A->comp[4][1] + A->comp[0][3]*A->comp[3][5]*A->comp[4][1] + A->comp[0][5]*A->comp[3][1]*A->comp[4][3] - 
           A->comp[0][1]*A->comp[3][5]*A->comp[4][3] - A->comp[0][3]*A->comp[3][1]*A->comp[4][5] + A->comp[0][1]*A->comp[3][3]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][3]*A->comp[5][1] + A->comp[1][3]*A->comp[3][5]*A->comp[5][1] + A->comp[2][5]*A->comp[4][3]*A->comp[5][1] - 
           A->comp[2][3]*A->comp[4][5]*A->comp[5][1] + A->comp[1][5]*A->comp[3][1]*A->comp[5][3] - A->comp[1][1]*A->comp[3][5]*A->comp[5][3] - 
           A->comp[2][5]*A->comp[4][1]*A->comp[5][3] + A->comp[2][1]*A->comp[4][5]*A->comp[5][3] - A->comp[1][3]*A->comp[3][1]*A->comp[5][5] +
           A->comp[1][1]*A->comp[3][3]*A->comp[5][5] + A->comp[2][3]*A->comp[4][1]*A->comp[5][5] - A->comp[2][1]*A->comp[4][3]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][3]*A->comp[6][1] + A->comp[2][3]*A->comp[3][5]*A->comp[6][1] - A->comp[1][5]*A->comp[4][3]*A->comp[6][1] + 
           A->comp[1][3]*A->comp[4][5]*A->comp[6][1] + A->comp[0][5]*A->comp[5][3]*A->comp[6][1] - A->comp[0][3]*A->comp[5][5]*A->comp[6][1] + 
           A->comp[2][5]*A->comp[3][1]*A->comp[6][3] - A->comp[2][1]*A->comp[3][5]*A->comp[6][3] + A->comp[1][5]*A->comp[4][1]*A->comp[6][3] - 
           A->comp[1][1]*A->comp[4][5]*A->comp[6][3] - A->comp[0][5]*A->comp[5][1]*A->comp[6][3] + A->comp[0][1]*A->comp[5][5]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][1]*A->comp[6][5] + A->comp[2][1]*A->comp[3][3]*A->comp[6][5] - A->comp[1][3]*A->comp[4][1]*A->comp[6][5] + 
           A->comp[1][1]*A->comp[4][3]*A->comp[6][5] + A->comp[0][3]*A->comp[5][1]*A->comp[6][5] - A->comp[0][1]*A->comp[5][3]*A->comp[6][5];

  /* 2,4,7 */
  test[21] = -(A->comp[0][6]*A->comp[1][3]*A->comp[2][1]) + A->comp[0][3]*A->comp[1][6]*A->comp[2][1] + A->comp[0][6]*A->comp[1][1]*A->comp[2][3] - 
           A->comp[0][1]*A->comp[1][6]*A->comp[2][3] - A->comp[0][3]*A->comp[1][1]*A->comp[2][6] + A->comp[0][1]*A->comp[1][3]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][3]*A->comp[4][1] + A->comp[0][3]*A->comp[3][6]*A->comp[4][1] + A->comp[0][6]*A->comp[3][1]*A->comp[4][3] - 
           A->comp[0][1]*A->comp[3][6]*A->comp[4][3] - A->comp[0][3]*A->comp[3][1]*A->comp[4][6] + A->comp[0][1]*A->comp[3][3]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][3]*A->comp[5][1] + A->comp[1][3]*A->comp[3][6]*A->comp[5][1] + A->comp[2][6]*A->comp[4][3]*A->comp[5][1] - 
           A->comp[2][3]*A->comp[4][6]*A->comp[5][1] + A->comp[1][6]*A->comp[3][1]*A->comp[5][3] - A->comp[1][1]*A->comp[3][6]*A->comp[5][3] - 
           A->comp[2][6]*A->comp[4][1]*A->comp[5][3] + A->comp[2][1]*A->comp[4][6]*A->comp[5][3] - A->comp[1][3]*A->comp[3][1]*A->comp[5][6] +
           A->comp[1][1]*A->comp[3][3]*A->comp[5][6] + A->comp[2][3]*A->comp[4][1]*A->comp[5][6] - A->comp[2][1]*A->comp[4][3]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][3]*A->comp[6][1] + A->comp[2][3]*A->comp[3][6]*A->comp[6][1] - A->comp[1][6]*A->comp[4][3]*A->comp[6][1] + 
           A->comp[1][3]*A->comp[4][6]*A->comp[6][1] + A->comp[0][6]*A->comp[5][3]*A->comp[6][1] - A->comp[0][3]*A->comp[5][6]*A->comp[6][1] + 
           A->comp[2][6]*A->comp[3][1]*A->comp[6][3] - A->comp[2][1]*A->comp[3][6]*A->comp[6][3] + A->comp[1][6]*A->comp[4][1]*A->comp[6][3] - 
           A->comp[1][1]*A->comp[4][6]*A->comp[6][3] - A->comp[0][6]*A->comp[5][1]*A->comp[6][3] + A->comp[0][1]*A->comp[5][6]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][1]*A->comp[6][6] + A->comp[2][1]*A->comp[3][3]*A->comp[6][6] - A->comp[1][3]*A->comp[4][1]*A->comp[6][6] + 
           A->comp[1][1]*A->comp[4][3]*A->comp[6][6] + A->comp[0][3]*A->comp[5][1]*A->comp[6][6] - A->comp[0][1]*A->comp[5][3]*A->comp[6][6];

  /* 2,5,6 */
  test[22] = -(A->comp[0][5]*A->comp[1][4]*A->comp[2][1]) + A->comp[0][4]*A->comp[1][5]*A->comp[2][1] + A->comp[0][5]*A->comp[1][1]*A->comp[2][4] - 
           A->comp[0][1]*A->comp[1][5]*A->comp[2][4] - A->comp[0][4]*A->comp[1][1]*A->comp[2][5] + A->comp[0][1]*A->comp[1][4]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][4]*A->comp[4][1] + A->comp[0][4]*A->comp[3][5]*A->comp[4][1] + A->comp[0][5]*A->comp[3][1]*A->comp[4][4] - 
           A->comp[0][1]*A->comp[3][5]*A->comp[4][4] - A->comp[0][4]*A->comp[3][1]*A->comp[4][5] + A->comp[0][1]*A->comp[3][4]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][4]*A->comp[5][1] + A->comp[1][4]*A->comp[3][5]*A->comp[5][1] + A->comp[2][5]*A->comp[4][4]*A->comp[5][1] - 
           A->comp[2][4]*A->comp[4][5]*A->comp[5][1] + A->comp[1][5]*A->comp[3][1]*A->comp[5][4] - A->comp[1][1]*A->comp[3][5]*A->comp[5][4] - 
           A->comp[2][5]*A->comp[4][1]*A->comp[5][4] + A->comp[2][1]*A->comp[4][5]*A->comp[5][4] - A->comp[1][4]*A->comp[3][1]*A->comp[5][5] +
           A->comp[1][1]*A->comp[3][4]*A->comp[5][5] + A->comp[2][4]*A->comp[4][1]*A->comp[5][5] - A->comp[2][1]*A->comp[4][4]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][4]*A->comp[6][1] + A->comp[2][4]*A->comp[3][5]*A->comp[6][1] - A->comp[1][5]*A->comp[4][4]*A->comp[6][1] + 
           A->comp[1][4]*A->comp[4][5]*A->comp[6][1] + A->comp[0][5]*A->comp[5][4]*A->comp[6][1] - A->comp[0][4]*A->comp[5][5]*A->comp[6][1] + 
           A->comp[2][5]*A->comp[3][1]*A->comp[6][4] - A->comp[2][1]*A->comp[3][5]*A->comp[6][4] + A->comp[1][5]*A->comp[4][1]*A->comp[6][4] - 
           A->comp[1][1]*A->comp[4][5]*A->comp[6][4] - A->comp[0][5]*A->comp[5][1]*A->comp[6][4] + A->comp[0][1]*A->comp[5][5]*A->comp[6][4] - 
           A->comp[2][4]*A->comp[3][1]*A->comp[6][5] + A->comp[2][1]*A->comp[3][4]*A->comp[6][5] - A->comp[1][4]*A->comp[4][1]*A->comp[6][5] + 
           A->comp[1][1]*A->comp[4][4]*A->comp[6][5] + A->comp[0][4]*A->comp[5][1]*A->comp[6][5] - A->comp[0][1]*A->comp[5][4]*A->comp[6][5];

  /* 2,5,7 */
  test[23] = -1.0 -(A->comp[0][6]*A->comp[1][4]*A->comp[2][1]) + A->comp[0][4]*A->comp[1][6]*A->comp[2][1] + A->comp[0][6]*A->comp[1][1]*A->comp[2][4] - 
           A->comp[0][1]*A->comp[1][6]*A->comp[2][4] - A->comp[0][4]*A->comp[1][1]*A->comp[2][6] + A->comp[0][1]*A->comp[1][4]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][4]*A->comp[4][1] + A->comp[0][4]*A->comp[3][6]*A->comp[4][1] + A->comp[0][6]*A->comp[3][1]*A->comp[4][4] - 
           A->comp[0][1]*A->comp[3][6]*A->comp[4][4] - A->comp[0][4]*A->comp[3][1]*A->comp[4][6] + A->comp[0][1]*A->comp[3][4]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][4]*A->comp[5][1] + A->comp[1][4]*A->comp[3][6]*A->comp[5][1] + A->comp[2][6]*A->comp[4][4]*A->comp[5][1] - 
           A->comp[2][4]*A->comp[4][6]*A->comp[5][1] + A->comp[1][6]*A->comp[3][1]*A->comp[5][4] - A->comp[1][1]*A->comp[3][6]*A->comp[5][4] - 
           A->comp[2][6]*A->comp[4][1]*A->comp[5][4] + A->comp[2][1]*A->comp[4][6]*A->comp[5][4] - A->comp[1][4]*A->comp[3][1]*A->comp[5][6] +
           A->comp[1][1]*A->comp[3][4]*A->comp[5][6] + A->comp[2][4]*A->comp[4][1]*A->comp[5][6] - A->comp[2][1]*A->comp[4][4]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][4]*A->comp[6][1] + A->comp[2][4]*A->comp[3][6]*A->comp[6][1] - A->comp[1][6]*A->comp[4][4]*A->comp[6][1] + 
           A->comp[1][4]*A->comp[4][6]*A->comp[6][1] + A->comp[0][6]*A->comp[5][4]*A->comp[6][1] - A->comp[0][4]*A->comp[5][6]*A->comp[6][1] + 
           A->comp[2][6]*A->comp[3][1]*A->comp[6][4] - A->comp[2][1]*A->comp[3][6]*A->comp[6][4] + A->comp[1][6]*A->comp[4][1]*A->comp[6][4] - 
           A->comp[1][1]*A->comp[4][6]*A->comp[6][4] - A->comp[0][6]*A->comp[5][1]*A->comp[6][4] + A->comp[0][1]*A->comp[5][6]*A->comp[6][4] - 
           A->comp[2][4]*A->comp[3][1]*A->comp[6][6] + A->comp[2][1]*A->comp[3][4]*A->comp[6][6] - A->comp[1][4]*A->comp[4][1]*A->comp[6][6] + 
           A->comp[1][1]*A->comp[4][4]*A->comp[6][6] + A->comp[0][4]*A->comp[5][1]*A->comp[6][6] - A->comp[0][1]*A->comp[5][4]*A->comp[6][6];

  /* 2,6,7 */
  test[24] = -(A->comp[0][6]*A->comp[1][5]*A->comp[2][1]) + A->comp[0][5]*A->comp[1][6]*A->comp[2][1] + A->comp[0][6]*A->comp[1][1]*A->comp[2][5] - 
           A->comp[0][1]*A->comp[1][6]*A->comp[2][5] - A->comp[0][5]*A->comp[1][1]*A->comp[2][6] + A->comp[0][1]*A->comp[1][5]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][5]*A->comp[4][1] + A->comp[0][5]*A->comp[3][6]*A->comp[4][1] + A->comp[0][6]*A->comp[3][1]*A->comp[4][5] - 
           A->comp[0][1]*A->comp[3][6]*A->comp[4][5] - A->comp[0][5]*A->comp[3][1]*A->comp[4][6] + A->comp[0][1]*A->comp[3][5]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][5]*A->comp[5][1] + A->comp[1][5]*A->comp[3][6]*A->comp[5][1] + A->comp[2][6]*A->comp[4][5]*A->comp[5][1] - 
           A->comp[2][5]*A->comp[4][6]*A->comp[5][1] + A->comp[1][6]*A->comp[3][1]*A->comp[5][5] - A->comp[1][1]*A->comp[3][6]*A->comp[5][5] - 
           A->comp[2][6]*A->comp[4][1]*A->comp[5][5] + A->comp[2][1]*A->comp[4][6]*A->comp[5][5] - A->comp[1][5]*A->comp[3][1]*A->comp[5][6] +
           A->comp[1][1]*A->comp[3][5]*A->comp[5][6] + A->comp[2][5]*A->comp[4][1]*A->comp[5][6] - A->comp[2][1]*A->comp[4][5]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][5]*A->comp[6][1] + A->comp[2][5]*A->comp[3][6]*A->comp[6][1] - A->comp[1][6]*A->comp[4][5]*A->comp[6][1] + 
           A->comp[1][5]*A->comp[4][6]*A->comp[6][1] + A->comp[0][6]*A->comp[5][5]*A->comp[6][1] - A->comp[0][5]*A->comp[5][6]*A->comp[6][1] + 
           A->comp[2][6]*A->comp[3][1]*A->comp[6][5] - A->comp[2][1]*A->comp[3][6]*A->comp[6][5] + A->comp[1][6]*A->comp[4][1]*A->comp[6][5] - 
           A->comp[1][1]*A->comp[4][6]*A->comp[6][5] - A->comp[0][6]*A->comp[5][1]*A->comp[6][5] + A->comp[0][1]*A->comp[5][6]*A->comp[6][5] - 
           A->comp[2][5]*A->comp[3][1]*A->comp[6][6] + A->comp[2][1]*A->comp[3][5]*A->comp[6][6] - A->comp[1][5]*A->comp[4][1]*A->comp[6][6] + 
           A->comp[1][1]*A->comp[4][5]*A->comp[6][6] + A->comp[0][5]*A->comp[5][1]*A->comp[6][6] - A->comp[0][1]*A->comp[5][5]*A->comp[6][6];

  /* 3,4,5 */
  test[25] = -(A->comp[0][4]*A->comp[1][3]*A->comp[2][2]) + A->comp[0][3]*A->comp[1][4]*A->comp[2][2] + A->comp[0][4]*A->comp[1][2]*A->comp[2][3] - 
           A->comp[0][2]*A->comp[1][4]*A->comp[2][3] - A->comp[0][3]*A->comp[1][2]*A->comp[2][4] + A->comp[0][2]*A->comp[1][3]*A->comp[2][4] - 
           A->comp[0][4]*A->comp[3][3]*A->comp[4][2] + A->comp[0][3]*A->comp[3][4]*A->comp[4][2] + A->comp[0][4]*A->comp[3][2]*A->comp[4][3] - 
           A->comp[0][2]*A->comp[3][4]*A->comp[4][3] - A->comp[0][3]*A->comp[3][2]*A->comp[4][4] + A->comp[0][2]*A->comp[3][3]*A->comp[4][4] - 
           A->comp[1][4]*A->comp[3][3]*A->comp[5][2] + A->comp[1][3]*A->comp[3][4]*A->comp[5][2] + A->comp[2][4]*A->comp[4][3]*A->comp[5][2] - 
           A->comp[2][3]*A->comp[4][4]*A->comp[5][2] + A->comp[1][4]*A->comp[3][2]*A->comp[5][3] - A->comp[1][2]*A->comp[3][4]*A->comp[5][3] - 
           A->comp[2][4]*A->comp[4][2]*A->comp[5][3] + A->comp[2][2]*A->comp[4][4]*A->comp[5][3] - A->comp[1][3]*A->comp[3][2]*A->comp[5][4] +
           A->comp[1][2]*A->comp[3][3]*A->comp[5][4] + A->comp[2][3]*A->comp[4][2]*A->comp[5][4] - A->comp[2][2]*A->comp[4][3]*A->comp[5][4] - 
           A->comp[2][4]*A->comp[3][3]*A->comp[6][2] + A->comp[2][3]*A->comp[3][4]*A->comp[6][2] - A->comp[1][4]*A->comp[4][3]*A->comp[6][2] + 
           A->comp[1][3]*A->comp[4][4]*A->comp[6][2] + A->comp[0][4]*A->comp[5][3]*A->comp[6][2] - A->comp[0][3]*A->comp[5][4]*A->comp[6][2] + 
           A->comp[2][4]*A->comp[3][2]*A->comp[6][3] - A->comp[2][2]*A->comp[3][4]*A->comp[6][3] + A->comp[1][4]*A->comp[4][2]*A->comp[6][3] - 
           A->comp[1][2]*A->comp[4][4]*A->comp[6][3] - A->comp[0][4]*A->comp[5][2]*A->comp[6][3] + A->comp[0][2]*A->comp[5][4]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][2]*A->comp[6][4] + A->comp[2][2]*A->comp[3][3]*A->comp[6][4] - A->comp[1][3]*A->comp[4][2]*A->comp[6][4] + 
           A->comp[1][2]*A->comp[4][3]*A->comp[6][4] + A->comp[0][3]*A->comp[5][2]*A->comp[6][4] - A->comp[0][2]*A->comp[5][3]*A->comp[6][4];

  /* 3,4,6 */
  test[26] = -(A->comp[0][5]*A->comp[1][3]*A->comp[2][2]) + A->comp[0][3]*A->comp[1][5]*A->comp[2][2] + A->comp[0][5]*A->comp[1][2]*A->comp[2][3] - 
           A->comp[0][2]*A->comp[1][5]*A->comp[2][3] - A->comp[0][3]*A->comp[1][2]*A->comp[2][5] + A->comp[0][2]*A->comp[1][3]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][3]*A->comp[4][2] + A->comp[0][3]*A->comp[3][5]*A->comp[4][2] + A->comp[0][5]*A->comp[3][2]*A->comp[4][3] - 
           A->comp[0][2]*A->comp[3][5]*A->comp[4][3] - A->comp[0][3]*A->comp[3][2]*A->comp[4][5] + A->comp[0][2]*A->comp[3][3]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][3]*A->comp[5][2] + A->comp[1][3]*A->comp[3][5]*A->comp[5][2] + A->comp[2][5]*A->comp[4][3]*A->comp[5][2] - 
           A->comp[2][3]*A->comp[4][5]*A->comp[5][2] + A->comp[1][5]*A->comp[3][2]*A->comp[5][3] - A->comp[1][2]*A->comp[3][5]*A->comp[5][3] - 
           A->comp[2][5]*A->comp[4][2]*A->comp[5][3] + A->comp[2][2]*A->comp[4][5]*A->comp[5][3] - A->comp[1][3]*A->comp[3][2]*A->comp[5][5] +
           A->comp[1][2]*A->comp[3][3]*A->comp[5][5] + A->comp[2][3]*A->comp[4][2]*A->comp[5][5] - A->comp[2][2]*A->comp[4][3]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][3]*A->comp[6][2] + A->comp[2][3]*A->comp[3][5]*A->comp[6][2] - A->comp[1][5]*A->comp[4][3]*A->comp[6][2] + 
           A->comp[1][3]*A->comp[4][5]*A->comp[6][2] + A->comp[0][5]*A->comp[5][3]*A->comp[6][2] - A->comp[0][3]*A->comp[5][5]*A->comp[6][2] + 
           A->comp[2][5]*A->comp[3][2]*A->comp[6][3] - A->comp[2][2]*A->comp[3][5]*A->comp[6][3] + A->comp[1][5]*A->comp[4][2]*A->comp[6][3] - 
           A->comp[1][2]*A->comp[4][5]*A->comp[6][3] - A->comp[0][5]*A->comp[5][2]*A->comp[6][3] + A->comp[0][2]*A->comp[5][5]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][2]*A->comp[6][5] + A->comp[2][2]*A->comp[3][3]*A->comp[6][5] - A->comp[1][3]*A->comp[4][2]*A->comp[6][5] + 
           A->comp[1][2]*A->comp[4][3]*A->comp[6][5] + A->comp[0][3]*A->comp[5][2]*A->comp[6][5] - A->comp[0][2]*A->comp[5][3]*A->comp[6][5];

  /* 3,4,7 */
  test[27] = -1.0 -(A->comp[0][6]*A->comp[1][3]*A->comp[2][2]) + A->comp[0][3]*A->comp[1][6]*A->comp[2][2] + A->comp[0][6]*A->comp[1][2]*A->comp[2][3] - 
           A->comp[0][2]*A->comp[1][6]*A->comp[2][3] - A->comp[0][3]*A->comp[1][2]*A->comp[2][6] + A->comp[0][2]*A->comp[1][3]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][3]*A->comp[4][2] + A->comp[0][3]*A->comp[3][6]*A->comp[4][2] + A->comp[0][6]*A->comp[3][2]*A->comp[4][3] - 
           A->comp[0][2]*A->comp[3][6]*A->comp[4][3] - A->comp[0][3]*A->comp[3][2]*A->comp[4][6] + A->comp[0][2]*A->comp[3][3]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][3]*A->comp[5][2] + A->comp[1][3]*A->comp[3][6]*A->comp[5][2] + A->comp[2][6]*A->comp[4][3]*A->comp[5][2] - 
           A->comp[2][3]*A->comp[4][6]*A->comp[5][2] + A->comp[1][6]*A->comp[3][2]*A->comp[5][3] - A->comp[1][2]*A->comp[3][6]*A->comp[5][3] - 
           A->comp[2][6]*A->comp[4][2]*A->comp[5][3] + A->comp[2][2]*A->comp[4][6]*A->comp[5][3] - A->comp[1][3]*A->comp[3][2]*A->comp[5][6] +
           A->comp[1][2]*A->comp[3][3]*A->comp[5][6] + A->comp[2][3]*A->comp[4][2]*A->comp[5][6] - A->comp[2][2]*A->comp[4][3]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][3]*A->comp[6][2] + A->comp[2][3]*A->comp[3][6]*A->comp[6][2] - A->comp[1][6]*A->comp[4][3]*A->comp[6][2] + 
           A->comp[1][3]*A->comp[4][6]*A->comp[6][2] + A->comp[0][6]*A->comp[5][3]*A->comp[6][2] - A->comp[0][3]*A->comp[5][6]*A->comp[6][2] + 
           A->comp[2][6]*A->comp[3][2]*A->comp[6][3] - A->comp[2][2]*A->comp[3][6]*A->comp[6][3] + A->comp[1][6]*A->comp[4][2]*A->comp[6][3] - 
           A->comp[1][2]*A->comp[4][6]*A->comp[6][3] - A->comp[0][6]*A->comp[5][2]*A->comp[6][3] + A->comp[0][2]*A->comp[5][6]*A->comp[6][3] - 
           A->comp[2][3]*A->comp[3][2]*A->comp[6][6] + A->comp[2][2]*A->comp[3][3]*A->comp[6][6] - A->comp[1][3]*A->comp[4][2]*A->comp[6][6] + 
           A->comp[1][2]*A->comp[4][3]*A->comp[6][6] + A->comp[0][3]*A->comp[5][2]*A->comp[6][6] - A->comp[0][2]*A->comp[5][3]*A->comp[6][6];

  /* 3,5,6 */
  test[28] = 1.0 -(A->comp[0][5]*A->comp[1][4]*A->comp[2][2]) + A->comp[0][4]*A->comp[1][5]*A->comp[2][2] + A->comp[0][5]*A->comp[1][2]*A->comp[2][4] - 
           A->comp[0][2]*A->comp[1][5]*A->comp[2][4] - A->comp[0][4]*A->comp[1][2]*A->comp[2][5] + A->comp[0][2]*A->comp[1][4]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][4]*A->comp[4][2] + A->comp[0][4]*A->comp[3][5]*A->comp[4][2] + A->comp[0][5]*A->comp[3][2]*A->comp[4][4] - 
           A->comp[0][2]*A->comp[3][5]*A->comp[4][4] - A->comp[0][4]*A->comp[3][2]*A->comp[4][5] + A->comp[0][2]*A->comp[3][4]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][4]*A->comp[5][2] + A->comp[1][4]*A->comp[3][5]*A->comp[5][2] + A->comp[2][5]*A->comp[4][4]*A->comp[5][2] - 
           A->comp[2][4]*A->comp[4][5]*A->comp[5][2] + A->comp[1][5]*A->comp[3][2]*A->comp[5][4] - A->comp[1][2]*A->comp[3][5]*A->comp[5][4] - 
           A->comp[2][5]*A->comp[4][2]*A->comp[5][4] + A->comp[2][2]*A->comp[4][5]*A->comp[5][4] - A->comp[1][4]*A->comp[3][2]*A->comp[5][5] +
           A->comp[1][2]*A->comp[3][4]*A->comp[5][5] + A->comp[2][4]*A->comp[4][2]*A->comp[5][5] - A->comp[2][2]*A->comp[4][4]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][4]*A->comp[6][2] + A->comp[2][4]*A->comp[3][5]*A->comp[6][2] - A->comp[1][5]*A->comp[4][4]*A->comp[6][2] + 
           A->comp[1][4]*A->comp[4][5]*A->comp[6][2] + A->comp[0][5]*A->comp[5][4]*A->comp[6][2] - A->comp[0][4]*A->comp[5][5]*A->comp[6][2] + 
           A->comp[2][5]*A->comp[3][2]*A->comp[6][4] - A->comp[2][2]*A->comp[3][5]*A->comp[6][4] + A->comp[1][5]*A->comp[4][2]*A->comp[6][4] - 
           A->comp[1][2]*A->comp[4][5]*A->comp[6][4] - A->comp[0][5]*A->comp[5][2]*A->comp[6][4] + A->comp[0][2]*A->comp[5][5]*A->comp[6][4] - 
           A->comp[2][4]*A->comp[3][2]*A->comp[6][5] + A->comp[2][2]*A->comp[3][4]*A->comp[6][5] - A->comp[1][4]*A->comp[4][2]*A->comp[6][5] + 
           A->comp[1][2]*A->comp[4][4]*A->comp[6][5] + A->comp[0][4]*A->comp[5][2]*A->comp[6][5] - A->comp[0][2]*A->comp[5][4]*A->comp[6][5];

  /* 3,5,7 */
  test[29] = -(A->comp[0][6]*A->comp[1][4]*A->comp[2][2]) + A->comp[0][4]*A->comp[1][6]*A->comp[2][2] + A->comp[0][6]*A->comp[1][2]*A->comp[2][4] - 
           A->comp[0][2]*A->comp[1][6]*A->comp[2][4] - A->comp[0][4]*A->comp[1][2]*A->comp[2][6] + A->comp[0][2]*A->comp[1][4]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][4]*A->comp[4][2] + A->comp[0][4]*A->comp[3][6]*A->comp[4][2] + A->comp[0][6]*A->comp[3][2]*A->comp[4][4] - 
           A->comp[0][2]*A->comp[3][6]*A->comp[4][4] - A->comp[0][4]*A->comp[3][2]*A->comp[4][6] + A->comp[0][2]*A->comp[3][4]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][4]*A->comp[5][2] + A->comp[1][4]*A->comp[3][6]*A->comp[5][2] + A->comp[2][6]*A->comp[4][4]*A->comp[5][2] - 
           A->comp[2][4]*A->comp[4][6]*A->comp[5][2] + A->comp[1][6]*A->comp[3][2]*A->comp[5][4] - A->comp[1][2]*A->comp[3][6]*A->comp[5][4] - 
           A->comp[2][6]*A->comp[4][2]*A->comp[5][4] + A->comp[2][2]*A->comp[4][6]*A->comp[5][4] - A->comp[1][4]*A->comp[3][2]*A->comp[5][6] +
           A->comp[1][2]*A->comp[3][4]*A->comp[5][6] + A->comp[2][4]*A->comp[4][2]*A->comp[5][6] - A->comp[2][2]*A->comp[4][4]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][4]*A->comp[6][2] + A->comp[2][4]*A->comp[3][6]*A->comp[6][2] - A->comp[1][6]*A->comp[4][4]*A->comp[6][2] + 
           A->comp[1][4]*A->comp[4][6]*A->comp[6][2] + A->comp[0][6]*A->comp[5][4]*A->comp[6][2] - A->comp[0][4]*A->comp[5][6]*A->comp[6][2] + 
           A->comp[2][6]*A->comp[3][2]*A->comp[6][4] - A->comp[2][2]*A->comp[3][6]*A->comp[6][4] + A->comp[1][6]*A->comp[4][2]*A->comp[6][4] - 
           A->comp[1][2]*A->comp[4][6]*A->comp[6][4] - A->comp[0][6]*A->comp[5][2]*A->comp[6][4] + A->comp[0][2]*A->comp[5][6]*A->comp[6][4] - 
           A->comp[2][4]*A->comp[3][2]*A->comp[6][6] + A->comp[2][2]*A->comp[3][4]*A->comp[6][6] - A->comp[1][4]*A->comp[4][2]*A->comp[6][6] + 
           A->comp[1][2]*A->comp[4][4]*A->comp[6][6] + A->comp[0][4]*A->comp[5][2]*A->comp[6][6] - A->comp[0][2]*A->comp[5][4]*A->comp[6][6];


  /* 3,6,7 */
  test[30] = -(A->comp[0][6]*A->comp[1][5]*A->comp[2][2]) + A->comp[0][5]*A->comp[1][6]*A->comp[2][2] + A->comp[0][6]*A->comp[1][2]*A->comp[2][5] - 
           A->comp[0][2]*A->comp[1][6]*A->comp[2][5] - A->comp[0][5]*A->comp[1][2]*A->comp[2][6] + A->comp[0][2]*A->comp[1][5]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][5]*A->comp[4][2] + A->comp[0][5]*A->comp[3][6]*A->comp[4][2] + A->comp[0][6]*A->comp[3][2]*A->comp[4][5] - 
           A->comp[0][2]*A->comp[3][6]*A->comp[4][5] - A->comp[0][5]*A->comp[3][2]*A->comp[4][6] + A->comp[0][2]*A->comp[3][5]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][5]*A->comp[5][2] + A->comp[1][5]*A->comp[3][6]*A->comp[5][2] + A->comp[2][6]*A->comp[4][5]*A->comp[5][2] - 
           A->comp[2][5]*A->comp[4][6]*A->comp[5][2] + A->comp[1][6]*A->comp[3][2]*A->comp[5][5] - A->comp[1][2]*A->comp[3][6]*A->comp[5][5] - 
           A->comp[2][6]*A->comp[4][2]*A->comp[5][5] + A->comp[2][2]*A->comp[4][6]*A->comp[5][5] - A->comp[1][5]*A->comp[3][2]*A->comp[5][6] +
           A->comp[1][2]*A->comp[3][5]*A->comp[5][6] + A->comp[2][5]*A->comp[4][2]*A->comp[5][6] - A->comp[2][2]*A->comp[4][5]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][5]*A->comp[6][2] + A->comp[2][5]*A->comp[3][6]*A->comp[6][2] - A->comp[1][6]*A->comp[4][5]*A->comp[6][2] + 
           A->comp[1][5]*A->comp[4][6]*A->comp[6][2] + A->comp[0][6]*A->comp[5][5]*A->comp[6][2] - A->comp[0][5]*A->comp[5][6]*A->comp[6][2] + 
           A->comp[2][6]*A->comp[3][2]*A->comp[6][5] - A->comp[2][2]*A->comp[3][6]*A->comp[6][5] + A->comp[1][6]*A->comp[4][2]*A->comp[6][5] - 
           A->comp[1][2]*A->comp[4][6]*A->comp[6][5] - A->comp[0][6]*A->comp[5][2]*A->comp[6][5] + A->comp[0][2]*A->comp[5][6]*A->comp[6][5] - 
           A->comp[2][5]*A->comp[3][2]*A->comp[6][6] + A->comp[2][2]*A->comp[3][5]*A->comp[6][6] - A->comp[1][5]*A->comp[4][2]*A->comp[6][6] + 
           A->comp[1][2]*A->comp[4][5]*A->comp[6][6] + A->comp[0][5]*A->comp[5][2]*A->comp[6][6] - A->comp[0][2]*A->comp[5][5]*A->comp[6][6];

  /* 4,5,6 */
  test[31] = -(A->comp[0][5]*A->comp[1][4]*A->comp[2][3]) + A->comp[0][4]*A->comp[1][5]*A->comp[2][3] + A->comp[0][5]*A->comp[1][3]*A->comp[2][4] - 
           A->comp[0][3]*A->comp[1][5]*A->comp[2][4] - A->comp[0][4]*A->comp[1][3]*A->comp[2][5] + A->comp[0][3]*A->comp[1][4]*A->comp[2][5] - 
           A->comp[0][5]*A->comp[3][4]*A->comp[4][3] + A->comp[0][4]*A->comp[3][5]*A->comp[4][3] + A->comp[0][5]*A->comp[3][3]*A->comp[4][4] - 
           A->comp[0][3]*A->comp[3][5]*A->comp[4][4] - A->comp[0][4]*A->comp[3][3]*A->comp[4][5] + A->comp[0][3]*A->comp[3][4]*A->comp[4][5] - 
           A->comp[1][5]*A->comp[3][4]*A->comp[5][3] + A->comp[1][4]*A->comp[3][5]*A->comp[5][3] + A->comp[2][5]*A->comp[4][4]*A->comp[5][3] - 
           A->comp[2][4]*A->comp[4][5]*A->comp[5][3] + A->comp[1][5]*A->comp[3][3]*A->comp[5][4] - A->comp[1][3]*A->comp[3][5]*A->comp[5][4] - 
           A->comp[2][5]*A->comp[4][3]*A->comp[5][4] + A->comp[2][3]*A->comp[4][5]*A->comp[5][4] - A->comp[1][4]*A->comp[3][3]*A->comp[5][5] +
           A->comp[1][3]*A->comp[3][4]*A->comp[5][5] + A->comp[2][4]*A->comp[4][3]*A->comp[5][5] - A->comp[2][3]*A->comp[4][4]*A->comp[5][5] - 
           A->comp[2][5]*A->comp[3][4]*A->comp[6][3] + A->comp[2][4]*A->comp[3][5]*A->comp[6][3] - A->comp[1][5]*A->comp[4][4]*A->comp[6][3] + 
           A->comp[1][4]*A->comp[4][5]*A->comp[6][3] + A->comp[0][5]*A->comp[5][4]*A->comp[6][3] - A->comp[0][4]*A->comp[5][5]*A->comp[6][3] + 
           A->comp[2][5]*A->comp[3][3]*A->comp[6][4] - A->comp[2][3]*A->comp[3][5]*A->comp[6][4] + A->comp[1][5]*A->comp[4][3]*A->comp[6][4] - 
           A->comp[1][3]*A->comp[4][5]*A->comp[6][4] - A->comp[0][5]*A->comp[5][3]*A->comp[6][4] + A->comp[0][3]*A->comp[5][5]*A->comp[6][4] - 
           A->comp[2][4]*A->comp[3][3]*A->comp[6][5] + A->comp[2][3]*A->comp[3][4]*A->comp[6][5] - A->comp[1][4]*A->comp[4][3]*A->comp[6][5] + 
           A->comp[1][3]*A->comp[4][4]*A->comp[6][5] + A->comp[0][4]*A->comp[5][3]*A->comp[6][5] - A->comp[0][3]*A->comp[5][4]*A->comp[6][5];

  /* 4,5,7 */
  test[32] = -(A->comp[0][6]*A->comp[1][4]*A->comp[2][3]) + A->comp[0][4]*A->comp[1][6]*A->comp[2][3] + A->comp[0][6]*A->comp[1][3]*A->comp[2][4] - 
           A->comp[0][3]*A->comp[1][6]*A->comp[2][4] - A->comp[0][4]*A->comp[1][3]*A->comp[2][6] + A->comp[0][3]*A->comp[1][4]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][4]*A->comp[4][3] + A->comp[0][4]*A->comp[3][6]*A->comp[4][3] + A->comp[0][6]*A->comp[3][3]*A->comp[4][4] - 
           A->comp[0][3]*A->comp[3][6]*A->comp[4][4] - A->comp[0][4]*A->comp[3][3]*A->comp[4][6] + A->comp[0][3]*A->comp[3][4]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][4]*A->comp[5][3] + A->comp[1][4]*A->comp[3][6]*A->comp[5][3] + A->comp[2][6]*A->comp[4][4]*A->comp[5][3] - 
           A->comp[2][4]*A->comp[4][6]*A->comp[5][3] + A->comp[1][6]*A->comp[3][3]*A->comp[5][4] - A->comp[1][3]*A->comp[3][6]*A->comp[5][4] - 
           A->comp[2][6]*A->comp[4][3]*A->comp[5][4] + A->comp[2][3]*A->comp[4][6]*A->comp[5][4] - A->comp[1][4]*A->comp[3][3]*A->comp[5][6] +
           A->comp[1][3]*A->comp[3][4]*A->comp[5][6] + A->comp[2][4]*A->comp[4][3]*A->comp[5][6] - A->comp[2][3]*A->comp[4][4]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][4]*A->comp[6][3] + A->comp[2][4]*A->comp[3][6]*A->comp[6][3] - A->comp[1][6]*A->comp[4][4]*A->comp[6][3] + 
           A->comp[1][4]*A->comp[4][6]*A->comp[6][3] + A->comp[0][6]*A->comp[5][4]*A->comp[6][3] - A->comp[0][4]*A->comp[5][6]*A->comp[6][3] + 
           A->comp[2][6]*A->comp[3][3]*A->comp[6][4] - A->comp[2][3]*A->comp[3][6]*A->comp[6][4] + A->comp[1][6]*A->comp[4][3]*A->comp[6][4] - 
           A->comp[1][3]*A->comp[4][6]*A->comp[6][4] - A->comp[0][6]*A->comp[5][3]*A->comp[6][4] + A->comp[0][3]*A->comp[5][6]*A->comp[6][4] - 
           A->comp[2][4]*A->comp[3][3]*A->comp[6][6] + A->comp[2][3]*A->comp[3][4]*A->comp[6][6] - A->comp[1][4]*A->comp[4][3]*A->comp[6][6] + 
           A->comp[1][3]*A->comp[4][4]*A->comp[6][6] + A->comp[0][4]*A->comp[5][3]*A->comp[6][6] - A->comp[0][3]*A->comp[5][4]*A->comp[6][6];

  /* 4,6,7 */
  test[33] = -(A->comp[0][6]*A->comp[1][5]*A->comp[2][3]) + A->comp[0][5]*A->comp[1][6]*A->comp[2][3] + A->comp[0][6]*A->comp[1][3]*A->comp[2][5] - 
           A->comp[0][3]*A->comp[1][6]*A->comp[2][5] - A->comp[0][5]*A->comp[1][3]*A->comp[2][6] + A->comp[0][3]*A->comp[1][5]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][5]*A->comp[4][3] + A->comp[0][5]*A->comp[3][6]*A->comp[4][3] + A->comp[0][6]*A->comp[3][3]*A->comp[4][5] - 
           A->comp[0][3]*A->comp[3][6]*A->comp[4][5] - A->comp[0][5]*A->comp[3][3]*A->comp[4][6] + A->comp[0][3]*A->comp[3][5]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][5]*A->comp[5][3] + A->comp[1][5]*A->comp[3][6]*A->comp[5][3] + A->comp[2][6]*A->comp[4][5]*A->comp[5][3] - 
           A->comp[2][5]*A->comp[4][6]*A->comp[5][3] + A->comp[1][6]*A->comp[3][3]*A->comp[5][5] - A->comp[1][3]*A->comp[3][6]*A->comp[5][5] - 
           A->comp[2][6]*A->comp[4][3]*A->comp[5][5] + A->comp[2][3]*A->comp[4][6]*A->comp[5][5] - A->comp[1][5]*A->comp[3][3]*A->comp[5][6] +
           A->comp[1][3]*A->comp[3][5]*A->comp[5][6] + A->comp[2][5]*A->comp[4][3]*A->comp[5][6] - A->comp[2][3]*A->comp[4][5]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][5]*A->comp[6][3] + A->comp[2][5]*A->comp[3][6]*A->comp[6][3] - A->comp[1][6]*A->comp[4][5]*A->comp[6][3] + 
           A->comp[1][5]*A->comp[4][6]*A->comp[6][3] + A->comp[0][6]*A->comp[5][5]*A->comp[6][3] - A->comp[0][5]*A->comp[5][6]*A->comp[6][3] + 
           A->comp[2][6]*A->comp[3][3]*A->comp[6][5] - A->comp[2][3]*A->comp[3][6]*A->comp[6][5] + A->comp[1][6]*A->comp[4][3]*A->comp[6][5] - 
           A->comp[1][3]*A->comp[4][6]*A->comp[6][5] - A->comp[0][6]*A->comp[5][3]*A->comp[6][5] + A->comp[0][3]*A->comp[5][6]*A->comp[6][5] - 
           A->comp[2][5]*A->comp[3][3]*A->comp[6][6] + A->comp[2][3]*A->comp[3][5]*A->comp[6][6] - A->comp[1][5]*A->comp[4][3]*A->comp[6][6] + 
           A->comp[1][3]*A->comp[4][5]*A->comp[6][6] + A->comp[0][5]*A->comp[5][3]*A->comp[6][6] - A->comp[0][3]*A->comp[5][5]*A->comp[6][6];
  /* 5,6,7 */
  test[34] = -(A->comp[0][6]*A->comp[1][5]*A->comp[2][3]) + A->comp[0][5]*A->comp[1][6]*A->comp[2][3] + A->comp[0][6]*A->comp[1][3]*A->comp[2][5] - 
           A->comp[0][3]*A->comp[1][6]*A->comp[2][5] - A->comp[0][5]*A->comp[1][3]*A->comp[2][6] + A->comp[0][3]*A->comp[1][5]*A->comp[2][6] - 
           A->comp[0][6]*A->comp[3][5]*A->comp[4][3] + A->comp[0][5]*A->comp[3][6]*A->comp[4][3] + A->comp[0][6]*A->comp[3][3]*A->comp[4][5] - 
           A->comp[0][3]*A->comp[3][6]*A->comp[4][5] - A->comp[0][5]*A->comp[3][3]*A->comp[4][6] + A->comp[0][3]*A->comp[3][5]*A->comp[4][6] - 
           A->comp[1][6]*A->comp[3][5]*A->comp[5][3] + A->comp[1][5]*A->comp[3][6]*A->comp[5][3] + A->comp[2][6]*A->comp[4][5]*A->comp[5][3] - 
           A->comp[2][5]*A->comp[4][6]*A->comp[5][3] + A->comp[1][6]*A->comp[3][3]*A->comp[5][5] - A->comp[1][3]*A->comp[3][6]*A->comp[5][5] - 
           A->comp[2][6]*A->comp[4][3]*A->comp[5][5] + A->comp[2][3]*A->comp[4][6]*A->comp[5][5] - A->comp[1][5]*A->comp[3][3]*A->comp[5][6] +
           A->comp[1][3]*A->comp[3][5]*A->comp[5][6] + A->comp[2][5]*A->comp[4][3]*A->comp[5][6] - A->comp[2][3]*A->comp[4][5]*A->comp[5][6] - 
           A->comp[2][6]*A->comp[3][5]*A->comp[6][3] + A->comp[2][5]*A->comp[3][6]*A->comp[6][3] - A->comp[1][6]*A->comp[4][5]*A->comp[6][3] + 
           A->comp[1][5]*A->comp[4][6]*A->comp[6][3] + A->comp[0][6]*A->comp[5][5]*A->comp[6][3] - A->comp[0][5]*A->comp[5][6]*A->comp[6][3] + 
           A->comp[2][6]*A->comp[3][3]*A->comp[6][5] - A->comp[2][3]*A->comp[3][6]*A->comp[6][5] + A->comp[1][6]*A->comp[4][3]*A->comp[6][5] - 
           A->comp[1][3]*A->comp[4][6]*A->comp[6][5] - A->comp[0][6]*A->comp[5][3]*A->comp[6][5] + A->comp[0][3]*A->comp[5][6]*A->comp[6][5] - 
           A->comp[2][5]*A->comp[3][3]*A->comp[6][6] + A->comp[2][3]*A->comp[3][5]*A->comp[6][6] - A->comp[1][5]*A->comp[4][3]*A->comp[6][6] + 
           A->comp[1][3]*A->comp[4][5]*A->comp[6][6] + A->comp[0][5]*A->comp[5][3]*A->comp[6][6] - A->comp[0][3]*A->comp[5][5]*A->comp[6][6];

  test_result = 0.0;
  for (i=0; i<35; i++) 
      {
      test_result +=test[i]*test[i];
      }
  test_result = sqrt(test_result);

  *cc=test_result;

  if(test_result<MIN_VALUE) return 1;
  else return 0;
  }


/* verify that the matrix is in G2: 
 return 1 if is in G2, 0 otherwise */
int check_G2(G2 const * const A, double *n, double *cc)
    {
    int unitarity, cubic_condition;
    G2 M, N;

    /* check unitarity */
    equal_G2(&M, A);  /* M=A */
    times_dag2_G2(&N, &M, &M); /* N=M*M^{dag} */
    one_G2(&M);
    minus_equal_G2(&N, &M); /* N=A*A^{dag}-1 */
    *n=norm_G2(&N);
    if(*n<MIN_VALUE) 
      {
      unitarity=1;
      }
    else 
      {
      unitarity=0;
      }

    /* check cubic condition */
    cubic_condition=cc_check_G2(A, cc);
    /* note that the cubic condition implies det=+1 */

    return unitarity*cubic_condition;
    }

#endif
