#ifndef MACRO_H
#define MACRO_H

#include"../../config.h" 

/* the choice of the group is done with configure
#define SuN_GROUP
*/

#ifdef G2_GROUP
  #define GAUGE_GROUP G2
  #define RAND_GAUGE_TRANSF
#endif

#ifdef Su2_GROUP
  #define GAUGE_GROUP Su2
#endif

#ifdef SuN_GROUP
  #define GAUGE_GROUP SuN
#endif

#ifdef SoN_GROUP
  #define GAUGE_GROUP SoN
#endif



#define INT_ALIGN 8
#define DOUBLE_ALIGN 16

#define MIN_VALUE 1.0e-13


#define PI 3.141592653589793238462643383279502884197169399375105820974944
#define PI2 6.283185307179586476925286766559005768394338798750211641949889
#define HALF_PI 1.570796326794896619231321691639751442098584699687552910487472

/* way to print a macro: if 
#define val1 val2 
then QUOTEME(val1) give the string "val2" */
#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)

#endif
