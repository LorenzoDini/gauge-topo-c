#ifndef GAUGE_INCLUDE_H
#define GAUGE_INCLUDE_H

#include"../config.h"

#ifdef ONE_FILE_MODE

  #include"./Const/const.h"
  #include"./Func_Point/function_pointers.h"
  #include"./G2/g2.h"
  #include"./G2/g2_check.h"
  #include"./G2/g2_un.h"
  #include"./G2/g2_upd.h"
  #include"./Gauge_Conf/gauge_conf.h"
  #include"./Geometry/geometry.h"
  #include"./Rng/random.h"
  #include"./Rng/dSFMT/dSFMT.h"
  #include"./Su2/su2.h"
  #include"./Su2/su2_upd.h"

  #include"./Const/const.c"
  #include"./Func_Point/function_pointers.c"
  #include"./G2/g2.c"
  #include"./G2/g2_check.c"
  #include"./G2/g2_un.c"
  #include"./G2/g2_upd.c"
  #include"./Gauge_Conf/gauge_conf_def.c"
  #include"./Gauge_Conf/gauge_conf_meas.c" 
  #include"./Gauge_Conf/gauge_conf_upd.c"
  #include"./Geometry/geometry.c"
  #include"./Rng/random.c"
  #include"./Rng/dSFMT/dSFMT.c"
  #include"./Su2/su2.c"
  #include"./Su2/su2_upd.c"

#else

  #include"./Const/const.h"
  #include"./Func_Point/function_pointers.h"
  #include"./Gauge_Conf/gauge_conf.h"
  #include"./Geometry/geometry.h"
  #include"./Macro/macro.h"
  #include"./Rng/random.h"
#endif

#endif
