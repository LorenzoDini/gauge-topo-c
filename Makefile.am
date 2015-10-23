bin_PROGRAMS = debug_g2 debug_rng debug_su2 debug_sun gauge_topo


if ONE_FILE_COMPILATION
  gauge_topo_SOURCES=./src/gauge_topo.c ./src/gauge_topo_include.h
  EXTRA_DIST = \
    ./src/input.in \
    ./src/Analisi/ANALISI_LOC.cc ./src/Analisi/analisi_loc.py \
    ./src/Analisi/ANALISI_TOP.cc ./src/Analisi/analisi_top.py \
    ./src/Analisi/AUTOCORR_SELF.cc ./src/Analisi/autocorr_self.py \
    ./src/Analisi/ordina.py \
    ./src/Analisi/random.cc \
    ./src/Analisi/unisci.py \
    ./src/Const/const.c ./src/Const/const.h \
    ./src/Endian/endianness.c ./src/Endian/endianness.h \
    ./src/Func_Point/function_pointers.c ./src/Func_Point/function_pointers.h \
    ./src/G2/g2.c ./src/G2/g2.h \
    ./src/G2/g2_check.c ./src/G2/g2_check.h \
    ./src/G2/g2_upd.c ./src/G2/g2_upd.h \
    ./src/Gauge_Conf/gauge_conf_def.c ./src/Gauge_Conf/gauge_conf_meas.c \
    ./src/Gauge_Conf/gauge_conf_upd.c ./src/Gauge_Conf/gauge_conf.h \
    ./src/Geometry/geometry.c ./src/Geometry/geometry.h \
    ./src/Macro/macro.h \
    ./src/Rng/random.c ./src/Rng/random.h ./src/Rng/myrngparam.h \
    $(dSFMT) \
    ./src/Su2/su2.c ./src/Su2/su2.h \
    ./src/Su2/su2_upd.c ./src/Su2/su2_upd.h \
    ./src/SuN/sun.c ./src/SuN/sun.h \
    ./src/SuN/sun_aux.c ./src/SuN/sun_aux.h \
    ./src/SuN/sun_upd.c ./src/SuN/sun_upd.h 

  DISTCHECK_CONFIGURE_FLAGS = --enable-one-file-compilation
else
  gauge_topo_SOURCES=./src/gauge_topo.c ./src/gauge_topo_include.h \
    ./src/Const/const.c ./src/Const/const.h \
    ./src/Endian/endianness.c ./src/Endian/endianness.h \
    ./src/Func_Point/function_pointers.c ./src/Func_Point/function_pointers.h \
    ./src/G2/g2.c ./src/G2/g2.h \
    ./src/G2/g2_check.c ./src/G2/g2_check.h \
    ./src/G2/g2_un.c ./src/G2/g2_un.h \
    ./src/G2/g2_upd.c ./src/G2/g2_upd.h \
    ./src/Gauge_Conf/gauge_conf_def.c ./src/Gauge_Conf/gauge_conf_meas.c \
    ./src/Gauge_Conf/gauge_conf_upd.c ./src/Gauge_Conf/gauge_conf.h \
    ./src/Geometry/geometry.c ./src/Geometry/geometry.h \
    ./src/Macro/macro.h \
    ./src/Rng/random.c ./src/Rng/random.h ./src/Rng/myrngparam.h \
    $(dSFMT) \
    ./src/Su2/su2.c ./src/Su2/su2.h \
    ./src/Su2/su2_upd.c ./src/Su2/su2_upd.h \
    ./src/SuN/sun.c ./src/SuN/sun.h \
    ./src/SuN/sun_aux.c ./src/SuN/sun_aux.h \
    ./src/SuN/sun_upd.c ./src/SuN/sun_upd.h 
  EXTRA_DIST = \
    ./src/input.in \
    ./src/Analisi/ANALISI_LOC.cc ./src/Analisi/analisi_loc.py \
    ./src/Analisi/ANALISI_TOP.cc ./src/Analisi/analisi_top.py \
    ./src/Analisi/AUTOCORR_SELF.cc ./src/Analisi/autocorr_self.py \
    ./src/Analisi/ordina.py \
    ./src/Analisi/random.cc \
    ./src/Analisi/unisci.py 
endif


debug_g2_SOURCES=./src/G2/debug_g2.c \
  ./src/Const/const.c ./src/Const/const.h \
  ./src/Endian/endianness.c ./src/Endian/endianness.h \
  ./src/G2/g2.c ./src/G2/g2.h \
  ./src/G2/g2_check.c ./src/G2/g2_check.h \
  ./src/G2/g2_un.c ./src/G2/g2_un.h \
  ./src/G2/g2_upd.c ./src/G2/g2_upd.h \
  ./src/Macro/macro.h \
  ./src/Su2/su2.c ./src/Su2/su2.h \
  ./src/Su2/su2_upd.c ./src/Su2/su2_upd.h \
  ./src/Rng/random.c ./src/Rng/random.h ./src/Rng/myrngparam.h \
  $(dSFMT)


debug_rng_SOURCES=./src/Rng/debug_rng.c ./src/Rng/random.c ./src/Rng/random.h ./src/Rng/myrngparam.h \
  $(dSFMT)


debug_su2_SOURCES=./src/Su2/debug_su2.c \
  ./src/Const/const.c ./src/Const/const.h \
  ./src/Endian/endianness.c ./src/Endian/endianness.h \
  ./src/Macro/macro.h \
  ./src/Su2/su2.c ./src/Su2/su2.h \
  ./src/Su2/su2_upd.c ./src/Su2/su2_upd.h \
  ./src/Rng/random.c ./src/Rng/random.h ./src/Rng/myrngparam.h \
  $(dSFMT)


debug_sun_SOURCES=./src/SuN/debug_sun.c \
  ./src/Const/const.c ./src/Const/const.h \
  ./src/Endian/endianness.c ./src/Endian/endianness.h \
  ./src/Macro/macro.h \
  ./src/Su2/su2.c ./src/Su2/su2.h \
  ./src/Su2/su2_upd.c ./src/Su2/su2_upd.h \
  ./src/SuN/sun.c ./src/SuN/sun.h \
  ./src/SuN/sun_aux.c ./src/SuN/sun_aux.h \
  ./src/SuN/sun_upd.c ./src/SuN/sun_upd.h \
  ./src/Rng/random.c ./src/Rng/random.h ./src/Rng/myrngparam.h \
  $(dSFMT)


dSFMT=./src/Rng/dSFMT/dSFMT.c ./src/Rng/dSFMT/dSFMT-common.h ./src/Rng/dSFMT/dSFMT.h \
  ./src/Rng/dSFMT/dSFMT-params11213.h ./src/Rng/dSFMT/dSFMT-params1279.h \
  ./src/Rng/dSFMT/dSFMT-params132049.h ./src/Rng/dSFMT/dSFMT-params19937.h \
  ./src/Rng/dSFMT/dSFMT-params216091.h ./src/Rng/dSFMT/dSFMT-params2203.h \
  ./src/Rng/dSFMT/dSFMT-params4253.h ./src/Rng/dSFMT/dSFMT-params44497.h \
  ./src/Rng/dSFMT/dSFMT-params521.h ./src/Rng/dSFMT/dSFMT-params86243.h \
  ./src/Rng/dSFMT/dSFMT-params.h


LDADD=-lm -lcrypto -lssl 


AM_CFLAGS = -std=c99 -O3 -Wall
