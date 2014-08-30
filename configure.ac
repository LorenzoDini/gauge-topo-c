AC_PREREQ([2.59])
AC_INIT([gauge_topo_c], [1.2], [bonati@df.unipi.it bonati@pi.infn.it])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([./src/gauge_topo.c])
AC_CONFIG_HEADERS([config.h])

# overwrite default values
if test x${ac_cv_env_CFLAGS_set} != xset; then
  CFLAGS=''
fi 
if test x${ac_cv_env_CPPFLAGS_set} != xset; then
  CPPFLAGS=''
fi 


# to enable one-file-compiation
AC_ARG_ENABLE([one-file-compilation], [  --enable-one-file-compilation     to enable compilation using only one file], [one_file_mode=yes], [one_file_mode=no])
if test x${one_file_mode} = xyes; then
  AC_DEFINE(ONE_FILE_MODE,, [Activate the one file compilation])
fi
AM_CONDITIONAL([ONE_FILE_COMPILATION], [test x${one_file_mode} = xyes])


# to disable simd check
AC_ARG_ENABLE([simd-check], [  --disable-simd-check     to disable simd check], [simd_check_mode=no], [simd_check_mode=yes])
AM_CONDITIONAL([SIMD_CHECK_COMPILATION], [test x${simd_check_mode} = xyes])


# Choice of the gauge group
AC_ARG_VAR(GAUGEGROUP, choices are Su2 SuN G2)
AC_ARG_VAR(N_c, the number of color for SuN)


# Checks for programs.
AC_PROG_CC


# Checks for simd extensions
AS_IF([test x${simd_check_mode} = xyes], [AX_EXT])


# Checks for header files.
AC_CHECK_HEADERS([complex.h math.h stdio.h stdlib.h time.h])


#check for posix_memalign, otherwise malloc is used
AC_CHECK_FUNCS([posix_memalign])

#define the group-related macro
if test x${GAUGEGROUP} = xSu2; then
  AC_DEFINE(Su2_GROUP,, [Uses the gauge group SU(2)])
  AC_DEFINE(NCOLOR, 2, [Number of colors])
else
  if test x${GAUGEGROUP} = xG2; then
    AC_DEFINE(G2_GROUP,, [Uses the gauge group G2])
    AC_DEFINE(NCOLOR, 7, [Number of colors])
  else
    if test x${GAUGEGROUP} = xSuN; then
      AC_DEFINE(SuN_GROUP,, [Uses the gauge group SU(N)])
      if test x${N_c} = x; then
        AC_MSG_WARN(No chosen N_c: using the default N_c=2)
        AC_DEFINE(NCOLOR, 2, [Number of colors])
      else
        AC_DEFINE_UNQUOTED(NCOLOR, ${N_c}, [Number of colors])
      fi
    else
      AC_MSG_WARN(No chosen GAUGEGROUP: using the default Su2)
      AC_DEFINE(Su2_GROUP,, [Uses the gauge group SU(2)])
      AC_DEFINE(NCOLOR, 2, [Number of colors])
    fi
  fi
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
