AC_PREREQ([2.59])
AC_INIT([g2_topo_c], [1.1], [bonati@df.unipi.it bonati@pi.infn.it])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([./src/g2_topo.c])
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


# Checks for programs.
AC_PROG_CC


# Checks for simd extensions
AS_IF([test x${simd_check_mode} = xyes], [AX_EXT])


# Checks for header files.
AC_CHECK_HEADERS([math.h stdio.h stdlib.h time.h])


#check for posix_memalign, otherwise malloc is used
AC_CHECK_FUNCS([posix_memalign])


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
