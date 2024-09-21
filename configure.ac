dnl Process this file with autoconf to produce a configure script.

dnl --- Initialize ----------------------------------------------------------- 
AC_PREREQ([2.71])
AC_INIT([Runuran],[0.39],[unuran@statmath.wu.ac.at])
AC_CONFIG_HEADERS([src/config.h])

dnl --- Configure flags ------------------------------------------------------
dnl UNU.RAN allows to write information about distribution and a log file.
dnl The default name of the logfile is set in 'src/unuran_config.h'.
dnl The logfile is enabled by the '--enable-logging' flag.
AC_ARG_ENABLE(logging,
	[AS_HELP_STRING([--enable-logging],
		[Debug: print informations about generator into logfile @<:@default=no@:>@])],
	[],
	[enable_logging=no])
AS_IF([test "x$enable_logging" = xyes],
	[AC_DEFINE([UNUR_ENABLE_LOGGING], [1], 
	           [Define to 1 if you want to use a logfile for logging informations about UNU.RAN objects].)], 
	[])

dnl UNU.RAN provides a function that returns a string with some information
dnl about a given generator object.
dnl This function is enabled by the '--enable-info' flag.
AC_ARG_ENABLE(info,
	[AS_HELP_STRING([--enable-info],
		[Info: provide function with information about generator objects @<:@default=yes@:>@])],
	[],
	[enable_info=yes])
AS_IF([test "x$enable_info" = xyes],
	[AC_DEFINE([UNUR_ENABLE_INFO], [1], 
	           [Define to 1 if you want to use the info routine for printing informations about UNU.RAN objects.])], 
	[])

dnl --- Compiler and compiler flags ------------------------------------------
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi

export CC=`"${R_HOME}/bin/R" CMD config CC`
## deprecated: export CPP=`"${R_HOME}/bin/R" CMD config CPP`
export CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
export CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP

dnl --- Enable Rmath library in UNU.RAN -------------------------------------
AH_BOTTOM([/* Use R functions in UNU.RAN. */
#define HAVE_LIBRMATH 1])

AH_BOTTOM([/* Tell UNU.RAN that we are within R. */
#define R_UNURAN 1])

dnl --- Checks for libraries -------------------------------------------------
dnl libm: the standard math library
AC_CHECK_LIB([m], [sqrt])

dnl --- Checks for header files ----------------------------------------------
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h unistd.h])
AC_CHECK_HEADERS_ONCE([sys/time.h])

dnl --- Checks for typedefs, structures, and compiler characteristics --------
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

dnl --- Checks for library functions -----------------------------------------
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strcasecmp strchr strtol strtoul])

AC_CHECK_DECLS(log1p,,,[#include <math.h>])
AH_BOTTOM([#if !HAVE_DECL_LOG1P
#define log1p _unur_log1p
#endif])

AC_CHECK_DECLS(hypot,,,[#include <math.h>])
AH_BOTTOM([#if !HAVE_DECL_HYPOT
#define hypot _unur_hypot
#endif])

AC_CHECK_DECLS([getopt, snprintf, vsnprintf])

dnl Check for C99 functions */
AC_CHECK_DECLS(isfinite,,,[#include <math.h>])
AC_CHECK_DECLS(isinf,,,[#include <math.h>])
AC_CHECK_DECLS(isnan,,,[#include <math.h>])
AC_CHECK_DECLS(INFINITY,,,[#include <math.h>])
AC_CHECK_DECLS(HUGE_VAL,,,[#include <math.h>])

dnl Check IEEE comparisons, whether "x != x" is true for NaNs
AC_CACHE_CHECK([for IEEE comparisons], ac_cv_c_ieee_comparisons,
[AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <math.h>
#include <stdlib.h>
int main (void) 
{ 
   int status; double nan;
   const double zero = 0.;
   nan = 0.0 / zero;
   status = (nan == nan);
   exit (status);
}]])],
[ac_cv_c_ieee_comparisons="yes"],
[ac_cv_c_ieee_comparisons="no"],
[ac_cv_c_ieee_comparisons="yes"])
])
if test "$ac_cv_c_ieee_comparisons" != no ; then
  AC_DEFINE([HAVE_IEEE_COMPARISONS], [1], 
      [Define to 1 if "x != x" is true for NaNs])
  AC_SUBST(HAVE_IEEE_COMPARISONS)
fi

dnl --- Make src/Makevars ----------------------------------------------------
UNURAN_SRC=`cd ./src; ls -1 unuran-src/*/*.c | tr '\012' ' '`
AC_SUBST(UNURAN_SRC)
AC_CONFIG_FILES([src/Makevars])

dnl --- Create output files --------------------------------------------------
AC_OUTPUT

dnl --------------------------------------------------------------------------
AC_MSG_RESULT([
*========================================================
*
* Configuration for Runuran:
*
*  Source directory:          ${srcdir}
*  C Compiler:                ${CC}
*
*  Enable logging:            ${enable_logging}
*  Enable info routine:       ${enable_info}
*
*========================================================
])

dnl --- End ------------------------------------------------------------------
