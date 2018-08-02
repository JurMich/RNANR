# autoconfigutation file
# provided by Ronny Lorenz

AC_INIT(RNANR,1.1.0)

AC_CONFIG_MACRO_DIR([.])

#Output this variables to the makefiles

AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_USE_SYSTEM_EXTENSIONS

AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

#Check headers
AC_CONFIG_HEADERS(config.h) 

AC_CANONICAL_HOST

#check for compilers c++ and c
AC_PROG_CC
AC_PROG_CPP

AC_CHECK_FUNCS([strchr strdup strtol])
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_TYPE_SIZE_T

AC_SEARCH_LIBS([exp], [m], [], [
  AC_MSG_ERROR([unable to find the exp() function])
])

# Search for the RNAlib2 and set variables for Makefiles if found.
# _CFLAGS and _LIBS can be overriden with the module name from
# PKG_CHECK_MODULES
PKG_CHECK_MODULES([VRNA], [RNAlib2],
                  [ AC_SUBST([VRNA_CFLAGS])
                    AC_SUBST([VRNA_LIBS])
                  ],
                  [])

AC_ARG_WITH([mpfr], AS_HELP_STRING([--without-mpfr], [do not build with MPFR support]))
AS_IF([ test "x$with_mpfr" != "xno"],[with_mpfr=yes])

if test "x$with_mpfr" != "xno"; then
    AC_CHECK_HEADER([mpfr.h], [], [
    AC_MSG_WARN([
==========================
Failed to find mpfr.h!

You probably need to install the mpfr-devel package or similar
==========================
    ])
    with_mpfr=no])
fi

if test "x$with_mpfr" != "xno"; then
    AC_SEARCH_LIBS([mpfr_sub], [mpfr], [
        GMP_LIBS="-lgmp"
        AC_SUBST([GMP_LIBS])
        AC_DEFINE([RNANR_WITH_MPFR], [1], [Use MPFR library support])
    ], [
    AC_MSG_WARN([
==========================
Failed to link against libmpfr!

You probably need to install the mpfr package or similar
==========================
    ])
    with_mpfr=no
    ])
fi

#    PKG_CHECK_MODULES([MPFR], [mpfr],
#                      [ AC_SUBST([MPFR_CFLAGS])
#                        AC_SUBST([MPFR_LIBS])
#                        AC_DEFINE([RNANR_WITH_MPFR], [1], [Use MPFR library support])
#                      ],
#                      [
#    AC_MSG_WARN([
#==========================
#MPFR library not found...
#==========================
#    ])
#    with_mpfr=no
#                      ])
#fi


AM_CONDITIONAL(RNANR_AM_SWITCH_MPFR, test "x$with_mpfr" != "xno")

AC_OUTPUT(Makefile src/Makefile)

