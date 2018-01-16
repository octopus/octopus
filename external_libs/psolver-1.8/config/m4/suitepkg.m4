# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Luigi Genovese, Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
# Usage AX_PACKAGE(NAME,version,soname,staticlinkline,cflags,compileprogram,linkroutine)
# (example of Psolver 
# AX_PACKAGE([PSOLVER],[1.8],[-lPSolver-1],[$LINALG_LIBS $LIB_FUTILE_LIBS],[$LIB_FLIB_CFLAGS],
#[program main
#    use psbase
#    use box
#    use iobox
#    use psbox
#    use poisson_solver
#  
#    write(*,*) PS_getVersion()
#  end program],
#[
#  use Poisson_solver
#  
#  type(coulomb_operator) :: kernel
#  real(dp), dimension(9) :: rhopot, potion
#  real(gp) :: eh
#  
#  call H_Potential("G", kernel, rhopot, potion, eh, 0._dp, .false.)]
#
AC_DEFUN([AX_PACKAGE],
[dnl Test for PSolver
  define([lcv],[translit([[$1]], [A-Z], [a-z])])
  AC_ARG_WITH(lcv-libs, AS_HELP_STRING([--with-lcv-libs], [Give the linker flags for an external lcv modules (default = None).]), [ax_$1_libs=$withval], [ax_$1_libs=])
  AC_ARG_WITH(lcv-incs, AS_HELP_STRING([--with-lcv-incs], [Give the compiler include flags for an external lcv library (default = None).]), [ax_$1_incdir=$withval], [ax_$1_incdir=])

  AC_LANG_PUSH(Fortran)
  dnl AC_REQUIRE([AC_PROG_FC])
  dnl AC_REQUIRE([PKG_PROG_PKG_CONFIG])

  AC_MSG_CHECKING([for $1 via PKG_CONFIG]) 
  dnl try first with pkg-config
  PKG_CHECK_MODULES([$1],
                    [lcv >= $2],
                    [ax_have_$1=yes],
                    [ax_have_$1=no])
  if test "$ax_have_$1" = "yes" ; then
    if test -z "${$1_CFLAGS// }" -a -n "$C_INCLUDE_PATH" ; then
      for path in ${C_INCLUDE_PATH//:/ }; do
        ax_$1_incdir="$ax_$1_incdir -I$path"
      done
      LIB_$1_CFLAGS=$ax_$1_incdir
    else
      LIB_$1_CFLAGS=$$1_CFLAGS
    fi
    LIB_$1_LIBS=$$1_LIBS
  fi
  AC_MSG_RESULT(
  "LIB_$1_LIBS= $LIB_$1_LIBS" 
  "LIB_$1_CFLAGS=$LIB_$1_CFLAGS" "have=" $ax_have_$1 "hello=" lcv)

  dnl try by hand search if failed
  if test "$ax_have_$1" != "yes" ; then
    dnl Test the modules for compilation
    
    dnl Test the modules for compilation
    AC_MSG_CHECKING([for lcv modules])
    FCFLAGS_SVG=$FCFLAGS
    if test -n "$ax_$1_incdir" ; then
      FCFLAGS="$FCFLAGS $ax_$1_incdir"
    elif test -n "$C_INCLUDE_PATH" ; then
      for path in ${C_INCLUDE_PATH//:/ }; do
        ax_$1_incdir="$ax_$1_incdir -I$path"
      done
      FCFLAGS="$FCFLAGS $ax_$1_incdir"
    fi
    FCFLAGS="$FCFLAGS $5"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE($6)], with$1mod=yes, with$1mod=no)
    AC_MSG_RESULT($with$1mod)
  
    dnl Test the library.
    AC_MSG_CHECKING([for lcv library])
    LIBS_SVG=$LIBS
    if test -z "$ax_$1_libs" ; then
      ax_$1_libs="$3"
    fi
    LIBS="$ax_$1_libs $LIBS_SVG"
    AC_LINK_IFELSE(
      AC_LANG_PROGRAM([], [[$7]]),
      [ax_have_$1=yes],
      [ax_have_$1=no])
    if test $ax_have_$1 != "yes" ; then
      dnl Static case, need to link with additional libs.
      ax_$1_libs="$ax_$1_libs $4"
      LIBS="$ax_$1_libs $LIBS_SVG"
      AC_LINK_IFELSE(
        AC_LANG_PROGRAM([], [[$7]]),
        [ax_have_$1=yes],
        [ax_have_$1=no])
    fi
    AC_MSG_RESULT($ax_have_$1)
  
    LIBS=$LIBS_SVG
    FCFLAGS=$FCFLAGS_SVG
    AC_LANG_POP(Fortran)
    
    if test "$ax_have_$1" = "yes" -a "$with$1mod" = "yes" ; then
      LIB_$1_CFLAGS=$ax_$1_incdir
      LIB_$1_LIBS=$ax_$1_libs
      ax_have_$1="yes"
    else
      ax_have_$1="no"
    fi
  fi

  dnl eventually control if the library is statically linked or not
  LIB_$1_DYNAMIC_LIBS=$3
  if test "$LIB_$1_LIBS" = "$LIB_$1_DYNAMIC_LIBS"; then
    ax_$1_static="no"
  else
    ax_$1_static="yes"
  fi

  AC_SUBST(LIB_$1_CFLAGS)
  AC_SUBST(LIB_$1_LIBS)
  AC_SUBST(LIB_$1_DYNAMIC_LIBS)   
])
