# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#
AC_DEFUN([AX_PSOLVER],
[dnl Test for PSolver
AC_REQUIRE([AX_FLIB])
AC_REQUIRE([AX_LINALG])
AC_REQUIRE([AX_MPI])
AX_PACKAGE([PSOLVER],[1.8],[-lPSolver-1],[$LINALG_LIBS $LIB_FUTILE_LIBS],[$LIB_FLIB_CFLAGS],
             [program main
    use psbase
    use box
    use iobox
    use psbox
    use poisson_solver
  
    write(*,*) PS_getVersion()
  end program],
         [use Poisson_solver
  type(coulomb_operator) :: kernel
  real(dp), dimension(9) :: rhopot, potion
  real(gp) :: eh
  
  call H_Potential("G", kernel, rhopot, potion, eh, 0._dp, .false.)
])
dnl   AC_ARG_WITH(psolver-libs, AS_HELP_STRING([--with-psolver-libs], [Give the linker flags for an external Psolver modules (default = None).]), ax_psolver_libs=$withval, ax_psolver_libs=)
dnl   AC_ARG_WITH(psolver-incs, AS_HELP_STRING([--with-psolver-incs], [Give the compiler include flags for an external Psolver library (default = None).]), ax_psolver_incdir=$withval, ax_psolver_incdir=)
dnl   
dnl   dnl try first with pkg-config
dnl   PKG_CHECK_MODULES([PSOLVER],
dnl                     [psolver >= 1.8],
dnl                     [ax_have_psolver=yes],
dnl                     [ax_have_psolver=no])
dnl   if test "$ax_have_psolver" = "yes" ; then
dnl     if test -z "${PSOLVER_CFLAGS// }" -a -n "$C_INCLUDE_PATH" ; then
dnl       for path in ${C_INCLUDE_PATH//:/ }; do
dnl         ax_psolver_incdir="$ax_psolver_incdir -I$path"
dnl       done
dnl       LIB_PSOLVER_CFLAGS=$ax_psolver_incdir
dnl     else
dnl       LIB_PSOLVER_CFLAGS=$PSOLVER_CFLAGS
dnl     fi
dnl     LIB_PSOLVER_LIBS=$PSOLVER_LIBS
dnl   fi
dnl 
dnl   dnl try by hand search if failed
dnl   if test "$ax_have_psolver" != "yes" ; then
dnl     dnl Test the modules for compilation
dnl     AC_LANG_PUSH(Fortran)
dnl     AC_REQUIRE([AC_PROG_FC])
dnl     AC_REQUIRE([AX_FLIB])
dnl     AC_REQUIRE([AX_LINALG])
dnl     AC_REQUIRE([AX_MPI])
dnl     
dnl     dnl Test the modules for compilation
dnl     AC_MSG_CHECKING([for PSolver modules])
dnl     FCFLAGS_SVG=$FCFLAGS
dnl     if test -n "$ax_psolver_incdir" ; then
dnl       FCFLAGS="$FCFLAGS $ax_psolver_incdir"
dnl     elif test -n "$C_INCLUDE_PATH" ; then
dnl       for path in ${C_INCLUDE_PATH//:/ }; do
dnl         ax_psolver_incdir="$ax_psolver_incdir -I$path"
dnl       done
dnl       FCFLAGS="$FCFLAGS $ax_psolver_incdir"
dnl     fi
dnl     FCFLAGS="$FCFLAGS $LIB_FLIB_CFLAGS"
dnl     AC_COMPILE_IFELSE([[program main
dnl     use psbase
dnl     use box
dnl     use iobox
dnl     use psbox
dnl     use poisson_solver
dnl   
dnl     write(*,*) PS_getVersion()
dnl   end program]], withpsolvermod=yes, withpsolvermod=no)
dnl     AC_MSG_RESULT($withpsolvermod)
dnl   
dnl     dnl Test the psolver library.
dnl     AC_MSG_CHECKING([for PSolver library])
dnl     LIBS_SVG=$LIBS
dnl     if test -z "$ax_psolver_libs" ; then
dnl       ax_psolver_libs="-lPSolver-1"
dnl     fi
dnl     LIBS="$ax_psolver_libs $LIBS_SVG"
dnl     AC_LINK_IFELSE(
dnl       AC_LANG_PROGRAM([], [[
dnl   use Poisson_solver
dnl   
dnl   type(coulomb_operator) :: kernel
dnl   real(dp), dimension(9) :: rhopot, potion
dnl   real(gp) :: eh
dnl   
dnl   call H_Potential("G", kernel, rhopot, potion, eh, 0._dp, .false.)
dnl   ]]),
dnl       [ax_have_psolver=yes],
dnl       [ax_have_psolver=no])
dnl     if test $ax_have_psolver != "yes" ; then
dnl       dnl Static case, need to link with additional libs.
dnl       ax_psolver_libs="$ax_psolver_libs $LINALG_LIBS $LIB_FUTILE_LIBS"
dnl       LIBS="$ax_psolver_libs $LIBS_SVG"
dnl       AC_LINK_IFELSE(
dnl         AC_LANG_PROGRAM([], [[
dnl   use Poisson_solver
dnl   
dnl   type(coulomb_operator) :: kernel
dnl   real(dp), dimension(9) :: rhopot, potion
dnl   real(gp) :: eh
dnl   
dnl   call H_Potential("G", kernel, rhopot, potion, eh, 0._dp, .false.)
dnl   ]]),
dnl         [ax_have_psolver=yes],
dnl         [ax_have_psolver=no])
dnl     fi
dnl     AC_MSG_RESULT($ax_have_psolver)
dnl   
dnl     LIBS=$LIBS_SVG
dnl     FCFLAGS=$FCFLAGS_SVG
dnl     AC_LANG_POP(Fortran)
dnl     
dnl     if test "$ax_have_psolver" = "yes" -a "$withpsolvermod" = "yes" ; then
dnl       LIB_PSOLVER_CFLAGS=$ax_psolver_incdir
dnl       LIB_PSOLVER_LIBS=$ax_psolver_libs
dnl       ax_have_psolver="yes"
dnl     else
dnl       ax_have_psolver="no"
dnl     fi
dnl   fi
dnl   
dnl   dnl LIB_XC_CFLAGS="-I/usr/include"
dnl   dnl   PKG_CHECK_MODULES(LIB_XC, psolver >= 2.0, ax_have_psolver="yes", ax_have_psolver="no")
dnl   
dnl   AC_SUBST(LIB_PSOLVER_CFLAGS)
dnl   AC_SUBST(LIB_PSOLVER_LIBS)
])
