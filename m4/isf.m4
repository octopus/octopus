## Copyright (C) 2013 J. Alberdi-Rodriguez
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$
##

AC_DEFUN([ACX_ISF], [
AC_REQUIRE([ACX_MPI])
acx_isf_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(isf-prefix, [AS_HELP_STRING([--with-isf-prefix=<lib>], [From version 1.7-dev.27, it is possible to compile only the Poisson Solver in a separate library, by using --disable-libbigdft in conjunction with --disable-binaries. The compilation will then enter in the PSolver/ subdirectory, compile libPSolver-1.a and stop. More information: http://bigdft.org])])
case $with_isf_prefix in
  yes | "") ;;
  no) acx_isf_ok=disable ;;
  *.a | *.so | *.so.* | *.o) LIBS_ISF=$with_isf_prefix ;
     xpath=${with_isf_prefix%/lib/*} 
     FCFLAGS_ISF="$ax_cv_f90_modflag$xpath/include";;  
  *) LIBS_ISF="-L$with_isf_prefix/lib"; 
     FCFLAGS_ISF="$ax_cv_f90_modflag$with_isf_prefix/include" ;;
esac

dnl The include dir must be specified when the library is given with a 
dnl specified file to be compiled static (i.e. *.a etc.)
AC_ARG_WITH(isf-include, [AS_HELP_STRING([--with-isf-include=DIR], [ISF Fortran include files directory])])
case $with_isf_include in
  "") if test "x$FCFLAGS_ISF" == x; then
        FCFLAGS_ISF="$ax_cv_f90_modflag/usr/include"
      fi;;
  *)  FCFLAGS_ISF="$ax_cv_f90_modflag$with_isf_include" ;;
esac

dnl We cannot use ISF if MPI is not found
if test "x$acx_mpi_ok" != xyes; then
  acx_isf_ok=nompi
fi

dnl Backup LIBS and FCFLAGS
acx_isf_save_LIBS="$LIBS"
acx_isf_save_FCFLAGS="$FCFLAGS"

FCFLAGS_ISF="$FCFLAGS_ISF"
FCFLAGS="$FCFLAGS_ISF $acx_isf_save_FCFLAGS"

# some symbols below will not be defined for version 1.0.4, making sure
# we have a version that is able to work in our code
testprogram="AC_LANG_PROGRAM([],[ 
    use Poisson_Solver

    implicit none
    type(coulomb_operator) :: pkernel
    call pkernel_set(pkernel,.true.)

  ])"


dnl First, check LIBS_ISF environment variable
if test x"$acx_isf_ok" = xno; then
  LIBS="$LIBS_ISF $LIBS_LAPACK $LIBS_BLAS $acx_isf_save_LIB"
  AC_MSG_CHECKING([for isf library])
  AC_LINK_IFELSE($testprogram, [acx_isf_ok=yes; LIBS_ISF="$LIBS_ISF "], [])
  if test $acx_isf_ok = no; then
    AC_MSG_RESULT([$acx_isf_ok])
  else
    AC_MSG_RESULT([$acx_isf_ok ($LIBS_ISF)])
  fi
fi


dnl Generic ISF library 
if test $acx_isf_ok = no; then
  AC_MSG_CHECKING([for isf library with -lPSolver-1])
  if test "$LIBS_ISF" = ""; then
    LIBS="-lPSolver-1 -lwrappers -lflib -labinit $LIBS_LAPACK $LIBS_BLAS $LIBS $acx_isf_save_LIB"
    AC_LINK_IFELSE($testprogram, [acx_isf_ok=yes; LIBS_ISF="-lPSolver-1 -lwrappers -lflib -labinit "], [])
  else
    LIBS="$LIBS_ISF -lPSolver-1 -lwrappers -lflib -labinit $LIBS_LAPACK $LIBS_BLAS $acx_isf_save_LIB"
    AC_LINK_IFELSE($testprogram, [acx_isf_ok=yes; 
                                  LIBS_ISF="$LIBS_ISF -lPSolver-1 -lwrappers -lflib -labinit "], [])  
  fi
  if test $acx_isf_ok = no; then
    AC_MSG_RESULT([$acx_isf_ok])
  else
    AC_MSG_RESULT([$acx_isf_ok ($LIBS_ISF)])
  fi
fi


dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_isf_ok" = xyes; then
  AC_DEFINE(HAVE_LIBISF,1,[Defined if you have ISF library.])
  $1
else
  AC_MSG_WARN([Could not find ISF library. 
               *** Will compile without ISF support])
  LIBS_ISF=""
  FCFLAGS_ISF=""
  $2
fi

AC_SUBST(LIBS_ISF)
AC_SUBST(FCFLAGS_ISF)
LIBS="$acx_isf_save_LIBS $LIBS_ISF"
FCFLAGS="$acx_isf_save_FCFLAGS $FCFLAGS_ISF"

])dnl ACX_ISF
