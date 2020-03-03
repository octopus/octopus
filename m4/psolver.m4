## Copyright (C) 2020 M. Oliveira
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
##

AC_DEFUN([ACX_PSOLVER], [
## AC_REQUIRE([ACX_MPI])
AC_REQUIRE([ACX_YAML])
AC_REQUIRE([ACX_FUTILE])
AC_REQUIRE([ACX_BLAS])
acx_psolver_ok=no

dnl Check if the library was given in the command line
dnl if not, use environment variables or defaults
AC_ARG_WITH(psolver-prefix, [AS_HELP_STRING([--with-psolver-prefix=DIR], [Directory where PSolver was installed.])])

# Set FCFLAGS_PSOLVER only if not set from environment
if test x"$FCFLAGS_PSOLVER" = x; then
  case $with_psolver_prefix in
    "") FCFLAGS_PSOLVER="$ax_cv_f90_modflag/usr/include" ;;
    *)  FCFLAGS_PSOLVER="$ax_cv_f90_modflag$with_psolver_prefix/include" ;;
  esac
fi

AC_ARG_WITH(psolver-include, [AS_HELP_STRING([--with-psolver-include=DIR], [Directory where PSolver Fortran headers were installed.])])
case $with_psolver_include in
  "") ;;
  *)  FCFLAGS_PSOLVER="$ax_cv_f90_modflag$with_psolver_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_psolver_save_LIBS="$LIBS"
acx_psolver_save_FCFLAGS="$FCFLAGS"


dnl The test
AC_MSG_CHECKING([for PSolver])

testprogram="AC_LANG_PROGRAM([],[
    use Poisson_Solver

    implicit none
    type(coulomb_operator) :: pkernel
    real(dp) :: eps(1,1,1)

    call pkernel_set(pkernel,eps=eps)

    call pkernel_free(pkernel)
  ])"


FCFLAGS="$FCFLAGS_PSOLVER $acx_psolver_save_FCFLAGS"

if test -z "$LIBS_PSOLVER"; then
  if test ! -z "$with_psolver_prefix"; then
    LIBS_PSOLVER="-L$with_psolver_prefix/lib"
  else
    LIBS_PSOLVER=""
  fi
  LIBS_PSOLVER="$LIBS_PSOLVER -lPSolver-1"
fi

LIBS="$LIBS_PSOLVER $LIBS_FUTILE $LIBS_LIBYAML $LIBS_BLAS $acx_psolver_save_LIBS"

AC_LINK_IFELSE($testprogram, [acx_psolver_ok=yes], [])

AC_MSG_RESULT([$acx_psolver_ok ($FCFLAGS_PSOLVER $LIBS_PSOLVER)])


dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_psolver_ok" = xyes; then
  AC_DEFINE(HAVE_PSOLVER, 1, [Defined if you have the PSolver library.])
else
  AC_MSG_WARN([Could not find PSolver library.])
  LIBS_PSOLVER=""
  FCFLAGS_PSOLVER=""
fi

AC_SUBST(LIBS_PSOLVER)
AC_SUBST(FCFLAGS_PSOLVER)
LIBS="$acx_psolver_save_LIBS"
FCFLAGS="$acx_psolver_save_FCFLAGS"

])dnl ACX_PSOLVER
