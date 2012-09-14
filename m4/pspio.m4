## Copyright (C) 2012 M. Oliveira
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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id: pspio.m4 754 2012-08-31 17:05:18Z micael $

AC_DEFUN([ACX_PSPIO], [
acx_pspio_ok=no
FCFLAGS_PSPIO=""
LIBS_PSPIO=""

dnl Backup LIBS and FCFLAGS
acx_pspio_save_LIBS="$LIBS"
acx_pspio_save_FCFLAGS="$FCFLAGS"

if test "$PKGCONFIG" != ""; then
  PSPIO_PREFIX=`$PKGCONFIG --variable=prefix pspio`
else
  PSPIO_PREFIX=/usr
fi
dnl Check if the library was given in the command line
AC_ARG_WITH(pspio-prefix, [AS_HELP_STRING([--with-pspio-prefix=DIR], [Directory where pspio was installed.])],[],[with_pspio_prefix=$PSPIO_PREFIX])
case $with_pspio_prefix in
  no ) acx_pspio_ok=disable ;;
  *) LIBS_PSPIO="-L$with_pspio_prefix/lib -lpspio"; FCFLAGS_PSPIO="$ax_cv_f90_modflag$with_pspio_prefix/include" ;;
esac

dnl The tests
if test "$acx_pspio_ok" = no; then
  AC_MSG_CHECKING([for pspio])
  # If the location has been passed with --with-pspio-prefix just test this
  if test "$LIBS_PSPIO"; then
    pspio_fcflags="$FCFLAGS_PSPIO"; pspio_libs="$LIBS_PSPIO"
    FCFLAGS="$pspio_fcflags $acx_pspio_save_FCFLAGS $GSL_CFLAGS"
    LIBS="$pspio_libs $acx_pspio_save_LIBS $GSL_LIBS"
    AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    use pspio_f90_types_m
    use pspio_f90_lib_m

    type(pspio_f90_pspdata_t) :: pspdata
    integer :: i
    i = pspio_f90_pspdata_free(pspdata)
]), [acx_pspio_ok=yes; FCFLAGS_PSPIO="$pspio_fcflags"; LIBS_PSPIO="$pspio_libs"], [])
  else
    pspio_libs="-lpspio"
    FCFLAGS="$pspio_fcflags $acx_pspio_save_FCFLAGS $GSL_CFLAGS"
    LIBS=" $acx_pspio_save_LIBS $pspio_libs $GSL_LIBS"
    AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    use pspio_f90_types_m
    use pspio_f90_lib_m

    type(pspio_f90_pspdata_t) :: pspdata
    integer :: i
    i = pspio_f90_pspdata_free(pspdata)
]), [acx_pspio_ok=yes; FCFLAGS_PSPIO="$pspio_fcflags"; LIBS_PSPIO="$pspio_libs"], [])
  fi
  AC_MSG_RESULT([$acx_pspio_ok ($LIBS_PSPIO)])
fi

AC_SUBST(FCFLAGS_PSPIO)
AC_SUBST(LIBS_PSPIO)
FCFLAGS="$acx_pspio_save_FCFLAGS"
LIBS="$acx_pspio_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pspio_ok" = xyes; then
  AC_DEFINE(HAVE_PSPIO,1,[Defined if you have the PSPIO library.])
  $1
else
  AC_DEFINE(HAVE_PSPIO,0,[Defined if you have the PSPIO library.])
  LIBS_PSPIO=""
  FCFLAGS_PSPIO=""
  $2
fi
])dnl ACX_PSPIO
