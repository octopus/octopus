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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id$
##

dnl NOT available from the GNU Autoconf Macro Archive at:
dnl https://trac.version.fz-juelich.de/libfm
dnl
AC_DEFUN([ACX_LIBFM], [
acx_libfm_ok=no

dnl We cannot use LIBFM if MPI is not found
if test "x$acx_mpi_ok" != xyes; then
  acx_libfm_ok=nompi
fi

dnl Get fortran linker name of LIBFM function to check for.
dnl if not compiling with fortran, convert the names
m4_if(_AC_LANG, Fortran, [fmm_func=fcs_init], [AC_FC_FUNC(fmm_func)])

dnl Check if the library was given in the command line
if test $acx_libfm_ok = no; then
  AC_ARG_WITH(libfm, [AS_HELP_STRING([--with-libfm=<lib>], [use LIBFM library <lib>])])
  case $with_libfm in
    yes | "") ;;
    no) acx_libfm_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_LIBFM="$with_libfm" ;;
    *) LIBS_LIBFM="-l$with_libfm" ;;
  esac
fi

dnl Backup LIBS 
acx_libfm_save_LIBS="$LIBS"
LIBS="$FLIBS $LIBS_LIBFM $LIBS"

echo $LIBS_LIBFM
dnl First, check LIBS_LIBFM environment variable
if test $acx_libfm_ok = no; then
  AC_MSG_CHECKING([for $fmm_func in $LIBS_LIBFM])
dnl AC_TRY_LINK_FUNC($fmm_func, [acx_libfm_ok=yes], [])
dnl if test $acx_libfm_ok = no; then
dnl   AC_MSG_RESULT([$acx_libfm_ok ($LIBS_LIBFM)])
dnl else
dnl   AC_MSG_RESULT([$acx_libfm_ok ($LIBS_LIBFM)])
dnl fi
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
  #include <fcs_fconfig.h>
    use fcs_module
    use iso_fortran_env
    use iso_c_binding
    implicit none
    type(c_ptr)                                     ::  handle
    type(c_ptr)                                     ::  ret
    character(len = 8)                              ::  method
    integer                                         ::  communicator
    ret = fcs_init(handle, trim(adjustl(method)) // c_null_char, communicator)
  ]), [acx_libfm_ok=yes], [])
fi

dnl Generic LIBFM library?
for libfm in fm_r64; do
  if test $acx_libfm_ok = no; then
    AC_CHECK_LIB($libfm , $fmm_func,
      [acx_libfm_ok=yes; LIBS_LIBFM="$LIBS_LIBFM -l$libfm"], [], [$FLIBS])
  fi
done

dnl Generic LIBFM library?
for libfm in fm_r64; do
  dnl if test x"$libfm" = xlibfm-openmpi; then       
  dnl   libfmCinit="libfmCinit-openmpi"
  dnl else
  dnl libfm="fm_r64"
  dnl fi
  if test $acx_libfm_ok = no; then
    AC_CHECK_LIB($libfm, $fmm_func,
      [acx_libfm_ok=yes; LIBS_LIBFM="$LIBS_LIBFM -l$libfm"], [], [$FLIBS])
  fi
done

AC_SUBST(LIBS_LIBFM)
LIBS="$acx_libfm_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libfm_ok" = xyes; then
  AC_DEFINE(HAVE_LIBFM,1,[Defined if you have LIBFM library.])
  $1
else
  AC_MSG_WARN([Could not find Libfm library. 
               *** Will compile without Libfm support])
  $2
fi
])dnl ACX_LIBFM
