## Copyright (C) 2015 D. Strubbe
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

AC_DEFUN([ACX_PATH_METIS], [

acx_external_metis=no
AC_MSG_CHECKING(for external METIS library)

# METIS is only useful in parallel.
if test x"$acx_mpi_ok" != xyes; then
  AC_MSG_RESULT([not used without MPI])
else
  acx_external_metis=no
  AC_REQUIRE([AC_PROG_CC])

  AC_ARG_WITH([metis-prefix],
    [AS_HELP_STRING([--with-metis-prefix=DIR],
    [Directory where external METIS library was installed (must be single-precision)])])

  case $with_metis_prefix in
    yes | "") ;;
    no ) acx_external_metis=disabled ;;
    *) LIBS_METIS="-L$with_metis_prefix/lib -lmetis";
       METIS_CFLAGS="-I$with_metis_prefix/include"
  esac

  if test x"$acx_external_metis" != xdisabled; then
  
    dnl Backup CFLAGS and LIBS
    acx_metis_save_CFLAGS="$CFLAGS"
    acx_metis_save_LIBS="$LIBS"

    if test "x${LIBS_METIS+set}" != xset ; then
      LIBS_METIS="-lmetis"
    fi

    CFLAGS="$CFLAGS $METIS_CFLAGS"
    LIBS="$LIBS $LIBS_METIS"

    AC_LANG_SAVE
    AC_LANG_C

    AC_LINK_IFELSE([AC_LANG_PROGRAM([
#include <metis.h>
#if defined(METIS_USE_DOUBLEPRECISION) || REALTYPEWIDTH == 64
  #error METIS must be compiled in single precision for Octopus.
#endif
],[
idx_t *options;
METIS_SetDefaultOptions(options);
    ])], [acx_external_metis=yes], [])

    AC_LANG_RESTORE
    AC_MSG_RESULT([$acx_external_metis ($METIS_CFLAGS $LIBS_METIS)])

    CFLAGS="$acx_metis_save_CFLAGS"
    LIBS="$acx_metis_save_LIBS"
  else
    AC_MSG_RESULT([disabled])
    acx_external_metis=no
  fi

  if test x"$acx_external_metis" = xno ; then
    dnl METIS was not found to link with, but is included in the distribution
  
    dnl We disable METIS support only if the user is requesting this explicitly
    AC_ARG_ENABLE(metis, AS_HELP_STRING([--disable-metis],
    			 [Do not compile with internal METIS domain-partitioning library.]),
			 [acx_internal_metis=$enableval],[acx_internal_metis=yes])
  
    AC_MSG_CHECKING([whether METIS included in Octopus is enabled])
  
    AC_MSG_RESULT([$acx_internal_metis])
  
    if test x"$acx_internal_metis" = xyes; then
      HAVE_METIS=1
      HAVE_COMP_METIS=1
      AC_DEFINE(HAVE_METIS, 1, [This is defined when we should compile with METIS support (default).])
      AC_DEFINE(HAVE_COMP_METIS, 1, [This is defined when we link with the internal METIS library (default).])
    else
      AC_MSG_WARN(Octopus will be compiled without METIS support)
    fi
  else
    acx_internal_metis=no
    AC_DEFINE(HAVE_METIS,1,[This is defined when we should compile with METIS support (default).])
  fi
fi

if test x"$acx_external_metis" = xno; then
  METIS_CFLAGS=""
  LIBS_METIS=""
fi

AC_SUBST(METIS_CFLAGS)
AC_SUBST(LIBS_METIS)

])dnl ACX_PATH_METIS
