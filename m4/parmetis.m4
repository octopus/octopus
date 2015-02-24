## Copyright (C) 2015 David Strubbe
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

AC_DEFUN([ACX_PATH_PARMETIS],[

AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([ACX_PATH_METIS])
acx_parmetis_ok=yes
AC_MSG_CHECKING(for ParMETIS library)

if test x"$acx_mpi_ok" != xyes; then
  AC_MSG_RESULT([requires MPI])
  acx_parmetis_ok=no
else
  if test x"$acx_external_metis" != xyes; then
    # To compile ParMETIS, you needed external METIS. So, no point in depending on internal METIS here.
    AC_MSG_RESULT([requires external METIS])
    acx_parmetis_ok=no
  fi
fi

if test x"$acx_parmetis_ok" != xno; then
  AC_ARG_WITH([parmetis-prefix],
    [AC_HELP_STRING([--with-parmetis-prefix],
    [Directory where ParMETIS library was installed])])

    # FIXME: accept setting by environment variables too
  case $with_parmetis_prefix in
    no ) acx_parmetis_ok=disabled ;;
    "") with_parmetis_prefix="/usr/include" ;;
  esac

  if test x"$acx_parmetis_ok" != xdisabled; then
  
    dnl Backup LIBS and FCFLAGS
    acx_parmetis_save_CFLAGS="$CFLAGS"
    acx_parmetis_save_LIBS="$LIBS"

    PARMETIS_CFLAGS="-I$with_parmetis_prefix/include"
    LIBS_PARMETIS="-L$with_parmetis_prefix/lib -lparmetis"

    CFLAGS="$CFLAGS $PARMETIS_CFLAGS $METIS_CFLAGS"
    LIBS="$LIBS $LIBS_PARMETIS $LIBS_METIS"

    AC_LANG_SAVE
    AC_LANG_C

    AC_LINK_IFELSE([AC_LANG_PROGRAM([
#include <parmetis.h>
#include <stdlib.h>
    ],[
idx_t *idx;
real_t *real;
MPI_Fint *fcomm;
ParMETIS_V3_PartKway(idx, idx, idx, NULL, NULL, idx, idx, idx, idx, real, real, idx, idx, idx, fcomm);
    ])], [acx_parmetis_ok=yes], [acx_parmetis_ok=no])

    AC_LANG_RESTORE
    AC_MSG_RESULT([$acx_parmetis_ok ($PARMETIS_CFLAGS $LIBS_PARMETIS)])

    CFLAGS="$acx_parmetis_save_CFLAGS"
    LIBS="$acx_parmetis_save_LIBS"
  fi

  if test x"$acx_parmetis_ok" = xyes; then
    HAVE_PARMETIS=1
    AC_DEFINE(HAVE_PARMETIS, 1, [Defined if you have the PARMETIS library.])
  else
    AC_MSG_WARN(Octopus will be compiled without PARMETIS support)
    PARMETIS_CFLAGS=""
    LIBS_PARMETIS=""
  fi

  AC_SUBST(PARMETIS_CFLAGS)
  AC_SUBST(LIBS_PARMETIS)
fi
])
