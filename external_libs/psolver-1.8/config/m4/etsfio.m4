# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_ETSF_IO],
[dnl Use ETSF_IO
  ax_have_etsf_io=no
  AC_ARG_WITH(etsf-io, AS_HELP_STRING([--with-etsf-io], [Link with ETSF_IO library (default = no).]), etsfio=$withval, etsfio="yes")
  AC_ARG_WITH(etsf-io-path, AS_HELP_STRING([--with-etsf-io-path], [Give the path of the ETSF_IO installation (default = /usr).]),
              ac_etsfio_dir=$withval, ac_etsfio_dir=)
  AC_ARG_WITH(netcdf-path, AS_HELP_STRING([--with-netcdf-path], [Give the path to NetCDF (required by ETSF_IO) (default = /usr).]),
              ac_netcdf_dir=$withval, ac_netcdf_dir=)
  AC_ARG_WITH(netcdf-libs, AS_HELP_STRING([--with-netcdf-libs], [Give the library to link with NetCDF (required by ETSF_IO) (default = -lnetcdff -lnetcdf).]),
              ac_netcdf_libs=$withval, ac_netcdf_libs="-lnetcdff -lnetcdf")
  if test "$etsfio" = "yes" ; then
     LDFLAGS_SVG="$LDFLAGS"
     LIBS_SVG="$LIBS"
     FCFLAGS_SVG="$FCFLAGS"
     if test -n "$ac_etsfio_dir" ; then
        LDFLAGS="$LDFLAGS -L$ac_etsfio_dir/lib"
        ac_etsfio_incs="-I$ac_etsfio_dir/include"
        FCFLAGS="$FCFLAGS $ac_etsfio_incs"
     elif test -n "$C_INCLUDE_PATH"; then
        for path in ${C_INCLUDE_PATH//:/ }; do
           ac_etsfio_incs="$ac_etsfio_incs -I$path"
        done
        FCFLAGS="$FCFLAGS $ac_etsfio_incs"
     fi
     if test -n "$ac_netcdf_dir" -a x"$ac_netcdf_dir" != x"$ac_etsfio_dir" ; then
        LDFLAGS="$LDFLAGS -L$ac_netcdf_dir/lib"
        FCFLAGS="$FCFLAGS -I$ac_netcdf_dir/include"
     fi
     LIBS="$LIBS -letsf_io $ac_netcdf_libs"
     AC_MSG_CHECKING([for ETSF_IO library])
     AC_LINK_IFELSE([[
  program main
    use etsf_io
    
    type(etsf_groups_flags) :: groups
    type(etsf_dims) :: dims
    logical :: lstat
    type(etsf_io_low_error) :: error_data
  
    call etsf_io_data_init("test", groups, dims, "test", "", lstat, error_data)
  end]], ax_have_etsf_io=yes, ax_have_etsf_io=no)
     AC_MSG_RESULT([$ax_have_etsf_io])
     LIBS="$LIBS_SVG"
     FCFLAGS="$FCFLAGS_SVG"
     if test "$ax_have_etsf_io" = "yes"; then
        if test -n "$ac_netcdf_dir" -a x"$ac_netcdf_dir" != x"$ac_etsfio_dir" ; then
           ac_etsfio_incs=$ac_etsfio_incs" -I$ac_netcdf_dir/include"
        fi
        AC_SUBST(LIBETSFIO_INCLUDE, $ac_etsfio_incs)
        LIBETSFIO_LIBS="-letsf_io_utils -letsf_io $ac_netcdf_libs"
     else
        ax_have_etsf_io="warn"
     fi
  fi
  AM_CONDITIONAL(HAVE_ETSF_IO, test "$ax_have_etsf_io" = "yes")
])