## Copyright (C) 2008 X. Andrade
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
## $Id: papi.m4 4621 2008-10-08 17:21:03Z xavier $
##
################################################

################################################################
# Check whether papi is installed
# ----------------------------------
AC_DEFUN([ACX_PAPI],
[
acx_papi_header=no
acx_papi_func=no
acx_papi=no

AC_CHECK_HEADER(papi.h, acx_papi_header=yes)
AC_CHECK_FUNC(PAPI_library_init, acx_papi_func=yes)

if test x"${acx_papi_header}" = x"yes" -a x"${acx_papi_func}" = x"yes"; then
acx_papi=yes
AC_DEFINE(HAVE_PAPI, 1, [PAPI is available])
AC_MSG_NOTICE([Octopus will use PAPI for performance measuring])
fi

AM_CONDITIONAL(COMPILE_PAPI, test x$acx_papi = xyes)

])
