## Copyright (C) 2006 The octopus developer team
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

AC_PREREQ(2.59)
AC_INIT([Octopus GUI], [0.0.1], [octopus-devel@tddft.org], octopus-gui)
AM_INIT_AUTOMAKE

dnl AM_MAINTAINER_MODE

AM_PATH_PYTHON(2.4)

dnl should the gui be built? It defaults to no.
AC_ARG_ENABLE([gui],
	      AS_HELP_STRING([--enable-gui], [enable gui (default is NO)]),
	      [ac_cv_build_gui=$enableval],
	      [ac_cv_build_gui=no])

if test $ac_cv_build_gui = yes; then
  AC_CACHE_CHECK([whether to build the gui], [ac_cv_build_gui], [ac_cv_build_gui=no])
  PKG_PROG_PKG_CONFIG
  PKG_CHECK_MODULES(PyGTK,
    [gtk+-2.0 >= 2.6
     pygtk-2.0 >= 2.6	    
     gnome-python-2.0 >= 2.12.4
     bonobo-activation-2.0 >= 2.14.0])
fi

AC_CONFIG_FILES([Makefile var2xml.py octopus-gui])
AC_OUTPUT
