## Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
## $Id: configure.ac 2402 2006-09-14 16:59:22Z acastro $

AC_INIT([string_f],[1.0],[octopus-devel@tddft.org],[string_f])
AC_CONFIG_SRCDIR([string_f.h])
AM_INIT_AUTOMAKE

# Installation prefix by default
AC_PREFIX_DEFAULT([/opt/etfs])

# who am i
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for header files.
AC_HEADER_STDC
AC_FUNC_ALLOCA

# Checks for functions
AC_CHECK_FUNCS([malloc])

AC_CONFIG_FILES([Makefile]) 
AC_OUTPUT
