## Copyright (C) 2010 The octopus team
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
## $Id: fc_command_line.m4 5986 2009-10-26 06:27:57Z dstrubbe $
##

################################################
# Generates macro definitions for:
#  LATEST_SVN : the svn revision number.
#  BUILD_TIME : when the configure script is launched.
#  FC : The "true" Fortran compiler (figuring out which compiler hides under the mpif90 disguise).
#  CC : The "true" C compiler (figuring out which compiler hides under the mpicc disguise).
#  CFLAGS : The flags passed to the C compiler.
#  FCFLAGS : The flags passed to the Fortran compiler.
# ----------------------------------
AC_DEFUN([ACX_OCTOPUS_COMPILATION_INFO],
[
AC_MSG_NOTICE([collecting compilation info...])

rev=$(sh ./build/svn_release_number.sh 2> /dev/null)
date=`date`
truecc=$(sh ./build/true_compiler.sh $CC 2> /dev/null)
truefc=$(sh ./build/true_compiler.sh $FC 2> /dev/null)

AC_DEFINE_UNQUOTED([LATEST_SVN], ["$rev"], [subversion revision number])
AC_DEFINE_UNQUOTED([BUILD_TIME], ["$date"], [date when configure was launched])
AC_DEFINE_UNQUOTED([CC], ["$CC $truecc"], [C compiler])
AC_DEFINE_UNQUOTED([FC], ["$FC $truefc"], [Fortran compiler])
# 132 characters is max in std fortran, and two are for the quotes
AC_DEFINE_UNQUOTED([CFLAGS], ["${CFLAGS:0:130}"], [C compiler])
AC_DEFINE_UNQUOTED([FCFLAGS], ["${FCFLAGS:0:130}"], [Fortran compiler])

]
)
