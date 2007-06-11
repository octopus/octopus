## Copyright (C) 2007 X. Andrade
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

#
# Check how to access command line arguments
# ----------------------------------

AC_DEFUN([ACX_FC_COMMAND_LINE_ARGUMENTS],
[
AC_MSG_CHECKING([how to access command line arguments])

#try Fortran 2003 API first

AC_LINK_IFELSE(
    AC_LANG_PROGRAM( [], [
    implicit none
    integer :: i
    character(len=32) :: arg
    i = command_argument_count()
    call get_command_argument(0, arg)
    ]),
    [acx_command_line_arguments="fortran 2003"
     AC_DEFINE(FC_COMMAND_LINE_ARGUMENTS, 2003, [compiler supports command line arguments])],
    [acx_command_line_arguments=none])

if test "$acx_command_line_arguments" == "none" ; then

#Fortran 77

#some compilers require a module or include
#                         NAS        NAG          DEC?   PGI
acx_command_line_modules="nas_system f90_unix_env dfport lib3f.h"

for acx_command_line_module in "fortran 77" $acx_command_line_modules ; do

if test "$acx_command_line_module" != "fortran 77" ; then
acx_cl_usemodule="use "$acx_command_line_module
fi

if test "$acx_command_line_module" == "lib3f.h" ; then
acx_cl_usemodule=""
acx_cl_includeheader="include 'lib3f.h'"
fi


AC_LINK_IFELSE(
    AC_LANG_PROGRAM( [], [
    $acx_cl_usemodule
    implicit none
    $acx_cl_includeheader

    integer :: i
    character(len=32) :: arg
    i = iargc()
    call getarg(0, arg)
    ]),
    [
    acx_command_line_arguments="$acx_command_line_module"
    AC_DEFINE(FC_COMMAND_LINE_ARGUMENTS, 77, [compiler supports command line arguments])
    if test "$acx_command_line_module" != "fortran 77" ; then
        if test "$acx_command_line_module" != "lib3f.h" ; then
	AC_DEFINE_UNQUOTED(FC_COMMAND_LINE_MODULE, $acx_command_line_module, [module required to get command line arguments])
	else
	AC_DEFINE_UNQUOTED(FC_COMMAND_LINE_INCLUDE, "$acx_command_line_module", [include required to get command line arguments])
	fi
    fi
    ]
    ,[])
done
fi

if test "$acx_command_line_arguments" == "none" ; then

#last resource, external
AC_LINK_IFELSE(
    AC_LANG_PROGRAM( [], [

    implicit none 

    interface iargc 
    integer function iargc() 
    end function iargc 
    end interface 

    interface getarg 
    subroutine getarg(c, a)
    integer :: c
    character(len=*) :: a
    end subroutine getarg
    end interface
 
    integer :: i 
    character(len=32) :: arg 
    i = iargc() 
    call getarg(0, arg) 
    ]),
    [
    acx_command_line_arguments=external
    AC_DEFINE(FC_COMMAND_LINE_ARGUMENTS, 77, [compiler supports command line arguments])
    AC_DEFINE(FC_COMMAND_LINE_EXTERNAL, 1, [iargc and getarg are external and must be declared])
    ]
    ,[])
fi

AC_MSG_RESULT($acx_command_line_arguments)

if test "$acx_command_line_arguments" == "none" ; then
    AC_MSG_WARN([Could not find how to access command line from Fortran.
                *** Some utilities will not work.])
fi

]
)
