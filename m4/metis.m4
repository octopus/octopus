## Copyright (C) DUNE project. http://users.dune-project.org/
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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$
##

dnl @synopsis IMMDX_LIB_METIS([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl
dnl This macro searches for the METIS library in the user specified
dnl location. The user may specify the location either by defining the
dnl environment variable METIS or by using the --with-metis option to
dnl configure. If the environment variable is defined it has precedent
dnl over everything else. If no location was specified then it searches
dnl in /usr/lib and /usr/local/lib for the library and in /usr/include
dnl and /usr/local/include for the header files. Upon sucessful
dnl completion the variables METIS_LIB and METIS_INCLUDE are set.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a METIS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_METIS. If ACTION-IF-NOT-FOUND is
dnl not specified then an error will be generated halting configure.
dnl
dnl @category InstalledPackages
dnl @author Ben Bergen <ben@cs.fau.de>
dnl @version 2003-01-19
dnl @license AllPermissive

AC_DEFUN([ACX_PATH_METIS], [
  AC_MSG_CHECKING(for METIS library)
  AC_REQUIRE([AC_PROG_CC])
  #
  # User hints...
  #
  AC_ARG_WITH([metis],
    [AC_HELP_STRING([--with-metis],
    [user defined path to METIS library])],
    [
      if test -n "$METIS" ; then
        AC_MSG_RESULT(yes)
	with_metis=$METIS
      elif test "$withval" != no ; then
	AC_MSG_RESULT(yes)
	with_metis=$withval
      else
	AC_MSG_RESULT(no)
      fi
    ],
    [
      if test -n "$METIS" ; then
	with_metis=$METIS
	AC_MSG_RESULT(yes)
      else
	with_metis=/usr
	if test ! -f "$with_metis/include/metis.h" ; then
	  if test ! -f "$with_metis/include/metis/metis.h" ; then 
            with_metis=/usr/local
            if test ! -f "$with_metis/include/metis.h" ; then
              with_metis=""
              AC_MSG_RESULT(no)
            else
              AC_MSG_RESULT(yes)
            fi
	  else
	    AC_MSG_RESULT(yes)
	  fi
	else
	  AC_MSG_RESULT(yes)
	fi
      fi
    ])
#
# locate METIS library
#
    if test -n "$with_metis" ; then
      old_CFLAGS=$CFLAGS
      old_LDFLAGS=$LDFLAGS

      if test -f "$with_metis/include/metis.h"; then
	lib_path="/lib"
	include_path="/include"
      fi

      if test -f "$with_metis/include/metis/metis.h"; then
	lib_path="/lib"
	include_path="/include/metis"
      fi

      if test -f "$with_metis/Lib/metis.h"; then
        # catch bad convention in the downloadable metis version
	lib_path=""
	include_path="/Lib"
      fi
      CFLAGS="-I$with_metis/$include_path"
      LDFLAGS="-L$with_metis/$lib_path"

      AC_LANG_SAVE
      AC_LANG_C

      AC_CHECK_LIB(metis, METIS_PartMeshDual,
			[metis_lib=yes], [metis_lib=no], [-lm])
			

      AC_CHECK_HEADER(metis.h, [metis_h=yes],
			[metis_h=no], [/* check */])

      AC_LANG_RESTORE

      CFLAGS=$old_CFLAGS 
      LDFLAGS=$old_LDFLAGS

      AC_MSG_CHECKING(METIS in $with_metis)
      if test "$metis_lib" = "yes" -a "$metis_h" = "yes" ; then
	AC_SUBST(METIS_INCLUDE, [-I$with_metis$include_path])
	AC_SUBST(METIS_LDFLAGS, [-L$with_metis$lib_path])
	AC_SUBST(METIS_LIB, [-lmetis])
	AC_MSG_RESULT(ok)
	LIBS_METIS_5="$METIS_LDFLAGS $METIS_LIB"
      else
	AC_MSG_RESULT(no)
	LIBS_METIS_5=""
      fi
    fi
		
    # tell automake
    AM_CONDITIONAL(METIS, test x$METIS_LIB = x1)
    if test x = x"$METIS_LIB" ; then
      with_metis=no
            
      dnl METIS was not found to link with, but is included in the distribution

      dnl We disable METIS support only if the user is requesting this explicitly
      AC_ARG_ENABLE(metis, AS_HELP_STRING([--disable-metis], [Do not compile with internal METIS domain-partitioning library.]),[acx_metis_ok=$enableval],[acx_metis_ok=yes])
      dnl Since v2.0, METIS ships with octopus, so we enable it by default

      AC_MSG_CHECKING([whether METIS included in Octopus is enabled])

      AC_MSG_RESULT([$acx_metis_ok])

      if test x"$acx_metis_ok" = xyes; then
        HAVE_METIS=1
        HAVE_COMP_METIS=1
        AC_DEFINE(HAVE_METIS, 1, [This is defined when we should compile with METIS support (default).])
        AC_DEFINE(HAVE_COMP_METIS, 1, [This is defined when we link with an external METIS library (default).])
      else
        AC_MSG_WARN(Octopus will be compiled without METIS support)
      fi
    else
      with_metis=yes
      AC_DEFINE(HAVE_METIS,1,[Define if you have METIS library])
    fi

    AC_SUBST(LIBS_METIS_5)
])dnl ACX_PATH_METIS
