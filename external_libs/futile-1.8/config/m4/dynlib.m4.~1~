# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

# AX_DYNAMIC_LIBRARIES([DEFAULT = "no"])
AC_DEFUN([AX_DYNAMIC_LIBRARIES],
[dnl Produce dynamic libraries and executables.
  AC_ARG_ENABLE(dynamic-libraries, AS_HELP_STRING([--enable-dynamic-libraries],
                                                 [Build dynamical libraries (disabled by default).]),
                [ax_build_dynamic=$enableval; ax_user_input="yes"],
                [ax_build_dynamic=m4_default([$1], ["no"]); ax_user_input="no"])
  dnl Test for library building tools.
  if test x"$ax_build_dynamic" = x"yes" ; then
    AC_REQUIRE([AX_FLAG_PIC])
    if test -z "$ax_flag_pic" ; then
      AC_MSG_WARN(["No available position-independent code flag, dynamic libraries disabled."])
      ax_build_dynamic=no
    else
      dnl try to link with dependencies to see if its working.
      dnl dependencies may not have been compiled with PIC.
      ax_dynamic_deps=yes
      ax_dyndeps_libs=m4_default([$3], [""])
      if test -n "$ax_dyndeps_libs" -a -n "$2" ; then
        AC_MSG_CHECKING([for dynamic linking of dependencies])
        AC_LANG_PUSH(Fortran)
        LDFLAGS_SVG=$LDFLAGS
        LDFLAGS="-shared $LDFLAGS"
        LIBS_SVG=$LIBS
        LIBS="$ax_dyndeps_libs $LIBS"
        FCFLAGS_SVG="$FCFLAGS"
        FCFLAGS="$ax_flag_pic $FCFLAGS"
        AC_LINK_IFELSE([[
subroutine testlib()
  call $2
end subroutine testlib
        ]], ax_dynamic_deps=yes, ax_dynamic_deps=no)
        AC_LANG_POP([Fortran])
        FCFLAGS=$FCFLAGS_SVG
        LIBS=$LIBS_SVG
        LDFLAGS=$LDFLAGS_SVG
        AC_MSG_RESULT([$ax_dynamic_deps])
      fi

      if test "$ax_dynamic_deps" = "yes" ; then
        case $CFLAGS   in *"$ax_flag_pic"*) ;; *) CFLAGS="$CFLAGS $ax_flag_pic";; esac
        case $CXXFLAGS in *"$ax_flag_pic"*) ;; *) CXXFLAGS="$CXXFLAGS $ax_flag_pic";; esac
        case $FCFLAGS  in *"$ax_flag_pic"*) ;; *) FCFLAGS="$FCFLAGS $ax_flag_pic";; esac
    
        eval "set x $ac_configure_args"
        shift
        ac_configure_args=
        for ac_arg ; do
          case $ac_arg in
            CFLAGS=* | CXXFLAGS=* | FCFLAGS=*) ;;
            *) ac_configure_args="$ac_configure_args '$ac_arg'"
          esac
        done
        ac_configure_args="$ac_configure_args 'CFLAGS=$CFLAGS' 'CXXFLAGS=$CXXFLAGS' 'FCFLAGS=$FCFLAGS'"
      else
        ax_build_dynamic="no"
      fi
    fi

    dnl Raise error if asked for dynamic but cannot proceed.
    if test "$ax_user_input" == "yes" -a "$ax_build_dynamic" != "yes" ; then
      AC_MSG_ERROR(["Dynamic build is not possible."])
    fi
  else
    ax_build_dynamic="no"
  fi

  AM_CONDITIONAL([BUILD_DYNAMIC_LIBS], [test "x$ax_build_dynamic" = "xyes"])
])
