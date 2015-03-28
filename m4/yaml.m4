## Copyright (C) YAML developers
##
## $Id$	
##

AC_DEFUN([ACX_YAML], [
  dnl Yaml input file support
  AC_ARG_WITH([yaml-path],
              AS_HELP_STRING([--with-yaml-path], [give a path to find libyaml.]),
              [ac_path_yaml=$withval])
  AC_ARG_ENABLE(internal-libyaml, AS_HELP_STRING([--disable-internal-libyaml], [Do not build and link with internal libyaml library (default = auto).]), ac_build_libyaml=$enableval, ac_build_libyaml="auto")
  
  ac_use_libyaml="no"
  if test x"$ac_path_yaml" == x"" ; then
     ac_path_yaml="/usr"
  fi
  dnl auto-detect if build = no or build = auto
  if test x"$ac_build_libyaml" != x"yes" ; then
     LDFLAGS_SVG="$LDFLAGS"
     AC_LANG_PUSH(C)
     LDFLAGS="$LDFLAGS -L$ac_path_yaml/lib"
     AC_CHECK_LIB([yaml], [yaml_parser_parse],
                  [ac_use_libyaml=yes], [ac_use_libyaml=no])
     if test x"$ac_use_libyaml" = x"yes"; then
        if test x"$ac_path_yaml" != x"/usr" ; then
           LIB_YAML_CFLAGS="-I$ac_path_yaml/include"
           LIB_YAML_LIBS="-L$ac_path_yaml/lib "
        fi
        LIB_YAML_LIBS=$LIB_YAML_LIBS"-lyaml"
     else
        AC_MSG_WARN([libyaml is not available, building internal one.])
     fi
     AC_LANG_POP(C)
     LDFLAGS="$LDFLAGS_SVG"
  fi
  dnl internal if yes or auto and not detected.
  if test x"$ac_use_libyaml" != x"yes" ; then
     ac_use_libyaml="yes"
     ac_build_libyaml="yes"
     LIB_YAML_CFLAGS="-I\$(top_srcdir)"/external_libs/yaml-0.1.4/include
     LIB_YAML_LIBS="\$(top_builddir)/external_libs/yaml-0.1.4/src/.libs/libyaml.a"
     dnl tar -xzf ${srcdir}/external_libs/PyYAML-3.10.tar.gz
  fi
  AC_DEFINE([HAVE_YAML], [], [If set, we can call yaml.h])
  AM_CONDITIONAL(BUILD_LIBYAML, test x"$ac_build_libyaml" = x"yes")
  AM_CONDITIONAL(HAVE_LIB_YAML, test x"$ac_use_libyaml" = x"yes")
  AC_SUBST(LIB_YAML_CFLAGS)
  AC_SUBST(LIB_YAML_LIBS)
])dnl ACX_YAML
