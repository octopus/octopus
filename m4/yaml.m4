## Copyright (C) YAML developers
##
## $Id$	
##

AC_DEFUN([ACX_YAML], [
  AC_MSG_CHECKING(for YAML library)


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
     LDFLAGS="-L$ac_path_yaml/lib"
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
  # Define the package version numbers and the bug reporting link of yaml.
  m4_define([YAML_MAJOR], 0)
  m4_define([YAML_MINOR], 1)
  m4_define([YAML_PATCH], 4)
  m4_define([YAML_BUGS], [http://pyyaml.org/newticket?component=libyaml])
  
  m4_define([YAML_RELEASE], 0)
  m4_define([YAML_CURRENT], 2)
  m4_define([YAML_REVISION], 2)
  m4_define([YAML_AGE], 0)
  
  # Define macro variables for the package version numbers of yaml.
  AC_DEFINE(YAML_VERSION_MAJOR, YAML_MAJOR, [Define the major version number.])
  AC_DEFINE(YAML_VERSION_MINOR, YAML_MINOR, [Define the minor version number.])
  AC_DEFINE(YAML_VERSION_PATCH, YAML_PATCH, [Define the patch version number.])
  AC_DEFINE(YAML_VERSION_STRING, "YAML_MAJOR.YAML_MINOR.YAML_PATCH", [Define the version string.])
])dnl ACX_YAML
