dnl Taken from Graphviz - Graph Visualization Software
dnl http://www.graphviz.org/
dnl Copyright: Common Public License Version 1.0

dnl -----------------------------------
dnl INCLUDES and LIBS for gd

AC_DEFUN([ACX_GDLIB],
[
  # We disable GD support only if the user is requesting this explicitly
  AC_ARG_ENABLE(gdlib, AS_HELP_STRING([--disable-gdlib], [Do not compile with GD image-processing library.]),[acx_gdlib_ok=$enableval],[acx_gdlib_ok=yes])
  # GD library is enabled by default

  if test x"$acx_gdlib_ok" = xyes; then
    AC_PATH_PROG(GDLIB_CONFIG, gdlib-config)

    if test -n "$GDLIB_CONFIG"; then
      if test "x$GD_CFLAGS" = x; then
	GD_CFLAGS=`$GDLIB_CONFIG --cflags`
      fi
       if test "x$GD_LIBS" = x; then
        GD_LIBS="-L`$GDLIB_CONFIG --libdir` -lgd `$GDLIB_CONFIG --ldflags` `$GDLIB_CONFIG --libs | awk '{if($NF=="@LIBICONV@"){$NF=""} print}'`"
      # Sometimes GD installation strangely leaves this token @LIBICONV@ in --libs, which must be removed
      fi
    else
      acx_gdlib_ok="no"
    fi
  fi

  if test x"$acx_gdlib_ok" = xyes; then
    acx_save_LIBS="$LIBS"
    LIBS="$LIBS $GD_LIBS"
    acx_save_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS $GD_CFLAGS"

    AC_MSG_CHECKING([whether gdlib can be linked])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
[#include <gd.h>],
[gdImagePtr im;]
)],[], [acx_gdlib_ok=no])
    AC_MSG_RESULT([$acx_gdlib_ok])

    LIBS="$acx_save_LIBS"
    if test x"$acx_gdlib_ok" = xno; then
      AC_MSG_WARN([GD library support has been disabled.
                   *** Some esoteric parts of octopus will not work.])
      CFLAGS="$acx_save_CFLAGS"
      GD_CFLAGS=""
      GD_LIBS=""
    else
      AC_DEFINE_UNQUOTED(HAVE_GDLIB, 1, [Define if libgd exists.])
      AC_SUBST(GD_CFLAGS)
      AC_SUBST(GD_LIBS)

      GD_MAJORVERSION=`$GDLIB_CONFIG --majorversion`
      GD_MINORVERSION=`$GDLIB_CONFIG --minorversion`
      GD_REVISION=`$GDLIB_CONFIG --revision`
      # we only use PNG, JPEG, GIF 
      for f in `$GDLIB_CONFIG --features` ; do
        case $f in
        GD_XPM )
          AC_DEFINE_UNQUOTED(HAVE_GD_XPM, 1, [Define if libgd supports xpm.])
          ;;
        GD_JPEG )
          AC_DEFINE_UNQUOTED(HAVE_GD_JPEG, 1, [Define if libgd supports jpeg.])
          ;;
        GD_FONTCONFIG )
          AC_DEFINE_UNQUOTED(HAVE_GD_FONTCONFIG, 1, [Define if libgd uses fontconfig.])
          ;;
        GD_FREETYPE )
          AC_DEFINE_UNQUOTED(HAVE_GD_FREETYPE, 1, [Define if libgd uses freetype.])
          ;;
        GD_PNG )
          AC_DEFINE_UNQUOTED(HAVE_GD_PNG, 1, [Define if libgd supports png.])
          ;;
        GD_GIF )
          AC_DEFINE_UNQUOTED(HAVE_GD_GIF, 1, [Define if libgd supports gif.])
          ;;
        esac
      done
    fi
  fi
])
