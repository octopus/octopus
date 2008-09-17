dnl Taken from Graphviz - Graph Visualization Software
dnl http://www.graphviz.org/
dnl Copyright: Common Public License Version 1.0

dnl -----------------------------------
dnl INCLUDES and LIBS for gd

AC_DEFUN([ACX_GDLIB],
[
  dnl We enable the GD library by default
  acx_gdlib_ok=yes

  dnl We disable GD support only if the user is requesting this explicitely
  AC_ARG_ENABLE(gdlib, AS_HELP_STRING([--disable-gdlib], [Do not compile with internal GD library support]), [acx_gdlib_ok=no])

  if test x"$acx_gdlib_ok" = xyes; then
    AC_PATH_PROG(GDLIB_CONFIG, gdlib-config)
    if test -n "$GDLIB_CONFIG"; then
      AC_DEFINE_UNQUOTED(HAVE_GDLIB, 1, [Define if libgd exists.])

      GD_CFLAGS=`$GDLIB_CONFIG --cflags`
      GD_LIBS="`$GDLIB_CONFIG --ldflags` -lgd `$GDLIB_CONFIG --libs`"
      GD_MAJORVERSION=`$GDLIB_CONFIG --majorversion`
      GD_MINORVERSION=`$GDLIB_CONFIG --minorversion`
      GD_REVISION=`$GDLIB_CONFIG --revision`
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

      # GD version check >= 2.0.33
dnl      gd_version_ok=yes
dnl      REQ_GD_MAJORVERSION=2
dnl      REQ_GD_MINORVERSION=0
dnl      REQ_GD_REVISION=33
dnl      if test $GD_MAJORVERSION -lt $REQ_GD_MAJORVERSION; then
dnl        gd_version_ok=no
dnl      else
dnl        if test $GD_MAJORVERSION -eq $REQ_GD_MAJORVERSION; then
dnl          if test $GD_MINORVERSION -lt $REQ_GD_MINORVERSION; then
dnl            gd_version_ok=no
dnl          else
dnl            if test $GD_MINORVERSION -eq $REQ_GD_MINORVERSION; then
dnl              if test $GD_REVISION -lt $REQ_GD_REVISION; then
dnl                gd_version_ok=no
dnl              fi
dnl            fi
dnl          fi
dnl        fi
dnl      fi
dnl      if test "x$gd_version_ok" = "xno"; then
dnl        AC_MSG_ERROR(GD library version < $REQ_GD_MAJORVERSION.$REQ_GD_MINORVERSION.$REQ_GD_REVISION)
dnl      fi

    else
      AC_MSG_WARN([GD library support has been disabled.
                    *** some esoteric parts of octopus will not work.])
    fi
  else
    AC_MSG_WARN([GD library support has been disabled.
                  *** some esoteric parts of octopus will not work.])
  fi

  AC_SUBST(GD_CFLAGS)
  AC_SUBST(GD_LIBS)

])
