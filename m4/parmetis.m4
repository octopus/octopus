dnl -*- mode: autoconf; tab-width: 8; indent-tabs-mode: nil; -*-
dnl vi: set et ts=8 sw=2 sts=2:
# $Id$
# searches for ParMetis headers and libs

AC_DEFUN([ACX_PATH_PARMETIS],[
  AC_MSG_CHECKING(for ParMETIS library)
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PATH_XTRA])
dnl joseba  AC_REQUIRE([DUNE_MPI])

  #
  # USer hints ...
  #
  AC_ARG_VAR([PARMETIS], [ParMETIS library location])
  AC_ARG_WITH([parmetis],
    [AC_HELP_STRING([--with-parmetis],[user defined path to ParMETIS library])],
    [
      # --with-parmetis supersedes $PARMETIS
      PARMETIS=""
      if test "$withval" != no ; then
          if test -d "$withval" ; then
	    # get absolute path
	    with_parmetis=`eval cd $withval 2>&1 && pwd`
            if test -f "$with_parmetis/include/parmetis.h" ; then
              lib_path=lib
              include_path=include
            fi
	    AC_MSG_RESULT(yes)
          else
            with_parmetis="no"
	    AC_MSG_RESULT(no)
      fi
	else
	    AC_MSG_RESULT(no)
	fi
	],
    [
    #echo with_parmetis1=$withparmetis PARMETIS=$PARMETIS
	if test -n "$PARMETIS" ; then
          if test -d "$PARMETIS" ; then
	    # get absolute path
	    with_parmetis=`eval cd $PARMETIS 2>&1 && pwd`
            PARMETIS=""
	    AC_MSG_RESULT(yes)
          else
            PARMETIS=""
            with_parmetis=no
	    AC_MSG_RESULT(no)
          fi
	else
	    with_parmetis=/usr/
	    include_path=include
	    lib_path=lib
	    if test ! -f "$with_parmetis/$include_path/parmetis.h" ; then
		with_parmetis=/usr/local/
		if test ! -f "$with_metis/$include_path/parmetis.h" ; then
		    with_parmetis="no"
		    AC_MSG_RESULT(failed)
		else
		    AC_MSG_RESULT(yes)
		fi
	    else
		AC_MSG_RESULT(yes)
	    fi
	fi
	])
  AC_ARG_WITH([parmetis-lib],
    [AC_HELP_STRING([--with-parmetis-lib],  [name of the parmetis libraries (default is parmetis)])],
    ,[with_parmetis_lib=parmetis])
  AC_ARG_WITH([metis-lib],
    [AC_HELP_STRING([--with-metis-lib],  [name of the metis libraries (default is metis)])],
    ,[with_metis_lib=metis])

  #echo with_parmetis=$with_parmetis
  # store old values
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"
     
  # call IMDX_LIB_METIS directly and not via AC_REQUIRE
  # because AC_REQUIRE support not allow parameters
  # without any parameters a missing METIS would halt configure
dnl joseba  IMMDX_LIB_METIS(,[true])
  
  ## do nothing if --without-parmetis is used
  if test x"$with_mpi" != x"no" && test x"$with_parmetis" != x"no" ; then

      # defaultpath
      PARMETIS_LIB_PATH="$with_parmetis/$lib_path"
      PARMETIS_INCLUDE_PATH="$with_parmetis/$include_path"
                  
      PARMETIS_LIBS="-L$PARMETIS_LIB_PATH -l$with_metis_lib $DUNEMPILIBS -lm"
      PARMETIS_LDFLAGS="$DUNEMPILDFLAGS"

      # set variables so that tests can use them
      CPPFLAGS="$CPPFLAGS -I$PARMETIS_INCLUDE_PATH $METIS_INCLUDE $DUNEMPICPPFLAGS"

      # check for central header
      AC_CHECK_HEADER([parmetis.h],[
	      PARMETIS_CPPFLAGS="-I$PARMETIS_INCLUDE_PATH $METIS_INCLUDE"
	      HAVE_PARMETIS="1"],[
	      HAVE_PARMETIS="0"
	      AC_MSG_WARN([parmetis.h not found in $PARMETIS_INCLUDE_PATH with $CPPFLAGS])]
      )

      PARMETIS_CPPFLAGS="${DUNEMPICPPFLAGS} ${PARMETIS_CPPFLAGS} -DENABLE_PARMETIS=1"

#      AC_LANG_PUSH([C++])
       #Joseba change
      acx_parmetis_save_FCFLAGS="$FCFLAGS"
      FCFLAGS="$PARMETIS_CPPFLAGS $acx_parmetis_save_FCFLAGS"
      acx_parmetis_save_CFLAGS="$CFLAGS"
      CFLAGS="$PARMETIS_CPPFLAGS $acx_parmetis_save_CFLAGS"
      
      # if header is found check for the libs

      LIBS="$DUNEMPILIBS -lm $LIBS"
      
      if test x$HAVE_PARMETIS = x1 ; then
	  DUNE_CHECK_LIB_EXT([$PARMETIS_LIB_PATH], [$with_metis_lib], [metis_partgraphkway],
              [
		  PARMETIS_LIBS="-L$PARMETIS_LIB_PATH -l$with_metis_lib $METIS_LDFLAGS $DUNEMPILIBS -lm"
		  LIBS="$PARMETIS_LIBS $ac_save_LIBS"
              ],[
		  HAVE_PARMETIS="0"
		  AC_MSG_WARN(libmetis not found!)
              ])
      fi

      if test x$HAVE_PARMETIS = x1 ; then
	  DUNE_CHECK_LIB_EXT([$PARMETIS_LIB_PATH], [$with_parmetis_lib], [parmetis_v3_partkway],
              [
		  PARMETIS_LIBS="-L$PARMETIS_LIB_PATH -l$with_parmetis_lib $PARMETIS_LIBS"
              ],[
		  HAVE_PARMETIS="0"
		  AC_MSG_WARN(libparmetis not found!)
              ])
      fi

#      AC_LANG_POP([C++])
      
      # pre-set variable for summary
      #with_parmetis="no"
      
      # did it work?
      AC_MSG_CHECKING(ParMETIS in $with_parmetis)
      if test x$HAVE_PARMETIS = x1 ; then
	  AC_SUBST(PARMETIS_LIBS, $PARMETIS_LIBS)
	  AC_SUBST(PARMETIS_LDFLAGS, $PARMETIS_LDFLAGS)
	  AC_SUBST(PARMETIS_CPPFLAGS, $PARMETIS_CPPFLAGS)
	  AC_DEFINE(HAVE_PARMETIS,ENABLE_PARMETIS,[Define if you have the Parmetis library.
		  This is only true if MPI was found by configure 
		  _and_ if the application uses the PARMETIS_CPPFLAGS])
	  AC_MSG_RESULT(ok)
	  
    # add to global list
	  DUNE_PKG_LIBS="$PARMETIS_LIBS $DUNE_PKG_LIBS"
	  DUNE_PKG_LDFLAGS="$DUNE_PKG_LDFLAGS $PARMETIS_LDFLAGS"
	  DUNE_PKG_CPPFLAGS="$DUNE_PKG_CPPFLAGS $PARMETIS_CPPFLAGS"
	  
          # Joseba change
          LIBS="$ac_save_LIBS $PARMETIS_LIBS"
          LDFLAGS="$ac_save_LDFLAGS $PARMETIS_LIBS"
          CPPFLAGS="$ac_save_CPPFLAGS $PARMETIS_CPPFLAGS"
          # AC_DEFINE(HAVE_METIS, 1, [This is defined when we should compile with METIS support (default).])
    # re-set variable correctly
	  with_parmetis="yes"
      else
	  with_parmetis="no"
	  AC_MSG_RESULT(failed)
      fi 
      
  # end of "no --without-parmetis"
  else
  	with_parmetis="no"
  fi
 
  # tell automake	
  AM_CONDITIONAL(PARMETIS, test x$HAVE_PARMETIS = x1)
  
  # restore variables. Joseba commented
  dnl LDFLAGS="$ac_save_LDFLAGS"
  dnl CPPFLAGS="$ac_save_CPPFLAGS"
  dnl LIBS="$ac_save_LIBS" # Joseba change
  dnl echo "PARMETIS found"
  dnl LIBS="$PARMETIS_LIBS $ac_save_LIBS"
  dnl LDFLAGS="$PARMETIS_LIBS $ac_save_LDFLAGS"
  dnl CPPFLAGS="$PARMETIS_CPPFLAGS $ac_save_CPPFLAGS"
  
  DUNE_ADD_SUMMARY_ENTRY([ParMETIS],[$with_parmetis])

])
