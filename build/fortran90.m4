## Written by Miguel Marques (based on fortran.m4)

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

## ----------------------- ##
## 1. Language selection.  ##
## ----------------------- ##


# ----------------------------- #
# 1d. The Fortran 90 language.  #
# ----------------------------- #


# AC_LANG(Fortran 90)
# -------------------
m4_define([AC_LANG(Fortran 90)],
[ac_ext=f90
ac_compile='$F90 -c $F90FLAGS conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_link='$F90 -o conftest$ac_exeext $F90FLAGS $LDFLAGS conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'
ac_compiler_gnu=$ac_cv_f90_compiler_gnu
])


# AC_LANG_FORTRAN90
# -----------------
AU_DEFUN([AC_LANG_FORTRAN90], [AC_LANG(Fortran 90)])


# _AC_LANG_ABBREV(Fortran 90)
# ---------------------------
m4_define([_AC_LANG_ABBREV(Fortran 90)], [f90])



## ---------------------- ##
## 2.Producing programs.  ##
## ---------------------- ##


# ------------------------ #
# 2d. Fortran 90 sources.  #
# ------------------------ #

# AC_LANG_SOURCE(Fortran 90)(BODY)
# --------------------------------
# FIXME: Apparently, according to former AC_TRY_COMPILER, the CPP
# directives must not be included.  But AC_TRY_RUN_NATIVE was not
# avoiding them, so?
m4_define([AC_LANG_SOURCE(Fortran 90)],
[$1])


# AC_LANG_PROGRAM(Fortran 90)([PROLOGUE], [BODY])
# -----------------------------------------------
# Yes, we discard the PROLOGUE.
m4_define([AC_LANG_PROGRAM(Fortran 90)],
[m4_ifval([$1],
       [m4_warn([syntax], [$0: ignoring PROLOGUE: $1])])dnl
      program main
$2
      end program main])


# AC_LANG_CALL(Fortran 90)(PROLOGUE, FUNCTION)
# --------------------------------------------
# FIXME: This is a guess, help!
m4_define([AC_LANG_CALL(Fortran 90)],
[AC_LANG_PROGRAM([$1],
[      call $2])])




## -------------------------------------------- ##
## 3. Looking for Compilers and Preprocessors.  ##
## -------------------------------------------- ##


# ----------------------------- #
# 3d. The Fortran 90 compiler.  #
# ----------------------------- #


# AC_LANG_PREPROC(Fortran 90)
# ---------------------------
# Find the Fortran 90 preprocessor.  Must be AC_DEFUN'd to be AC_REQUIRE'able.
AC_DEFUN([AC_LANG_PREPROC(Fortran 90)],
[
  # this should not be hardwired
  if test -z "$F90CPP"; then F90CPP="/lib/cpp -C"; fi
	AC_SUBST(F90CPP)
])


# AC_LANG_COMPILER(Fortran 90)
# ----------------------------
# Find the Fortran 90 compiler.  Must be AC_DEFUN'd to be
# AC_REQUIRE'able.
AC_DEFUN([AC_LANG_COMPILER(Fortran 90)],
[AC_REQUIRE([AC_PROG_F90])])


# ac_cv_prog_g90
# --------------
# We used to name the cache variable this way.
AU_DEFUN([ac_cv_prog_g90],
[ac_cv_f90_compiler_gnu])


# AC_PROG_F90([COMPILERS...])
# ---------------------------
# COMPILERS is a space separated list of Fortran 90 compilers to search
# for.
#
# Compilers are ordered by
#  1. F90, F95
#  2. Good/tested native compilers, bad/untested native compilers
AC_DEFUN([AC_PROG_F90],
[AC_LANG_PUSH(Fortran 90)dnl
AC_ARG_VAR([F90],    [Fortran 90 compiler command])dnl
AC_ARG_VAR([F90FLAGS], [Fortran 90 compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl

if test x"$2" = x"MPI" ; then
	AC_CHECK_TOOLS(F90, [m4_default([$1],
                  [mpif90 mpxlf90])])
else
	AC_CHECK_TOOLS(F90, [m4_default([$1],
                  [f90 abf90 ifc xlf90 pgf90 epcf90 f95 fort xlf95 lf95 g95])])
fi

# Provide some information about the compiler.
echo "$as_me:__oline__:" \
     "checking for _AC_LANG compiler version" >&AS_MESSAGE_LOG_FD
ac_compiler=`set X $ac_compile; echo $[2]`
_AC_EVAL([$ac_compiler --version </dev/null >&AS_MESSAGE_LOG_FD])
_AC_EVAL([$ac_compiler -v </dev/null >&AS_MESSAGE_LOG_FD])
_AC_EVAL([$ac_compiler -V </dev/null >&AS_MESSAGE_LOG_FD])

m4_expand_once([_AC_COMPILER_EXEEXT])[]dnl
m4_expand_once([_AC_COMPILER_OBJEXT])[]dnl
# If we don't use `.F' as extension, the preprocessor is not run on the
# input file.
ac_save_ext=$ac_ext
ac_ext=F90
_AC_LANG_COMPILER_GNU
ac_ext=$ac_save_ext
G90=`test $ac_compiler_gnu = yes && echo yes`
dnl _AC_PROG_F90_G

dnl let us see if we know the fortran compiler
AC_REQUIRE([AC_CANONICAL_HOST])

if test -z "${F90FLAGS}"; then
	case "${host}" in
	i?86*linux*)
		case "${F90}" in
		pgf90*)
  		F90FLAGS="$F90FLAGS -O2 -fast -Munroll -Mnoframe -Mdalign"
		;;
		abf90*)
			F90FLAGS="$F90FLAGS -O -YEXT_NAMES=LCS -YEXT_SFX=_"
		;;
		ifc*)
			F90FLAGS="$F90FLAGS -nbs -lowercase -align -zero -r16 -fpp1 -pad -unroll -tpp6 -O3 -rcd -xW -mp1"
		;;
		*)
			F90FLAGS="$F90FLAGS -O"
		esac
	;;
	alphaev*)
		F90FLAGS="$F90FLAGS -align dcommons -fast -tune host -arch host -noautomatic"
		;;
	powerpc-ibm*)
		F90FLAGS="$F90FLAGS -bmaxdata:0x80000000 -qmaxmem=-1 -qsuffix=f=f90 -Q -O5 -qstrict -qtune=auto -qarch=auto -qhot -qipa"
	  ;;
	mips-sgi-irix*)
		F90FLAGS="$F90FLAGS -O3 -INLINE -n32 -LANG:recursive=on"
		;;
	*)
		F90FLAGS="$F90FLAGS -O"
	esac
fi

AC_LANG_POP(Fortran 90)dnl
])# AC_PROG_F90


# _AC_PROG_F90_G
# --------------
# Check whether -g works, even if F90FLAGS is set, in case the package
# plays around with F90FLAGS (such as to build both debugging and normal
# versions of a library), tasteless as that idea is.
m4_define([_AC_PROG_F90_G],
[ac_test_F90FLAGS=${F90FLAGS+set}
ac_save_F90FLAGS=$F90FLAGS
F90FLAGS=
AC_CACHE_CHECK(whether $F90 accepts -g, ac_cv_prog_f90_g,
[F90FLAGS=-g
_AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_f90_g=yes],
[ac_cv_prog_f90_g=no])
])
if test "$ac_test_F90FLAGS" = set; then
  F90FLAGS=$ac_save_F90FLAGS
elif test $ac_cv_prog_f90_g = yes; then
  F90FLAGS="-g"
else
  F90FLAGS=
fi[]dnl
])# _AC_PROG_F90_G


# AC_PROG_F90_C_O
# ---------------
# Test if the Fortran 90 compiler accepts the options `-c' and `-o'
# simultaneously, and define `F90_NO_MINUS_C_MINUS_O' if it does not.
#
# The usefulness of this macro is questionable, as I can't really see
# why anyone would use it.  The only reason I include it is for
# completeness, since a similar test exists for the C compiler.
AC_DEFUN([AC_PROG_F90_C_O],
[AC_REQUIRE([AC_PROG_F90])dnl
AC_CACHE_CHECK([whether $F90 understand -c and -o together],
               [ac_cv_prog_f90_c_o],
[AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
# We test twice because some compilers refuse to overwrite an existing
# `.o' file with `-o', although they will create one.
ac_try='$F90 $F90FLAGS -c conftest.$ac_ext -o conftest.$ac_objext >&AS_MESSAGE_LOG_FD'
if AC_TRY_EVAL(ac_try) &&
     test -f conftest.$ac_objext &&
     AC_TRY_EVAL(ac_try); then
  ac_cv_prog_f90_c_o=yes
else
  ac_cv_prog_f90_c_o=no
fi
rm -f conftest*])
if test $ac_cv_prog_f90_c_o = no; then
  AC_DEFINE(F90_NO_MINUS_C_MINUS_O, 1,
            [Define to 1 if your Fortran 90 compiler doesn't accept
             -c and -o together.])
fi
])# AC_PROG_F90_C_O





## ------------------------------- ##
## 4. Compilers' characteristics.  ##
## ------------------------------- ##


# ---------------------------------------- #
# 4d. Fortran 90 compiler characteristics. #
# ---------------------------------------- #


# _AC_PROG_F90_V_OUTPUT([FLAG = $ac_cv_prog_f90_v])
# -------------------------------------------------
# Link a trivial Fortran program, compiling with a verbose output FLAG
# (which default value, $ac_cv_prog_f90_v, is computed by
# _AC_PROG_F90_V), and return the output in $ac_f90_v_output.  This
# output is processed in the way expected by AC_F90_LIBRARY_LDFLAGS,
# so that any link flags that are echoed by the compiler appear as
# space-separated items.
AC_DEFUN([_AC_PROG_F90_V_OUTPUT],
[AC_REQUIRE([AC_PROG_F90])dnl
AC_LANG_PUSH(Fortran 90)dnl

AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])

# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran 90 compiler in order to get
# "verbose" output that we can then parse for the Fortran 90 linker
# flags.
ac_save_F90FLAGS=$F90FLAGS
F90FLAGS="$F90FLAGS m4_default([$1], [$ac_cv_prog_f90_v])"
(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
ac_f90_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_f90_v_output" >&AS_MESSAGE_LOG_FD
F90FLAGS=$ac_save_F90FLAGS

rm -f conftest*
AC_LANG_POP(Fortran 90)dnl

# If we are using xlf then replace all the commas with spaces.
if echo $ac_f90_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_f90_v_output=`echo $ac_f90_v_output | sed 's/,/ /g'`
fi

# On HP/UX there is a line like: "LPATH is: /foo:/bar:/baz" where
# /foo, /bar, and /baz are search directories for the Fortran linker.
# Here, we change these into -L/foo -L/bar -L/baz (and put it first):
ac_f90_v_output="`echo $ac_f90_v_output | 
	grep 'LPATH is:' |
	sed 's,.*LPATH is\(: *[[^ ]]*\).*,\1,;s,: */, -L/,g'` $ac_f90_v_output"

# If we are using Cray Fortran then delete quotes.
# Use "\"" instead of '"' for font-lock-mode.
# FIXME: a more general fix for quoted arguments with spaces?
if echo $ac_f90_v_output | grep cft90 >/dev/null 2>&1; then
  ac_f90_v_output=`echo $ac_f90_v_output | sed "s/\"//g"`
fi[]dnl
])# _AC_PROG_F90_V_OUTPUT


# _AC_PROG_F90_V
# --------------
#
# Determine the flag that causes the Fortran 90 compiler to print
# information of library and object files (normally -v)
# Needed for AC_F90_LIBRARY_FLAGS
# Some compilers don't accept -v (Lahey: -verbose, xlf: -V, Fujitsu: -###)
AC_DEFUN([_AC_PROG_F90_V],
[AC_CACHE_CHECK([how to get verbose linking output from $F90],
                [ac_cv_prog_f90_v],
[AC_LANG_ASSERT(Fortran 90)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_f90_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose --verbose -V -\#\#\#; do
  _AC_PROG_F90_V_OUTPUT($ac_verb)
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_f90_v_output; do
     case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
          ac_cv_prog_f90_v=$ac_verb
          break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_f90_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $F90])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AC_PROG_F90_V


# AC_F90_LIBRARY_LDFLAGS
# ----------------------
#
# Determine the linker flags (e.g. "-L" and "-l") for the Fortran 90
# intrinsic and run-time libraries that are required to successfully
# link a Fortran 90 program or shared library.  The output variable
# F90LIBS is set to these flags.
#
# This macro is intended to be used in those situations when it is
# necessary to mix, e.g. C++ and Fortran 90, source code into a single
# program or shared library.
#
# For example, if object files from a C++ and Fortran 90 compiler must
# be linked together, then the C++ compiler/linker must be used for
# linking (since special C++-ish things need to happen at link time
# like calling global constructors, instantiating templates, enabling
# exception support, etc.).
#
# However, the Fortran 90 intrinsic and run-time libraries must be
# linked in as well, but the C++ compiler/linker doesn't know how to
# add these Fortran 90 libraries.  Hence, the macro
# "AC_F90_LIBRARY_LDFLAGS" was created to determine these Fortran 90
# libraries.
#
# This macro was packaged in its current form by Matthew D. Langston.
# However, nearly all of this macro came from the "OCTAVE_FLIBS" macro
# in "octave-2.0.13/aclocal.m4", and full credit should go to John
# W. Eaton for writing this extremely useful macro.  Thank you John.
AC_DEFUN([AC_F90_LIBRARY_LDFLAGS],
[AC_LANG_PUSH(Fortran 90)dnl
_AC_PROG_F90_V
AC_CACHE_CHECK([for Fortran 90 libraries], ac_cv_flibs,
[if test "x$F90LIBS" != "x"; then
  ac_cv_flibs="$F90LIBS" # Let the user override the test.
else

_AC_PROG_F90_V_OUTPUT

ac_cv_flibs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_f90_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f90libs, ,
              ac_cv_f90libs="$ac_cv_f90libs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f90libs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_f90libs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt0.o | -lc | -lgcc | -libmil | -LANG:=*)
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_f90libs="$ac_cv_f90libs $ac_arg"
          ;;
        -[[LRuY]])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
	  case $[2] in
             "" | -*);;
             *)
		ac_arg="$ac_arg$[2]"
		shift; shift
		set X $ac_arg "$[@]"
		;;
	  esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_f90libs, ,
                               [ac_arg="$ac_arg $ac_j"
                               ac_cv_f90libs="$ac_cv_f90libs $ac_j"])
          done
          ;;
				*.link) # necessary for the ifc compiler
					ac_cv_f90libs="$ac_cv_f90libs $ac_arg"
					;;
				-lowercase) # necessary for the ifc compiler
					;;
				-link) # necessary for xlf90
					;;
        -[[lLR]]*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f90libs, ,
                             ac_cv_f90libs="$ac_cv_f90libs $ac_arg")
          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_f90_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_f90libs)
      ;;
esac
fi # test "x$F90LIBS" = "x"
])
F90LIBS="$ac_cv_f90libs"
AC_SUBST(F90LIBS)
AC_LANG_POP(Fortran 90)dnl
])# AC_F90_LIBRARY_LDFLAGS


# AC_F90_DUMMY_MAIN([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------
#
# Detect name of dummy main routine required by the Fortran libraries,
# (if any) and define F90_DUMMY_MAIN to this name (which should be
# used for a dummy declaration, if it is defined).  On some systems,
# linking a C program to the Fortran library does not work unless you
# supply a dummy function called something like MAIN__.
#
# Execute ACTION-IF-NOT-FOUND if no way of successfully linking a C
# program with the F90 libs is found; default to exiting with an error
# message.  Execute ACTION-IF-FOUND if a dummy routine name is needed
# and found or if it is not needed (default to defining F90_DUMMY_MAIN
# when needed).
#
# What is technically happening is that the Fortran libraries provide
# their own main() function, which usually initializes Fortran I/O and
# similar stuff, and then calls MAIN__, which is the entry point of
# your program.  Usually, a C program will override this with its own
# main() routine, but the linker sometimes complain if you don't
# provide a dummy (never-called) MAIN__ routine anyway.
#
# Of course, programs that want to allow Fortran subroutines to do
# I/O, etcetera, should call their main routine MAIN__() (or whatever)
# instead of main().  A separate autoconf test (AC_F90_MAIN) checks
# for the routine to use in this case (since the semantics of the test
# are slightly different).  To link to e.g. purely numerical
# libraries, this is normally not necessary, however, and most C/C++
# programs are reluctant to turn over so much control to Fortran.  =)
#
# The name variants we check for are (in order):
#   MAIN__ (g77, MAIN__ required on some systems; IRIX, MAIN__ optional)
#   MAIN_, __main (SunOS)
#   MAIN _MAIN __MAIN main_ main__ _main (we follow DDD and try these too)
AC_DEFUN([AC_F90_DUMMY_MAIN],
[AC_REQUIRE([AC_F90_LIBRARY_LDFLAGS])dnl
m4_define([_AC_LANG_PROGRAM_C_F90_HOOKS],
[#ifdef F90_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F90_DUMMY_MAIN() { return 1; }
#endif
])
AC_CACHE_CHECK([for dummy main to link with Fortran 90 libraries],
               ac_cv_f90_dummy_main,
[AC_LANG_PUSH(C)dnl
 ac_f90_dm_save_LIBS=$LIBS
 LIBS="$LIBS $F90LIBS"

 # First, try linking without a dummy main:
 AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])],
                [ac_cv_f90_dummy_main=none],
                [ac_cv_f90_dummy_main=unknown])

 if test $ac_cv_f90_dummy_main = unknown; then
   for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[@%:@define F90_DUMMY_MAIN $ac_func]])],
                    [ac_cv_f90_dummy_main=$ac_func; break])
   done
 fi
 rm -f conftest*
 LIBS=$ac_f90_dm_save_LIBS
 AC_LANG_POP(C)dnl
])
F90_DUMMY_MAIN=$ac_cv_f90_dummy_main
AS_IF([test "$F90_DUMMY_MAIN" != unknown],
      [m4_default([$1],
[if test $F90_DUMMY_MAIN != none; then
  AC_DEFINE_UNQUOTED([F90_DUMMY_MAIN], $F90_DUMMY_MAIN,
                     [Define to dummy `main' function (if any) required to
                      link to the Fortran 90 libraries.])
fi])],
      [m4_default([$2],
                [AC_MSG_ERROR([linking to Fortran libraries from C fails])])])
])# AC_F90_DUMMY_MAIN


# AC_F90_MAIN
# -----------
# Define F90_MAIN to name of alternate main() function for use with
# the Fortran libraries.  (Typically, the libraries may define their
# own main() to initialize I/O, etcetera, that then call your own
# routine called MAIN__ or whatever.)  See AC_F90_DUMMY_MAIN, above.
# If no such alternate name is found, just define F90_MAIN to main.
#
AC_DEFUN([AC_F90_MAIN],
[AC_REQUIRE([AC_F90_LIBRARY_LDFLAGS])dnl
AC_CACHE_CHECK([for alternate main to link with Fortran 90 libraries],
               ac_cv_f90_main,
[AC_LANG_PUSH(C)dnl
 ac_f90_m_save_LIBS=$LIBS
 LIBS="$LIBS $F90LIBS"
 ac_cv_f90_main="main" # default entry point name

 for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
   AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@undef F90_DUMMY_MAIN
@%:@define main $ac_func])],
                  [ac_cv_f90_main=$ac_func; break])
 done
 rm -f conftest*
 LIBS=$ac_f90_m_save_LIBS
 AC_LANG_POP(C)dnl
])
AC_DEFINE_UNQUOTED([F90_MAIN], $ac_cv_f90_main,
                   [Define to alternate name for `main' routine that is
                    called from a `main' in the Fortran libraries.])
])# AC_F90_MAIN


# _AC_F90_NAME_MANGLING
# ---------------------
# Test for the name mangling scheme used by the Fortran 90 compiler.
#
# Sets ac_cv_f90_mangling. The value contains three fields, separated
# by commas:
#
# lower case / upper case:
#    case translation of the Fortran 90 symbols
# underscore / no underscore:
#    whether the compiler appends "_" to symbol names
# extra underscore / no extra underscore:
#    whether the compiler appends an extra "_" to symbol names already
#    containing at least one underscore
#
AC_DEFUN([_AC_F90_NAME_MANGLING],
[AC_REQUIRE([AC_F90_LIBRARY_LDFLAGS])dnl
AC_REQUIRE([AC_F90_DUMMY_MAIN])dnl
AC_CACHE_CHECK([for Fortran 90 name-mangling scheme],
               ac_cv_f90_mangling,
[AC_LANG_PUSH(Fortran 90)dnl
AC_COMPILE_IFELSE(
[      subroutine foobar()
      return
      end subroutine foobar
      subroutine foo_bar()
      return
      end subroutine foo_bar],
[mv conftest.$ac_objext cf90_test.$ac_objext

  AC_LANG_PUSH(C)dnl

  ac_save_LIBS=$LIBS
  LIBS="cf90_test.$ac_objext $LIBS $F90LIBS"

  ac_success=no
  for ac_foobar in foobar FOOBAR; do
    for ac_underscore in "" "_"; do
      ac_func="$ac_foobar$ac_underscore"
      AC_TRY_LINK_FUNC($ac_func,
         [ac_success=yes; break 2])
    done
  done

  if test "$ac_success" = "yes"; then
     case $ac_foobar in
        foobar)
           ac_case=lower
           ac_foo_bar=foo_bar
           ;;
        FOOBAR)
           ac_case=upper
           ac_foo_bar=FOO_BAR
           ;;
     esac

     ac_success_extra=no
     for ac_extra in "" "_"; do
        ac_func="$ac_foo_bar$ac_underscore$ac_extra"
        AC_TRY_LINK_FUNC($ac_func,
        [ac_success_extra=yes; break])
     done

     if test "$ac_success_extra" = "yes"; then
	ac_cv_f90_mangling="$ac_case case"
        if test -z "$ac_underscore"; then
           ac_cv_f90_mangling="$ac_cv_f90_mangling, no underscore"
	else
           ac_cv_f90_mangling="$ac_cv_f90_mangling, underscore"
        fi
        if test -z "$ac_extra"; then
           ac_cv_f90_mangling="$ac_cv_f90_mangling, no extra underscore"
	else
           ac_cv_f90_mangling="$ac_cv_f90_mangling, extra underscore"
        fi
      else
	ac_cv_f90_mangling="unknown"
      fi
  else
     ac_cv_f90_mangling="unknown"
  fi

  LIBS=$ac_save_LIBS
  AC_LANG_POP(C)dnl
  rm -f cf90_test* conftest*])
AC_LANG_POP(Fortran 90)dnl
])
])# _AC_F90_NAME_MANGLING

# The replacement is empty.
AU_DEFUN([AC_F90_NAME_MANGLING], [])


# AC_F90_WRAPPERS
# ---------------
# Defines C macros F90_FUNC(name,NAME) and F90_FUNC_(name,NAME) to
# properly mangle the names of C identifiers, and C identifiers with
# underscores, respectively, so that they match the name mangling
# scheme used by the Fortran 90 compiler.
AC_DEFUN([AC_F90_WRAPPERS],
[AC_REQUIRE([_AC_F90_NAME_MANGLING])dnl
AH_TEMPLATE([F90_FUNC],
    [Define to a macro mangling the given C identifier (in lower and upper
     case), which must not contain underscores, for linking with Fortran.])dnl
AH_TEMPLATE([F90_FUNC_],
    [As F90_FUNC, but for C identifiers containing underscores.])dnl
case $ac_cv_f90_mangling in
  "lower case, no underscore, no extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [name])
          AC_DEFINE([F90_FUNC_(name,NAME)], [name]) ;;
  "lower case, no underscore, extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [name])
          AC_DEFINE([F90_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, no extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE([F90_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE([F90_FUNC_(name,NAME)], [name ## __]) ;;
  "upper case, no underscore, no extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [NAME])
          AC_DEFINE([F90_FUNC_(name,NAME)], [NAME]) ;;
  "upper case, no underscore, extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [NAME])
          AC_DEFINE([F90_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, no extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE([F90_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE([F90_FUNC_(name,NAME)], [NAME ## __]) ;;
  *)
          AC_MSG_WARN([unknown Fortran 90 name-mangling scheme])
          ;;
esac
])# AC_F90_WRAPPERS


# AC_F90_FUNC(NAME, [SHELLVAR = NAME])
# ------------------------------------
# For a Fortran subroutine of given NAME, define a shell variable
# $SHELLVAR to the Fortran-90 mangled name.  If the SHELLVAR
# argument is not supplied, it defaults to NAME.
AC_DEFUN([AC_F90_FUNC],
[AC_REQUIRE([_AC_F90_NAME_MANGLING])dnl
case $ac_cv_f90_mangling in
  upper*) ac_val="m4_toupper([$1])" ;;
  lower*) ac_val="m4_tolower([$1])" ;;
  *)      ac_val="unknown" ;;
esac
case $ac_cv_f90_mangling in *," underscore"*) ac_val="$ac_val"_ ;; esac
m4_if(m4_index([$1],[_]),-1,[],
[case $ac_cv_f90_mangling in *," extra underscore"*) ac_val="$ac_val"_ ;; esac
])
m4_default([$2],[$1])="$ac_val"
])# AC_F90_FUNC

################################################
# AC_LANG_FUNC_LINK_TRY(C)(FUNCTION)
# ----------------------------------
# Don't include <ctype.h> because on OSF/1 3.0 it includes
# <sys/types.h> which includes <sys/select.h> which contains a
# prototype for select.  Similarly for bzero.
m4_define([AC_LANG_FUNC_LINK_TRY(Fortran 90)],
[AC_LANG_PROGRAM([], [call [$1]])])
