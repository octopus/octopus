AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([Poisson solver],[1.8],[Damien.Caliste@cea.fr],[psolver])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([config/m4])

dnl Init something, don't know exactly what...
AM_INIT_AUTOMAKE([-Wno-portability])
dnl To keep the object files in the same structure as the source files (pb does not work with libABINIT)
dnl AM_INIT_AUTOMAKE([--warnings=no-portability, subdir-objects])

dnl Versioning
PSOLVER_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
PSOLVER_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
PSOLVER_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
if test -z "$PSOLVER_MICRO_VERSION" ; then
  PSOLVER_MICRO_VERSION=0
fi
PSOLVER_NANO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f4)
AC_SUBST(PSOLVER_MAJOR_VERSION)
AC_SUBST(PSOLVER_MINOR_VERSION)
AC_SUBST(PSOLVER_MICRO_VERSION)
AC_SUBST(PSOLVER_NANO_VERSION)

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl We defined CC for the bindings and CXX for FFTW.
AC_PROG_CC()
AC_PROG_CXX()
AM_PROG_CC_C_O()

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])

dnl Fortran 77 (to be removed)
AC_PROG_F77([], [Fortran 77])

dnl Defining Fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

AC_FC_WRAPPERS()
AX_FC_DEFS()
AX_FC_OPENMP()

dnl Set if the module name are capitalized or not with the FC compiler. 
AX_FC_MOD()
AM_CONDITIONAL(CAPITALIZE, test "$ax_fc_mod_capitalize" = "yes")
AC_SUBST(MODULE_EXT, $ax_fc_mod_ext)
dnl default installation directories
AX_FC_MODULEDIR()

dnl Set the FCFLAGS and FFLAGS variable
if test -z "$FFLAGS" ; then
  AC_SUBST(FFLAGS, "-O2")
fi
if test -n "$FCFLAGS_f90" ; then
  FCFLAGS="$FCFLAGS $FCFLAGS_f90"
fi
if test -z "$FCFLAGS" ; then
  AC_SUBST(FCFLAGS, "-O2")
fi

dnl Look for ranlib and canonicalize it with the $host variable
AC_CHECK_TOOL(RANLIB, ranlib, :)
if test "$ranlib" = ":" ; then
  AC_MSG_ERROR(["No 'ranlib' program found."])
fi
dnl Look for ar and canonicalize it with the $host variable
AC_CHECK_TOOL(AR, ar, :)
if test "$ar" = ":" ; then
  AC_MSG_ERROR(["No 'ar' program found."])
fi

dnl Enable or not MPI.
AX_MPI()

dnl Python YAML support.
AX_PYYAML()

AX_FLIB()

dnl Produce dynamic excecutables.
AX_DYNAMIC_LIBRARIES([], [f_lib_initialize], [$LIB_FUTILE_LIBS])

dnl Test accelerators
AX_ACC_CUDA()

AX_LINALG()
AX_ETSF_IO()

dnl Agregate all the optional libs.
LIB_PSOLVER_DEPS="$LIBETSFIO_LIBS $LIB_FUTILE_LIBS"
dnl if the futile libraries are dynamic add the dependencies
if test $ax_FUTILE_static != "yes" ; then
  LIB_PSOLVER_DEPS="$LIB_PSOLVER_DEPS $LIBCUDA_LIBS $LINALG_LIBS $MPI_LDFLAGS $MPI_LIBS"
fi
LIB_PSOLVER_PKG=""
if test $ax_build_dynamic != "yes" ; then
  LIB_PSOLVER_PKG="$LIB_PSOLVER_DEPS"
fi
AC_SUBST(LIB_PSOLVER_PKG)
AC_SUBST(LIB_PSOLVER_DEPS)
dnl Test for the developers' documentation.
AC_CHECK_PROG(ac_have_doxygen, [doxygen], [yes], [no])
AC_ARG_ENABLE(devel-doc,
              AS_HELP_STRING([--enable-devel-doc],
                             [compile documentation for developers using Doxygen.]),
              ac_devel_doc=$enableval, ac_devel_doc="no")
AM_CONDITIONAL(BUILD_DEVDOC, test x"$ac_have_doxygen" = x"yes" -a x"$ac_devel_doc" = x"yes")
dnl Test for the documentation.
AC_CHECK_PROG(ac_have_pdflatex, [pdflatex], [yes], [no])
AM_CONDITIONAL(BUILD_PDFDOC, test x"$ac_have_pdflatex" = x"yes" -a x"$ac_devel_doc" = x"yes")

dnl LIB_PSOLVER_DEPS="$LIBETSFIO_LIBS $LIBCUDA_LIBS $LINALG_LIBS $MPI_LDFLAGS $MPI_LIBS"
dnl LIB_PSOLVER_PKG=""
dnl if test $ax_build_dynamic != "yes" ; then
dnl   LIB_PSOLVER_PKG="$LIB_PSOLVER_DEPS"
dnl fi
dnl AC_SUBST(LIB_PSOLVER_PKG)
dnl LIB_PSOLVER_DEPS="$LIB_PSOLVER_DEPS $LIB_FUTILE_LIBS"
dnl AC_SUBST(LIB_PSOLVER_DEPS)


dnl default installation directories
docdir="${docdir}"
AC_SUBST(docdir)

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
psolver.pc
config/Doxyfile
doc/Makefile
src/configure.inc
src/Makefile
tests/Makefile
])

AC_OUTPUT

echo "
Poisson-solver configure summary:
================================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 flags:           $FCFLAGS
  Fortran77 compiler:        $F77
  Fortran77 flags:           $FFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIB_PSOLVER_DEPS $LIBS

Build:
  Dynamic library:           $ax_build_dynamic

Options:
  Max. Memory per task (GB): $MEMORY_LIMIT
  With OpenMP:               $ax_fc_openmp_msg
  With MPI:                  $ac_use_mpi
   | include dir.:           $MPI_INCLUDE
   | linker flags:           $MPI_LDFLAGS
   | linked libraries:       $MPI_LIBS
   | MPI_INIT_THREAD avail:  $ac_with_mpinithread
   | MPI2 support:           $ac_with_mpi2
  With Cuda support:         $ax_have_acc_cuda
   | NVidia Cuda Compiler:   $NVCC
   | Cuda flags:             $CUDA_FLAGS
  With ETSF_IO:              $ax_have_etsf_io
   | include dir.:           $LIBETSFIO_INCLUDE

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Static libraries:          $libdir
  Fortran modules:           $moduledir
  Documentation:             $docdir
"

if test x"$ac_use_mpi" = x"warn" -o x"$ax_have_etsf_io" = x"warn" -o x"$withscalapack" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
  if test x"$ax_have_etsf_io" = x"warn" ; then
    echo "  ETSF_IO:                   disabled by tests"
  fi
  if test x"$withscalapack" = x"warn" ; then
    echo "  Scalapack binding:         disabled by tests"
  fi

fi

dnl Do not touch this line. On the Cray machine in Basel, the preprocessor
dnl in Fortran sources does not allow C comments (/*  */), so we remove them.
$SED -n '1h;1!H;${;g;s%/\*[[^*]]*\*/%%g;p;}' config.h > config.inc
