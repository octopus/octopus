AC_DEFUN([ACX_LIBVDWXC],
[

  acx_libvdwxc_ok=no
  acx_libvdwxc_mpi_ok=no

  dnl BACKUP LIBS AND FCFLAGS
  acx_libvdwxc_save_LIBS="$LIBS"
  acx_libvdwxc_save_FCFLAGS="$FCFLAGS"

  dnl Check if the library was given in the command line
  AC_ARG_WITH(libvdwxc-prefix, [AS_HELP_STRING([--with-libvdwxc-prefix=DIR], [Directory where libvdwxc is installed.])])

  if test x"$FCFLAGS_LIBVDWXC" = x; then
    case $with_libvdwxc_prefix in
      "") FCFLAGS_LIBVDWXC="-I/usr/include" ;;
      *)  FCFLAGS_LIBVDWXC="-I$with_libvdwxc_prefix/include" ;;
    esac
  fi

  AC_MSG_CHECKING([for libvdwxc])

  libvdwxc_program_serial="AC_LANG_PROGRAM([],[
    implicit none

    include 'vdwxcfort.f90'

    integer, pointer :: vdw

    call vdwxc_new(1, vdw)
    call vdwxc_print(vdw)
    call vdwxc_init_serial(vdw)
    call vdwxc_finalize(vdw)
  ])"

  libvdwxc_program_mpi="AC_LANG_PROGRAM([],[
    implicit none

    include 'vdwxcfort.f90'

    integer, pointer :: vdw
    integer          :: has_mpi

    call vdwxc_new(1, vdw)
    ! We are supposed to pass a communicator like MPI_COMM_WORLD below,
    ! but I do not know how to get that here.  We pass 42 instead.
    ! The build system is happy as long as it compiles.
    call vdwxc_init_mpi(vdw, 42)
    call vdwxc_finalize(vdw)
  ])"

  FCFLAGS="$FCFLAGS_LIBVDWXC $acx_libvdwxc_save_FCFLAGS"

  if test ! -z "$with_libvdwxc_prefix"; then
    LIBS_LIBVDWXC="-L$with_libvdwxc_prefix/lib -lvdwxcfort"
  else
    LIBS_LIBVDWXC="-lvdwxcfort"
  fi

  LIBS="$LIBS_LIBVDWXC $acx_libvdwxc_save_LIBS"
  AC_LINK_IFELSE($libvdwxc_program_serial, [acx_libvdwxc_ok=yes], [acx_libvdwxc_ok=no])
  AC_LINK_IFELSE($libvdwxc_program_mpi, [acx_libvdwxc_mpi_ok=yes], [acx_libvdwxc_mpi_ok=no])


  if test x"$acx_libvdwxc_ok" = xyes; then
    AC_DEFINE(HAVE_LIBVDWXC, 1, [Define if libvdwxc is available])
    AC_MSG_RESULT([$acx_libvdwxc_ok ($FCFLAGS_LIBVDWXC $LIBS_LIBVDWXC)])
    if test x"$acx_libvdwxc_mpi_ok" = xyes; then
      AC_DEFINE(HAVE_LIBVDWXC_MPI, 1, [Define if libvdwxc has MPI support])
      AC_MSG_RESULT([$acx_libvdwxc_mpi_ok])
    else
      AC_MSG_WARN([libvdwxc not compiled with MPI, cannot be used with domain decomposition])
    fi
  else
    LIBS_LIBVDWXC=""
    FCFLAGS_LIBVDWXC=""
    AC_MSG_RESULT([$acx_libvdwxc_ok])
    AC_MSG_WARN([Could not find libvdwxc library.
                 *** Will compile without libvdwxc support])
  fi

  AC_SUBST(FCFLAGS_LIBVDWXC)
  AC_SUBST(LIBS_LIBVDWXC)

  FCFLAGS="$acx_libvdwxc_save_FCFLAGS"
  LIBS="$acx_libvdwxc_save_LIBS"
])
