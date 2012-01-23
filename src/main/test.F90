!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: $

#include "global.h"

program oct_test
  use calc_mode_m
  use command_line_m
  use datasets_m
  use derivatives_m
  use fft_m
  use global_m
  use hamiltonian_m
  use io_m
  use loct_m
  use messages_m
  use mpi_m
  use multicomm_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_calc_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  integer :: test_type
  integer :: test_mode
  integer :: ierr

  character(len=256) :: sys_name

  integer, parameter ::              &
    HARTREE_TEST       =   1,        &
    DER_TEST           =   2,        &
    ORT_TEST           =   3

  integer, parameter :: &
    TEST_REAL    = 1,   &
    TEST_COMPLEX = 2,   &
    TEST_ALL     = 3

  call getopt_init(ierr)
  if(ierr .eq. 0) call getopt_octopus()
  call getopt_end()

  call global_init()
  call calc_mode_init()
  call parser_init()
  call messages_init()

  call messages_obsolete_variable('WhichTest', 'TestMode')

  !%Variable TestMode
  !%Type integer
  !%Default hartree_test
  !%Section Utilities::oct-test
  !%Description
  !% Decides what kind of test should be performed.
  !%Option hartree_test 1
  !% Tests the various Hartree solvers.
  !%Option derivatives 2
  !% Tests the implementation of the finite-difference operators, used to calculate derivatives.
  !%Option orthogonalization 3
  !% Tests the implementation of the orthogonalization routines.
  !%End
  call parse_integer('TestMode', HARTREE_TEST, test_mode)
  call datasets_init(test_mode)

  call messages_obsolete_variable('TestDerivatives', 'TestType')
  call messages_obsolete_variable('TestOrthogonalization', 'TestType')
  
  !%Variable TestType
  !%Type integer
  !%Default all
  !%Section Utilities::oct-test
  !%Description
  !% Decides what on what type of values the test should be performed.
  !%Option real 1
  !% Tests derivatives for real functions.
  !%Option complex 2
  !% Tests derivatives for complex functions.
  !%Option all 3
  !% Tests derivatives for both real and complex functions.
  !%End
  call parse_integer('TestType', TEST_ALL, test_type)

  call io_init()
  call profiling_init()

  ! Let us print our logo
  if(mpi_grp_is_root(mpi_world)) then
    call io_dump_file(stdout, trim(trim(conf%share) // '/logo'))
  end if

  ! Let us print the version
  message(1) = ""
  message(2) = str_center("Running octopus", 70)
  message(3) = ""
  call messages_info(3)

  message(1) = &
       "Version                : " // trim(conf%version)
  message(2) = &
       "Revision               : "// trim(conf%latest_svn)
  message(3) = &
       "Build time             : "// trim(conf%build_time)
  call messages_info(3)

  write(message(1), '(a, i1)') &
       'Configuration options  : max-dim=', MAX_DIM
!!$
#ifdef HAVE_OPENMP
  message(1) = trim(message(1))//' openmp'
#endif
#ifdef HAVE_MPI
  message(1) = trim(message(1))//' mpi'
#endif
#ifdef HAVE_OPENCL
  message(1) = trim(message(1))//' opencl'
#endif
#ifdef HAVE_M128D
  message(1) = trim(message(1))//' sse2'
#endif
#ifdef HAVE_M256D
  message(1) = trim(message(1))//' avx'
#endif
#ifdef HAVE_BLUE_GENE
  message(1) = trim(message(1))//' bluegene'
#endif

  message(2) = &
       'Optional libraries     :'
#ifdef HAVE_MPI2
  message(2) = trim(message(2))//' mpi2'
#endif
#ifdef HAVE_NETCDF
  message(2) = trim(message(2))//' netcdf'
#endif
#ifdef HAVE_METIS
  message(2) = trim(message(2))//' metis'
#endif
#ifdef HAVE_GDLIB
  message(2) = trim(message(2))//' gdlib'
#endif
#ifdef HAVE_PAPI
  message(2) = trim(message(2))//' papi'
#endif
#ifdef HAVE_SPARSKIT
  message(2) = trim(message(2))//' sparskit'
#endif
#ifdef HAVE_ETSF_IO
  message(2) = trim(message(2))//' etsf_io'
#endif
#ifdef HAVE_BERKELEYGW
  message(2) = trim(message(2))//' berkeleygw'
#endif
#ifdef HAVE_PFFT
  message(2) = trim(message(2))//' pfft'
#endif
#ifdef HAVE_NFFT
  message(2) = trim(message(2))//' nfft'
#endif
#ifdef HAVE_SCALAPACK
  message(2) = trim(message(2))//' scalapack'
#endif
#ifdef HAVE_LIBFM
  message(2) = trim(message(2))//' libfm'
#endif

  message(3) = &
       'Architecture           : '// TOSTRING(OCT_ARCH)
  call messages_info(3)

  message(1) = &
       "C compiler             : "//trim(conf%cc)
  message(2) = &
       "C compiler flags       : "//trim(conf%cflags)
  message(3) = &
       "Fortran compiler       : "//trim(conf%fc)
  message(4) = &
       "Fortran compiler flags : "//trim(conf%fcflags)
  call messages_info(4)

  message(1) = ""
  call messages_info(1)

  ! Let us print where we are running
  call loct_sysname(sys_name)
  write(message(1), '(a)') str_center("The octopus is swimming in " // trim(sys_name), 70)
  message(2) = ""
  call messages_info(2)

#if defined(HAVE_MPI)
  call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

  call print_date("Calculation started on ")

  if(no_datasets > 1) then
    message(1) = 'Info: Multi-Dataset Mode'
    message(2) = 'Info: Running dataset "'//trim(current_label)//'"'
    call messages_info(2, stress = .true.)
  end if

  call messages_print_stress(stdout, "Test mode")
  call messages_print_var_option(stdout, "TestMode", test_mode)
  call messages_print_stress(stdout)

  call fft_all_init()
  call unit_system_init()

  select case(test_mode)
  case(HARTREE_TEST)
    call test_hartree()
  case(DER_TEST)
    call test_derivatives()
  case(ORT_TEST)
    call test_orthogonalization()
  end select

  call fft_all_end()
  call profiling_output()
  call profiling_end()
  call io_end()
  call print_date("Calculation ended on ")
  call datasets_end()
  call messages_end()
  call parser_end()
  call calc_mode_end()
  call global_end()

  contains

! ---------------------------------------------------------
  subroutine test_hartree
    type(system_t) :: sys

    PUSH_SUB(test_hartree)

    call calc_mode_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call system_init(sys)
    call poisson_test(sys%gr%mesh)
    call system_end(sys)

    POP_SUB(test_hartree)
  end subroutine test_hartree


! ---------------------------------------------------------
  subroutine test_derivatives()
    type(system_t) :: sys

    PUSH_SUB(test_derivatives)

    call system_init(sys)

    message(1) = 'Info: Testing the finite-differences derivatives.'
    message(2) = ''
    call messages_info(2)

    if(test_type == TEST_ALL .or. test_type == TEST_REAL) then
      call dderivatives_test(sys%gr%der)
    end if

    if(test_type == TEST_ALL .or. test_type == TEST_COMPLEX) then
      call zderivatives_test(sys%gr%der)
    end if

    call system_end(sys)

    POP_SUB(test_derivatives)
  end subroutine test_derivatives

  ! ---------------------------------------------------------

  subroutine test_orthogonalization()
    type(system_t) :: sys

    PUSH_SUB(test_orthogonalization)

    call calc_mode_set_parallelization(P_STRATEGY_STATES, default = .false.)
    call calc_mode_set_scalapack_compat()

    call system_init(sys)

    message(1) = 'Info: Testing orthogonalization.'
    message(2) = ''
    call messages_info(2)

    if(test_type == TEST_ALL .or. test_type == TEST_REAL) then
      message(1) = 'Info: Real wave-functions.'
      call messages_info(1)
      call dstates_calc_orth_test(sys%st, sys%mc, sys%gr%mesh)
    end if

    if(test_type == TEST_ALL .or. test_type == TEST_COMPLEX) then
      message(1) = 'Info: Complex wave-functions.'
      call messages_info(1)
      call zstates_calc_orth_test(sys%st, sys%mc, sys%gr%mesh)
    end if

    call system_end(sys)

    POP_SUB(test_orthogonalization)
  end subroutine test_orthogonalization

end program oct_test

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
