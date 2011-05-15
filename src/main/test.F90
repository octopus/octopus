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
  use pfft_m
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

  call messages_print_stress(stdout, "Test mode")
  call messages_print_var_option(stdout, "TestMode", test_mode)
  call messages_print_stress(stdout)

  call fft_all_init()
#ifdef HAVE_PFFT
  call pfft_all_init()
#endif
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
