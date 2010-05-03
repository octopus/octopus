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
  use command_line_m
  use datasets_m
  use derivatives_m
  use fft_m
  use global_m
  use io_m
  use loct_m
  use parser_m
  use messages_m
  use mpi_m
  use multicomm_m
  use poisson_m
  use profiling_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  integer :: which_test
  integer, parameter ::              &
    HARTREE_TEST       =   1,        &
    DER_TEST           =   2

  call global_init()
  call parser_init()

  conf%devel_version = .true.

  call parse_integer('DebugLevel', 0, conf%debug_level)
  if(conf%debug_level>0) then
    in_debug_mode = .true.
  else
    in_debug_mode = .false.
  end if

  !%Variable WhichTest
  !%Type integer
  !%Default hartree_test
  !%Section Utilities::oct-test
  !%Description
  !% Decides what kind of test should be performed.
  !%Option hartree_test 1
  !% Tests the various Hartree solvers.
  !%Option derivatives 2
  !% Tests the implementation of the finite difference operators, used to calculate derivatives.
  !%End
  call parse_integer('WhichTest', HARTREE_TEST, which_test)
  !if(.not.varinfo_valid_option('CalculationMode', calc_mode)) call input_error('CalculationMode')
  call datasets_init(which_test)

  call parse_integer(datasets_check('Dimensions'), 3, calc_dim)
  if( calc_dim > 3 .or. calc_dim < 1) call input_error('Dimensions')

  call io_init()
  call profiling_init()

  call messages_print_stress(stdout, "Which Test")
  call messages_print_var_option(stdout, "WhichTest", which_test)
  call messages_print_stress(stdout)

  call fft_all_init()
  call unit_system_init()

  select case(which_test)
  case(HARTREE_TEST); call test_hartree()
  case(DER_TEST); call test_derivatives()
  end select

  call fft_all_end()
  call profiling_end()
  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()

  contains

  subroutine test_hartree
    type(system_t) :: sys

    call system_init(sys)
    call poisson_test(sys%gr%mesh)
    call system_end(sys)

  end subroutine test_hartree

  subroutine test_derivatives()
    type(system_t) :: sys

    call system_init(sys)

    message(1) = 'Info: Testing the finite-differences derivatives'
    message(2) = ''
    call write_info(2)

    call dderivatives_test(sys%gr%der)
    call zderivatives_test(sys%gr%der)

    call system_end(sys)
  end subroutine test_derivatives

end program oct_test

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
