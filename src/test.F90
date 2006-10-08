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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: $

#include "global.h"

program oct_test
  use string_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use profiling_m
  use varinfo_m
  use fft_m
  use units_m
  use mpi_m
  use multicomm_m
  use system_m
  use poisson_m

  implicit none

  integer :: which_test
  integer, parameter :: HARTREE_TEST = 1

  call global_init()
  call parser_init()

  conf%devel_version = .true.

  call loct_parse_int('DebugLevel', 0, conf%debug_level)
  if(conf%debug_level>0) then
    in_debug_mode = .true.
  else
    in_debug_mode = .false.
  end if

  !%Variable WhichTest
  !%Type integer
  !%Default hartree_test
  !%Section Generalities
  !%Description
  !% Decides what kind of test should be performed.
  !%Option hartree_test 1
  !% Tests the various Hartree solvers.
  !%End
  call loct_parse_int('WhichTest', HARTREE_TEST, which_test)
  !if(.not.varinfo_valid_option('CalculationMode', calc_mode)) call input_error('CalculationMode')
  call datasets_init(which_test)

  call loct_parse_int(check_inp('Dimensions'), 3, calc_dim)
  if( calc_dim > 3 .or. calc_dim < 1) call input_error('Dimensions')

  call io_init()

  call messages_print_stress(stdout, "Which Test")
  call messages_print_var_option(stdout, "WhichTest", which_test)
  call messages_print_stress(stdout)

#ifdef HAVE_FFT
  call fft_all_init()
#endif
  call units_init()

  select case(which_test)
  case(HARTREE_TEST); call test_hartree
  end select

  call fft_all_end()
  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()

  contains

  subroutine test_hartree
    integer :: parallel_mask
    type(system_t) :: sys

    parallel_mask = 0
    parallel_mask = ibset(parallel_mask, P_STRATEGY_DOMAINS - 1) ! all modes are parallel in domains

    call system_init(sys, parallel_mask)
    call poisson_test(sys%gr)
    call system_end(sys)

  end subroutine test_hartree

end program oct_test
