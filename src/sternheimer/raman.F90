!! Copyright (C) 2004 Xavier Andrade
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
!! $Id: em_resp.F90 2686 2007-02-03 22:10:51Z xavier $

#include "global.h"
#define RESTART_DIR "raman/"

module raman_m
  use datasets_m
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use geometry_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use loct_parser_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use h_sys_output_m
  use phonons_lr_m
  use pert_m
  use restart_m
  use states_m
  use sternheimer_m
  use string_m
  use system_m
  use units_m

  implicit none

  private
  public :: &
       raman_run
  
contains

  ! ---------------------------------------------------------
  subroutine raman_run(sys, h, fromscratch)
    type(system_t),         intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(in)    :: fromscratch

    type(sternheimer_t) :: sh
    type(lr_t)          :: lr(1:1)

    !CONSTRUCT

    call push_sub('raman.raman')

    call read_wfs(sys%st, sys%gr, sys%geo, .false.)

    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)

    call system_h_setup(sys, h)
    call sternheimer_init(sh, sys, h, "VM")

    call lr_init(lr(1))
    call lr_allocate(lr(1), sys%st, sys%gr%m)

    message(1) = "HELLO WORLD!"
    call write_info(1)
    
    call lr_dealloc(lr(1))
    call sternheimer_end(sh)
    call states_deallocate_wfns(sys%st)

    call pop_sub()
  end subroutine raman_run

end module raman_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
