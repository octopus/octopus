! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module td_write_base_oct_m
  use global_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                     &
    td_write_prop_t,            &
    td_write_print_header_init, &
    td_write_print_header_end

  type td_write_prop_t
    type(c_ptr)          :: handle
    type(c_ptr), pointer :: mult_handles(:)
    type(mpi_grp_t)      :: mpi_grp
    integer              :: hand_start
    integer              :: hand_end
    logical              :: write = .false.
    logical              :: resolve_states = .false.    !< Whether to resolve output by state
  end type td_write_prop_t

contains

  ! ---------------------------------------------------------
  subroutine td_write_print_header_init(out)
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(td_write_print_header_init)

    call write_iter_clear(out)
    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)
    call write_iter_string(out,'# HEADER')
    call write_iter_nl(out)

    POP_SUB(td_write_print_header_init)
  end subroutine td_write_print_header_init


  ! ---------------------------------------------------------
  subroutine td_write_print_header_end(out)
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(td_write_print_header_end)

    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)

    POP_SUB(td_write_print_header_end)
  end subroutine td_write_print_header_end

end module td_write_base_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
