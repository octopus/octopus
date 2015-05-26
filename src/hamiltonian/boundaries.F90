!! Copyright (C) 2015 U. De Giovannini
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

#include "global.h"

module boundaries_m
  use boundaries_abs_m
  use boundaries_ecs_m
  use io_m
  use io_function_m
  use io_m
  use cube_function_m
  use geometry_m
  use global_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use simul_box_m
  use varinfo_m

  implicit none

  private
  public ::               &
    bc_init,              &
    bc_end,               &
    bc_write_info,        &
    bc_t

  type bc_t
    logical         :: use_ab
    type(ab_t)      :: ab    !< Absorbing boundaries (mask & cap)
    logical         :: use_ecs
    type(ecs_t)     :: ecs   !< Exterior complex scaling 
  end type bc_t    

  integer, parameter :: &
    BC_NONE = 0, &
    BC_AB   = 2, &
    BC_ECS  = 4

contains

  ! ---------------------------------------------------------
  subroutine bc_init(this, mesh, sb, geo)
    type(bc_t),               intent(out) :: this
    type(mesh_t),             intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    type(geometry_t),         intent(in)  :: geo

    integer :: bc_flags
    
    PUSH_SUB(bc_init)

    !%Variable Boundaries
    !%Type flag
    !%Default none
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% To improve the quality of the spectra by avoiding the formation of
    !% standing density waves, one can make the boundaries of the simulation
    !% box absorbing and use exterior complex scaling.
    !%Option none 0
    !% Reflecting boundaries.
    !%Option absorbing 2
    !% Absorbing boundaries with a mask function or a complex absorbing potential.
    !%Option exterior 4
    !% Exterior complex scaling.
    !%End
    call parse_variable('Boundaries', BC_NONE, bc_flags)
    if(.not.varinfo_valid_option('Boundaries', bc_flags, is_flag = .true.)) then
      call messages_input_error('Boundaries')
    end if

    this%use_ab  = iand(bc_flags, BC_AB)  /= 0
    this%use_ecs = iand(bc_flags, BC_ECS) /= 0

    if(this%use_ab) call ab_init(this%ab, mesh, sb, geo)
    if(this%use_ecs) then
      call messages_not_implemented('Exterior complex scaling')
    end if

    POP_SUB(bc_init)
  end subroutine bc_init

  ! ---------------------------------------------------------
  subroutine bc_end(this)
    type(bc_t),   intent(inout) :: this
    PUSH_SUB(bc_end)


    POP_SUB(bc_end)
  end subroutine bc_end

  ! ---------------------------------------------------------
  subroutine bc_write_info(this)
    type(bc_t),   intent(in) :: this
    PUSH_SUB(bc_write_info)


    POP_SUB(bc_write_info)
  end subroutine bc_write_info

end module boundaries_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
