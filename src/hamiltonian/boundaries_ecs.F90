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

! This module should implement Exterior Complex Scaling.
! In principle this should work both for static and time dependent 
! calculations also in conjunction with the complex scaling module 
! (and thus allow to calculate shape-resonances energies and lifetimes)


module boundaries_ecs_m
  use io_function_m
  use io_m
  use cube_function_m
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
    ecs_t,                &
    ecs_init,             &
    ecs_end,              &
    ecs_write_info

  type ecs_t
    integer       :: type !< ECS type 
  end type ecs_t    


contains

  subroutine ecs_init(this, mesh, sb)
    type(ecs_t),              intent(out) :: this
    type(mesh_t),             intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    
    PUSH_SUB(ecs_init)


    POP_SUB(ecs_init)
  end subroutine ecs_init

  subroutine ecs_end(this)
    type(ecs_t),   intent(inout) :: this
    PUSH_SUB(ecs_end)


    POP_SUB(ecs_end)
  end subroutine ecs_end

  subroutine ecs_write_info(this)
    type(ecs_t),   intent(in) :: this
    PUSH_SUB(ecs_write_info)


    POP_SUB(ecs_write_info)
  end subroutine ecs_write_info




end module boundaries_ecs_m



