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
    type(ab_t)      :: ab    !< Absorbing boundaries
    type(ecs_t)     :: ecs   !< Exterior complex scaling 
  end type bc_t    

contains

  subroutine bc_init(this, mesh, sb)
    type(bc_t),               intent(out) :: this
    type(mesh_t),             intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    
    PUSH_SUB(bc_init)


    POP_SUB(bc_init)
  end subroutine bc_init

  subroutine bc_end(this)
    type(bc_t),   intent(inout) :: this
    PUSH_SUB(bc_end)


    POP_SUB(bc_end)
  end subroutine bc_end

  subroutine bc_write_info(this)
    type(bc_t),   intent(in) :: this
    PUSH_SUB(bc_write_info)


    POP_SUB(bc_write_info)
  end subroutine bc_write_info




end module boundaries_m



