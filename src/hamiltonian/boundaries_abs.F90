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

! This module should implement absorbing boundaries under the form of 
! mask-function, suitable only for the TD Shroedinger equation, and
! complex absorbing potential (CAPs), suitable also for the static case

module boundaries_abs_m
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
    ab_t,                 &  
    ab_mask_t,            &
    ab_cap_t,             &
    ab_init,              &
    ab_end,               &
    ab_write_info

  type ab_mask_t
    CMPLX, pointer          :: mf(:)     !< The mask-function on the mesh
    type(cube_function_t)   :: cf        !< The mask-function on the cube
  end type ab_mask_t

  type ab_cap_t
    CMPLX, pointer          :: mf(:)     !< The CAP on the mesh
    type(cube_function_t)   :: cf        !< The CAP on the cube
  end type ab_cap_t
  
  type ab_t
    type(ab_mask_t) :: mask
    type(ab_cap_t)  :: cap
  end type ab_t    

contains

  subroutine ab_init(this, mesh, sb)
    type(ab_t),               intent(out) :: this
    type(mesh_t),             intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    
    PUSH_SUB(ab_init)


    POP_SUB(ab_init)
  end subroutine ab_init

  subroutine ab_end(this)
    type(ab_t),   intent(inout) :: this
    PUSH_SUB(ab_end)


    POP_SUB(ab_end)
  end subroutine ab_end

  subroutine ab_write_info(this)
    type(ab_t),   intent(in) :: this
    PUSH_SUB(ab_write_info)


    POP_SUB(ab_write_info)
  end subroutine ab_write_info




end module boundaries_abs_m



