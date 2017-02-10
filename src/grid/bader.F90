!! Copyright (C) 2017 J.Jornet-Somoza
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

module bader_oct_m
  use global_oct_m
  use index_oct_m
  use mesh_oct_m
  use messages_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public :: &
     bader_t,       &
     bader_init,    &
     bader_end,     &
     bader_analyze, &
     bader_write

  type bader_t
    integer, allocatable :: path(:,:)  
    integer, allocatable :: volnum(:,:,:), known
    integer, allocatable :: nnion(:)  
    integer          :: nvols, pnum, bnum, pdim, bdim, refine_edge_itrs
    integer, pointer :: position(:)
    FLOAT,   pointer :: val(:)
    FLOAT,   pointer :: volume(:)

    FLOAT,   pointer :: population(:)
  end type bader_t

contains

  !----------------------------------------------------------------
  subroutine bader_init(this, mesh)
    type(bader_t), intent(out) :: this
    type(mesh_t),   intent(in)  :: mesh
    
    PUSH_SUB(bader_init)

    if(mesh%parallel_in_domains) &
      call messages_experimental("Bader bader parallel in domains")
    
    SAFE_ALLOCATE(this%map(1:mesh%np))
    this%map(1:mesh%np) = -1

    POP_SUB(bader_init)
  end subroutine bader_init


  !----------------------------------------------------------------
  subroutine bader_end(this)
    type(bader_t), intent(inout) :: this

    PUSH_SUB(bader_end)

    ASSERT(associated(this%map))
    SAFE_DEALLOCATE_P(this%map)

    if(associated(this%position)) then
      SAFE_DEALLOCATE_P(this%position)
      SAFE_DEALLOCATE_P(this%val)
      SAFE_DEALLOCATE_P(this%volume)
      SAFE_DEALLOCATE_P(this%population)
    end if
    
    POP_SUB(bader_end)
  end subroutine bader_end


  !----------------------------------------------------------------
  subroutine bader_analyze(this, mesh, f, rho, threshold)
    type(bader_t), intent(inout) :: this
