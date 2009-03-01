!! Copyright (C) 2009 X. Andrade
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
!! $Id: cube_function.F90 4907 2009-01-31 11:21:32Z xavier $

#include "global.h"

module cubic_mesh_type_m
  implicit none
  public

  type cubic_mesh_t
    private
    integer, pointer :: p
  end type cubic_mesh_t

end module cubic_mesh_type_m

module cubic_mesh_m
  use cubic_mesh_type_m
  use global_m
  use mesh_m
  use messages_m
  use profiling_m
  use simul_box_m
  
  implicit none
  private

  public ::               &
       cubic_mesh_t,      &
       cubic_mesh_init,   &
       cubic_mesh_end

  interface
    !implemented in cubic_mesh_low.c

    subroutine cubic_mesh_init_c(this, nx, ny, nz, ox, oy, oz, np, np_part, lxyz)
      use cubic_mesh_type_m
      type(cubic_mesh_t), intent(out) :: this
      integer,            intent(in)  :: nx, ny, nz, ox, oy, oz
      integer,            intent(in)  :: np, np_part, lxyz
    end subroutine cubic_mesh_init_c
    
    subroutine cubic_mesh_end(this)
      use cubic_mesh_type_m
      type(cubic_mesh_t), intent(inout) :: this
    end subroutine cubic_mesh_end

    subroutine cubic_mesh_from_mesh(this, ff)
      use cubic_mesh_type_m
      type(cubic_mesh_t), intent(inout) :: this
      FLOAT,              intent(in)    :: ff
    end subroutine cubic_mesh_from_mesh

    subroutine cubic_mesh_to_mesh(this, ff)
      use cubic_mesh_type_m
      type(cubic_mesh_t), intent(in)  :: this
      FLOAT,              intent(out) :: ff
    end subroutine cubic_mesh_to_mesh

  end interface

contains

  subroutine cubic_mesh_init(this, mesh, padding)
    type(cubic_mesh_t), intent(out) :: this
    type(mesh_t),       intent(in)  :: mesh
    integer, optional,  intent(in)  :: padding(1:MAX_DIM)

    integer :: nn(1:MAX_DIM)
    integer :: oo(1:MAX_DIM)
    integer :: pad(1:MAX_DIM)

    pad(1:MAX_DIM) = 0
    if(present(padding)) then
      pad(1:MAX_DIM) = padding(1:MAX_DIM)
    end if

    nn(1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) - mesh%idx%nr(1, 1:MAX_DIM) + 1
    oo(1:MAX_DIM) = mesh%idx%nr(1, 1:MAX_DIM)

    call cubic_mesh_init_c(this,                         &
         nn(1), nn(2), nn(3), oo(1) - pad(1), oo(2) - pad(2), oo(3) - pad(3), &
         mesh%np, mesh%np_part, mesh%idx%Lxyz(1, 1))

  end subroutine cubic_mesh_init

end module cubic_mesh_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
