!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! $Id: poisson_no.F90 14221 2015-06-05 16:37:56Z xavier $

#include "global.h"

module poisson_no_oct_m
  use cube_oct_m
  use geometry_oct_m
  use global_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none
  private
  public ::                  &
    poisson_no_t,           &
    poisson_no_init,        &
    poisson_no_end,         &
    poisson_no_solve

  type poisson_no_t
    !integer         :: all_nodes_comm
    FLOAT                    :: qq(MAX_DIM) !< q-point for exchange in periodic system
  end type poisson_no_t
contains

  subroutine poisson_no_init(this, mesh, cube)
    type(poisson_no_t),  intent(out)   :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube
! may need to add these later for housekeeping in no poisson case. Otherwise delete these 2 lines and 
!  type member above
!    integer,             intent(in)    :: all_nodes_comm
!    logical, optional,   intent(in)    :: init_world


    PUSH_SUB(poisson_no_init)

    !this%all_nodes_comm = all_nodes_comm
    this%qq = M_ZERO

    POP_SUB(poisson_no_init)
  end subroutine poisson_no_init

  !-----------------------------------------------------------------
  subroutine poisson_no_end(this)
    type(poisson_no_t), intent(inout) :: this

    PUSH_SUB(poisson_no_end)

! nothing to do - only integer objects in poisson_no_t

    POP_SUB(poisson_no_end)
  end subroutine poisson_no_end

  !-----------------------------------------------------------------

  subroutine poisson_no_solve(this, mesh, cube, pot, rho)
    type(poisson_no_t),             intent(in)    :: this
    type(mesh_t),                   intent(in)    :: mesh
    type(cube_t),                   intent(in)    :: cube
    FLOAT,                          intent(out)   :: pot(:)
    FLOAT,                          intent(in)    :: rho(:)

    PUSH_SUB(poisson_no_solve)

    pot(1:mesh%np) = M_ZERO 

    POP_SUB(poisson_no_solve)
  end subroutine poisson_no_solve

end module poisson_no_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
