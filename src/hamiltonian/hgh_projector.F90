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
!! $Id: epot.F90 2648 2007-01-09 19:08:10Z lorenzen $

#include "global.h"

module hgh_projector_m
  use atom_m
  use global_m
  use grid_m
  use lalg_basic_m
  use mesh_m
  use messages_m
  use simul_box_m
  use profiling_m
  use ps_m
  use species_m
  use submesh_m
  use geometry_m
  use mpi_m

  implicit none

  private
  public :: &
       hgh_projector_t,              &
       hgh_projector_null,           &
       hgh_projector_init,           &
       dhgh_project, zhgh_project,   &
       dhgh_project_bra,             &
       zhgh_project_bra,             &
       dhgh_project_ket,             &
       zhgh_project_ket,             &
       hgh_projector_end

  

  type hgh_projector_t
    private
    integer        :: n_s         !< number of points inside the sphere
    FLOAT, pointer :: p(:, :)     !< projectors
    FLOAT, pointer :: lp(:, :, :) !< angular momentum times projectors
    FLOAT          :: h(3, 3)     !< parameters
    FLOAT          :: k(3, 3)     !< spin-orbit parameters
  end type hgh_projector_t


contains

  ! ---------------------------------------------------------
  subroutine hgh_projector_null(hgh_p)
    type(hgh_projector_t), intent(out) :: hgh_p

    PUSH_SUB(hgh_projector_null)

    nullify(hgh_p%p)
    nullify(hgh_p%lp)
    hgh_p%h = M_ZERO
    hgh_p%k = M_ZERO

    POP_SUB(hgh_projector_null)
  end subroutine hgh_projector_null

  ! ---------------------------------------------------------
  subroutine hgh_projector_init(hgh_p, sm, gr, a, l, lm, so_strength)
    type(hgh_projector_t), intent(inout) :: hgh_p
    type(submesh_t),       intent(in)    :: sm
    type(grid_t),          intent(in)    :: gr
    type(atom_t),          intent(in)    :: a
    integer,               intent(in)    :: l, lm
    FLOAT,                 intent(in)    :: so_strength

    integer :: is, i
    FLOAT :: v, dv(1:3), x(1:3)
    type(ps_t), pointer :: ps

    PUSH_SUB(hgh_projector_init)

    hgh_p%n_s = sm%np
    SAFE_ALLOCATE(hgh_p%p (1:hgh_p%n_s, 1:3))
    SAFE_ALLOCATE(hgh_p%lp(1:hgh_p%n_s, 1:3, 1:3))

    x = M_ZERO
    do is = 1, hgh_p%n_s
      x(1:3) = sm%x(is, 1:3)
      
      do i = 1, 3
        call species_real_nl_projector(a%spec, x, l, lm, i, v, dv)
        hgh_p%p (is, i) = v
        hgh_p%lp(is, 1, i) = x(2)*dv(3) - x(3)*dv(2)
        hgh_p%lp(is, 2, i) = x(3)*dv(1) - x(1)*dv(3)
        hgh_p%lp(is, 3, i) = x(1)*dv(2) - x(2)*dv(1)
      end do
    end do

    ps => species_ps(a%spec)
    hgh_p%h(:, :) = ps%h(l, :, :)
    hgh_p%k(:, :) = ps%k(l, :, :)*so_strength
    nullify(ps)

    POP_SUB(hgh_projector_init)
  end subroutine hgh_projector_init

  ! ---------------------------------------------------------
  subroutine hgh_projector_end(hgh_p)
    type(hgh_projector_t), intent(inout) :: hgh_p

    PUSH_SUB(hgh_projector_end)

    SAFE_DEALLOCATE_P(hgh_p%p)
    SAFE_DEALLOCATE_P(hgh_p%lp)

    POP_SUB(hgh_projector_end)
  end subroutine hgh_projector_end

#include "undef.F90"
#include "real.F90"
#include "hgh_projector_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hgh_projector_inc.F90"

end module hgh_projector_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
