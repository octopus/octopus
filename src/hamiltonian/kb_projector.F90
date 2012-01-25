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

module kb_projector_m
  use double_grid_m
  use geometry_m
  use global_m
  use grid_m
  use lalg_basic_m
  use mesh_m
  use messages_m
  use profiling_m
  use ps_m
  use mpi_m
  use multicomm_m
  use simul_box_m
  use species_m
  use submesh_m

  implicit none

  private
  public :: &
       kb_projector_t,             &
       kb_projector_null,          &
       kb_projector_init,          &
       dkb_project, zkb_project,   &
       dkb_project_bra,            &
       zkb_project_bra,            &
       dkb_project_ket,            &
       zkb_project_ket,            &
       dkb_mul_energies,           &
       zkb_mul_energies,           &
       kb_projector_end

  type kb_projector_t
    integer          :: n_s       !< number of points inside the sphere
    integer          :: n_c       !< number of components per projector
    FLOAT,   pointer :: p(:,:)    !< projectors
    FLOAT            :: e(2)      !< KB energies
  end type kb_projector_t


contains

  ! ---------------------------------------------------------
  subroutine kb_projector_null(kb_p)
    type(kb_projector_t), intent(out) :: kb_p

    PUSH_SUB(kb_projector_null)

    nullify(kb_p%p)

    POP_SUB(kb_projector_null)
  end subroutine kb_projector_null

  ! ---------------------------------------------------------
  subroutine kb_projector_init(kb_p, sm, gr, a, l, lm)
    type(kb_projector_t), intent(inout) :: kb_p
    type(submesh_t),      intent(in)    :: sm
    type(grid_t),         intent(in)    :: gr
    type(atom_t),         intent(in)    :: a
    integer,              intent(in)    :: l, lm

    integer :: n_c, ic
    type(ps_t), pointer :: ps

    PUSH_SUB(kb_projector_init)

    ps => species_ps(a%spec)

    kb_p%n_s = sm%np
    if (l == 0 .or. ps%kbc == 1) then
      n_c = 1
    else ! we have j-dependent projectors
      n_c = 2
    end if
    kb_p%n_c = n_c

    SAFE_ALLOCATE(kb_p%p (1:kb_p%n_s, 1:2))
    kb_p%p = M_ZERO
    kb_p%e = M_ZERO
    
    do ic = 1, n_c
      call double_grid_apply_non_local(gr%dgrid, a%spec, sm%mesh, sm, a%x, kb_p%p(:, ic), l, lm, ic)
      kb_p%e(ic) = ps%h(l, ic, ic)
    end do

    if (n_c == 2) then
      ! We need to weight the projectors.
      ! The weights will be included in the KB energies
      kb_p%e(1) = kb_p%e(1)*real(l+1, REAL_PRECISION)/real(2*l+1, REAL_PRECISION)
      kb_p%e(2) = kb_p%e(2)*real(l,   REAL_PRECISION)/real(2*l+1, REAL_PRECISION)
    end if

    nullify(ps)
    POP_SUB(kb_projector_init)
  end subroutine kb_projector_init

  ! ---------------------------------------------------------
  subroutine kb_projector_end(kb_p)
    type(kb_projector_t), intent(inout) :: kb_p

    PUSH_SUB(kb_projector_end)

    SAFE_DEALLOCATE_P(kb_p%p)

    POP_SUB(kb_projector_end)
  end subroutine kb_projector_end

#include "undef.F90"
#include "real.F90"
#include "kb_projector_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "kb_projector_inc.F90"

end module kb_projector_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
