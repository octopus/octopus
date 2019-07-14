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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module hgh_projector_oct_m
  use atom_oct_m
  use blas_oct_m
  use comm_oct_m
  use global_oct_m
  use hardware_oct_m
  use mesh_oct_m
  use messages_oct_m
  use profiling_oct_m
  use ps_oct_m
  use species_oct_m
  use submesh_oct_m

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
    integer                :: n_s         !< number of points inside the sphere
    FLOAT, pointer, public :: dp(:, :)    !< projectors
    CMPLX, pointer         :: zp(:, :)
    FLOAT,          public :: h(3, 3)     !< parameters
    FLOAT                  :: k(3, 3)     !< spin-orbit parameters
  end type hgh_projector_t


contains

  ! ---------------------------------------------------------
  subroutine hgh_projector_null(hgh_p)
    type(hgh_projector_t), intent(out) :: hgh_p

    PUSH_SUB(hgh_projector_null)

    nullify(hgh_p%dp)
    nullify(hgh_p%zp)
    hgh_p%h = M_ZERO
    hgh_p%k = M_ZERO

    POP_SUB(hgh_projector_null)
  end subroutine hgh_projector_null

  ! ---------------------------------------------------------
  subroutine hgh_projector_init(hgh_p, sm, reltyp, a, l, lm, so_strength)
    type(hgh_projector_t), intent(inout) :: hgh_p
    type(submesh_t),       intent(in)    :: sm
    integer,               intent(in)    :: reltyp
    type(atom_t), target,  intent(in)    :: a
    integer,               intent(in)    :: l, lm
    FLOAT,                 intent(in)    :: so_strength

    integer :: is, i
    FLOAT :: v, x(1:3)
    type(ps_t), pointer :: ps

    PUSH_SUB(hgh_projector_init)

    hgh_p%n_s = sm%np
    if(reltyp == 0) then
      SAFE_ALLOCATE(hgh_p%dp(1:hgh_p%n_s, 1:3))
      x = M_ZERO
      do is = 1, hgh_p%n_s
        x(1:3) = sm%x(is, 1:3)

        do i = 1, 3
          call species_real_nl_projector(a%species, x, l, lm, i, v)
          hgh_p%dp (is, i) = v
        end do
      end do
    else
      SAFE_ALLOCATE(hgh_p%zp(1:hgh_p%n_s, 1:3))
      do i = 1, 3
        call species_nl_projector(a%species, hgh_p%n_s, sm%x(:, 0:3), l, lm, i, hgh_p%zp(:, i)) 
      end do
    end if

    ps => species_ps(a%species)
    hgh_p%h(:, :) = ps%h(l, :, :)
    hgh_p%k(:, :) = ps%k(l, :, :)*so_strength
    nullify(ps)

    POP_SUB(hgh_projector_init)
  end subroutine hgh_projector_init

  ! ---------------------------------------------------------
  subroutine hgh_projector_end(hgh_p)
    type(hgh_projector_t), intent(inout) :: hgh_p

    PUSH_SUB(hgh_projector_end)

    SAFE_DEALLOCATE_P(hgh_p%dp)
    SAFE_DEALLOCATE_P(hgh_p%zp)

    POP_SUB(hgh_projector_end)
  end subroutine hgh_projector_end


#include "undef.F90"
#include "real.F90"
#include "hgh_projector_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hgh_projector_inc.F90"

end module hgh_projector_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
