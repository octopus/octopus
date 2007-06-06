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
  use global_m
  use grid_m
  use mesh_m
  use messages_m
  use simul_box_m
  use ps_m
  use specie_m
  use specie_pot_m
  use submesh_m
  use geometry_m
  use mpi_m
  use mpi_debug_m

  implicit none

  private
  public :: &
       hgh_projector_t,              &
       hgh_projector_null,           &
       hgh_projector_init,           &
       dhgh_project, zhgh_project,   &
       dhgh_dproject, zhgh_dproject, &
       hgh_projector_end

  

  type hgh_projector_t
    private
    integer        :: n_s         ! number of points inside the sphere
    FLOAT, pointer :: p(:, :)     ! projectors
    FLOAT, pointer :: dp(:, :, :) ! projectors derivatives
    FLOAT, pointer :: lp(:, :, :) ! angular momentum times projectors
    FLOAT          :: h(3, 3)     ! parameters
    FLOAT          :: k(3, 3)     ! spin-orbit parameters
  end type hgh_projector_t


contains

  ! ---------------------------------------------------------
  subroutine hgh_projector_null(hgh_p)
    type(hgh_projector_t), intent(out) :: hgh_p

    call push_sub('hgh_projector.hgh_projector_null')

    nullify(hgh_p%p)
    nullify(hgh_p%dp)
    nullify(hgh_p%lp)
    hgh_p%h = M_ZERO
    hgh_p%k = M_ZERO

    call pop_sub()
  end subroutine hgh_projector_null

  ! ---------------------------------------------------------
  subroutine hgh_projector_init(hgh_p, sm, gr, a, l, lm)
    type(hgh_projector_t), intent(inout) :: hgh_p
    type(submesh_t),       intent(in)    :: sm
    type(grid_t),          intent(in)    :: gr
    type(atom_t),          intent(in)    :: a
    integer,               intent(in)    :: l, lm
 
    integer :: j, k, i
    FLOAT :: v, dv(3), x(3), x_in(3), r

    call push_sub('hgh_projector.hgh_projector_init')

    hgh_p%n_s = sm%ns
    ALLOCATE(hgh_p%p (hgh_p%n_s, 3),    hgh_p%n_s*3)
    ALLOCATE(hgh_p%dp(hgh_p%n_s, 3, 3), hgh_p%n_s*3*3)
    ALLOCATE(hgh_p%lp(hgh_p%n_s, 3, 3), hgh_p%n_s*3*3)

    do j = 1, hgh_p%n_s

      do k = 1, 3**gr%sb%periodic_dim
        x_in(:) = gr%m%x(sm%jxyz(j), :) - gr%sb%shift(k,:)
        x(:) = x_in(:) - a%x
        r = sqrt(sum(x*x))
        if (r > a%spec%ps%rc_max + gr%m%h(1)) cycle

        do i = 1, 3
          call specie_real_nl_projector(a%spec, x, l, lm, i, v, dv(1:3))
          hgh_p%p(j, i) = v
          hgh_p%dp(j, :, i) = dv(:)
          hgh_p%lp(j, 1, i) = x(2)*dv(3) - x(3)*dv(2)
          hgh_p%lp(j, 2, i) = x(3)*dv(1) - x(1)*dv(3)
          hgh_p%lp(j, 3, i) = x(1)*dv(2) - x(2)*dv(1)
        end do
      end do
    end do

    hgh_p%h(:, :) = a%spec%ps%h(l, :, :)
    hgh_p%k(:, :) = a%spec%ps%k(l, :, :)

    call pop_sub()
  end subroutine hgh_projector_init

  ! ---------------------------------------------------------
  subroutine hgh_projector_end(hgh_p)
    type(hgh_projector_t), intent(inout) :: hgh_p

    call push_sub('hgh_projector.hgh_projector_end')

    if (associated(hgh_p%p))      deallocate(hgh_p%p)
    if (associated(hgh_p%dp))     deallocate(hgh_p%dp)
    if (associated(hgh_p%lp))   deallocate(hgh_p%lp)

    call pop_sub()
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
