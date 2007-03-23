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
  use global_m
  use grid_m
  use mesh_m
  use messages_m
  use simul_box_m
  use ps_m
  use specie_m
  use specie_pot_m
  use geometry_m
  use mpi_m
  use mpi_debug_m

  implicit none

  private
  public :: &
       kb_projector_t,             &
       kb_projector_null,          &
       kb_projector_init,          &
       dkb_project, zkb_project,   &
       dkb_dproject, zkb_dproject, &
       kb_projector_end

  type kb_projector_t
    private
    integer          :: n_s       ! number of points inside the sphere
    integer          :: n_c       ! number of components per projector
    FLOAT,   pointer :: p(:,:)    ! projectors
    FLOAT,   pointer :: dp(:,:,:) ! projectors derivatives
    FLOAT,   pointer :: e(:)      ! KB energies
  end type kb_projector_t


contains

  ! ---------------------------------------------------------
  subroutine kb_projector_null(kb_p)
    type(kb_projector_t), intent(out) :: kb_p

    call push_sub('kb_projector.kb_projector_null')

    nullify(kb_p%p)
    nullify(kb_p%dp)
    nullify(kb_p%e)

    call pop_sub()
  end subroutine kb_projector_null

  ! ---------------------------------------------------------
  subroutine kb_projector_init(kb_p, n_s, jxyz, gr, a, l, lm)
    type(kb_projector_t), intent(inout) :: kb_p
    integer,              intent(in)    :: n_s
    integer,              intent(in)    :: jxyz(:)
    type(grid_t),         intent(in)    :: gr
    type(atom_t),         intent(in)    :: a
    integer,              intent(in)    :: l, lm
 
    integer :: n_c, j, k, i
    FLOAT :: v, dv(3), r, x(3), x_in(3)

    call push_sub('kb_projector.kb_projector_init')

    kb_p%n_s = n_s
    if (l == 0 .or. a%spec%ps%kbc == 1) then
      n_c = 1
    else ! we have j-dependent projectors
      n_c = 2
    end if
    kb_p%n_c = n_c
    ALLOCATE(kb_p%p (n_s, n_c),    n_s*n_c)
    ALLOCATE(kb_p%dp(n_s, 3, n_c), n_s*3*n_c)
    ALLOCATE(kb_p%e (n_c),         n_c)
    
    if (gr%sb%periodic_dim == 0) then 

      do i = 1, n_c
        call double_grid_apply_non_local(gr%dgrid, a%spec, gr%m, a%x, n_s, jxyz, l, lm, i, &
             kb_p%p(:, i), kb_p%dp(:,:,i))
      end do

    else 

      do j = 1, kb_p%n_s
        do k = 1, 3**gr%sb%periodic_dim
          x_in(:) = gr%m%x(jxyz(j), :) - gr%sb%shift(k,:)
          x(:) = x_in(:) - a%x
          r = sqrt(sum(x*x))
          if (r > a%spec%ps%rc_max + gr%m%h(1)) cycle
          
          do i = 1, n_c
            call specie_real_nl_projector(a%spec, a%x, x_in, l, lm, i, v, dv)
            kb_p%p(j, i) = v
            kb_p%dp(j, :, i) = dv
          end do
        end do
      end do

    end if

    do i = 1, n_c
      kb_p%e(i) = a%spec%ps%h(l, i, i)
    end do

    if (n_c == 2) then
      ! We need to weight the projectors.
      ! The weights will be included in the KB energies
      kb_p%e(1) = kb_p%e(1)*real(l+1, REAL_PRECISION)/real(2*l+1, REAL_PRECISION)
      kb_p%e(2) = kb_p%e(2)*real(l,   REAL_PRECISION)/real(2*l+1, REAL_PRECISION)
    end if

    call pop_sub()
  end subroutine kb_projector_init

  ! ---------------------------------------------------------
  subroutine kb_projector_end(kb_p)
    type(kb_projector_t), intent(inout) :: kb_p

    call push_sub('kb_projector.kb_projector_end')

    if (associated(kb_p%p))  deallocate(kb_p%p)
    if (associated(kb_p%dp)) deallocate(kb_p%dp)
    if (associated(kb_p%e))  deallocate(kb_p%e)

    call pop_sub()
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
