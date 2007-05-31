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
  use submesh_m
  use geometry_m
  use mpi_m
  use mpi_debug_m
  use multicomm_m

  implicit none

  private
  public :: &
       kb_projector_t,             &
       kb_projector_null,          &
       kb_projector_init,          &
#ifdef HAVE_MPI
       kb_projector_broadcast,     &
#endif
       dkb_project, zkb_project,   &
       dkb_dproject, zkb_dproject, &
       kb_projector_end

  type kb_projector_t
    private
    integer          :: n_s       ! number of points inside the sphere
    integer          :: n_c       ! number of components per projector
    FLOAT,   pointer :: p(:,:)    ! projectors
    FLOAT,   pointer :: dp(:,:,:) ! projectors derivatives
    FLOAT            :: e(2)      ! KB energies
  end type kb_projector_t


contains

  ! ---------------------------------------------------------
  subroutine kb_projector_null(kb_p)
    type(kb_projector_t), intent(out) :: kb_p

    call push_sub('kb_projector.kb_projector_null')

    nullify(kb_p%p)
    nullify(kb_p%dp)

    call pop_sub()
  end subroutine kb_projector_null

  ! ---------------------------------------------------------
  subroutine kb_projector_init(kb_p, sm, gr, a, l, lm, gen_grads)
    type(kb_projector_t), intent(inout) :: kb_p
    type(submesh_t),      intent(in)    :: sm
    type(grid_t),         intent(in)    :: gr
    type(atom_t),         intent(in)    :: a
    integer,              intent(in)    :: l, lm
    logical,              intent(in)    :: gen_grads

    integer :: n_c, j, k, i
    FLOAT :: v, dv(3), r, x(3), x_in(3)

    call push_sub('kb_projector.kb_projector_init')

    kb_p%n_s = sm%ns
    if (l == 0 .or. a%spec%ps%kbc == 1) then
      n_c = 1
    else ! we have j-dependent projectors
      n_c = 2
    end if
    kb_p%n_c = n_c

    ALLOCATE(kb_p%p (kb_p%n_s, n_c),    kb_p%n_s*n_c)
    ALLOCATE(kb_p%dp(kb_p%n_s, 3, n_c), kb_p%n_s*3*n_c)
    kb_p%p = M_ZERO
    kb_p%dp = M_ZERO
    
    if (gr%sb%periodic_dim == 0) then 

      do i = 1, n_c
        call double_grid_apply_non_local(gr%dgrid, a%spec, gr%m, a%x, sm, l, lm, i, kb_p%p(:, i))
        if (gen_grads) then
          call double_grid_apply_gnon_local(gr%dgrid, a%spec, gr%m, a%x, sm, l, lm, i, kb_p%dp(:,:,i))
        end if
      end do

    else 

      do j = 1, kb_p%n_s
        do k = 1, 3**gr%sb%periodic_dim
          x_in(:) = gr%m%x(sm%jxyz(j), :) - gr%sb%shift(k,:)
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
#ifdef HAVE_MPI
  subroutine kb_projector_broadcast(kb_p, sm, gr, mc, a, l, lm, gen_grads, root)
    type(kb_projector_t), intent(inout) :: kb_p
    type(submesh_t),      intent(in)    :: sm
    type(grid_t),         intent(in)    :: gr
    type(multicomm_t),    intent(in)    :: mc
    type(atom_t),         intent(in)    :: a
    integer,              intent(in)    :: l, lm
    logical,              intent(in)    :: gen_grads
    integer,              intent(in)    :: root

    integer :: n_c, i, rank, mpi_err

    call push_sub('kb_projector.kb_projector_init')

    rank = mc%who_am_i(P_STRATEGY_STATES) - 1

    if (root /= rank ) then
      
      kb_p%n_s = sm%ns
      if (l == 0 .or. a%spec%ps%kbc == 1) then
        n_c = 1
      else ! we have j-dependent projectors
        n_c = 2
      end if
      kb_p%n_c = n_c
      ALLOCATE(kb_p%p (kb_p%n_s, n_c),    kb_p%n_s*n_c)
      ALLOCATE(kb_p%dp(kb_p%n_s, 3, n_c), kb_p%n_s*3*n_c)
      
      do i = 1, n_c
        kb_p%e(i) = a%spec%ps%h(l, i, i)
      end do
      
      if (n_c == 2) then
        ! We need to weight the projectors.
        ! The weights will be included in the KB energies
        kb_p%e(1) = kb_p%e(1)*real(l+1, REAL_PRECISION)/real(2*l+1, REAL_PRECISION)
        kb_p%e(2) = kb_p%e(2)*real(l,   REAL_PRECISION)/real(2*l+1, REAL_PRECISION)
      end if

    end if

    call MPI_Bcast(kb_p%p, kb_p%n_s*kb_p%n_c, MPI_FLOAT, root, mc%group_comm(P_STRATEGY_STATES), mpi_err)

    if (gen_grads) then 
      call MPI_Bcast(kb_p%dp, kb_p%n_s*3*kb_p%n_c, MPI_FLOAT, root, mc%group_comm(P_STRATEGY_STATES), mpi_err)
    end if

    call pop_sub()
  end subroutine kb_projector_broadcast
#endif
  ! ---------------------------------------------------------
  subroutine kb_projector_end(kb_p)
    type(kb_projector_t), intent(inout) :: kb_p

    call push_sub('kb_projector.kb_projector_end')

    if (associated(kb_p%p))  deallocate(kb_p%p)
    if (associated(kb_p%dp)) deallocate(kb_p%dp)

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
