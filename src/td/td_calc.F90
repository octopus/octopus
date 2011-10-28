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
!! $Id$

#include "global.h"

module td_calc_m
  use forces_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_base_m
  use hamiltonian_m
  use lasers_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use profiling_m
  use states_calc_m
  use states_m

  implicit none

  private
  public ::       &
    td_calc_tacc, &
    td_calc_tvel

contains

! ---------------------------------------------------------
! Electronic acceleration (to calculate harmonic spectrum...)
! It is calculated as:
!
! d2<x>/dt2 = d<p>/dt + i<[H,[V_nl,x]]> =
!           = i<[V_l,p]> + i<[V_nl,p]> - E(t)N + i<[H,[V_nl,x]]>
!
! WARNING: This subroutine only works if ions are not
!          allowed to move
! ---------------------------------------------------------
subroutine td_calc_tacc(gr, geo, st, hm, acc, time)
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(inout) :: geo
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  FLOAT,               intent(in)    :: time
  FLOAT,               intent(out)   :: acc(MAX_DIM)

  FLOAT :: field(MAX_DIM), x(MAX_DIM)
  CMPLX, allocatable :: hzpsi(:,:), hhzpsi(:,:), xzpsi(:,:,:), vnl_xzpsi(:,:)
  integer  :: j, k, i, ik, ist, idim

#if defined(HAVE_MPI)
  FLOAT   :: y(MAX_DIM)
#endif

  PUSH_SUB(td_calc_tacc)

  ! The term i<[V_l,p]> + i<[V_nl,p]> may be considered as equal but opposite to the
  ! force exerted by the electrons on the ions. COMMENT: This has to be thought about.
  ! Maybe we are forgetting something....
  call total_force_calculate(gr, geo, hm%ep, st, acc)

  ! Adds the laser contribution : i<[V_laser, p]>
  ! WARNING: this ignores the possibility of non-electric td external fields.
  field = M_ZERO
  do j = 1, hm%ep%no_lasers
    call laser_electric_field(hm%ep%lasers(j), field(1:gr%sb%dim), time, CNST(0.001))
    acc(1:gr%mesh%sb%dim) = acc(1:gr%mesh%sb%dim) - st%qtot*field(1:gr%mesh%sb%dim)
  end do

  if(.not. hm%ep%non_local) then
    POP_SUB(td_calc_tacc)
    return
  end if

  ! And now, i<[H,[V_nl,x]]>
  x = M_ZERO
  SAFE_ALLOCATE(hzpsi (1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(hhzpsi(1:3, 1:gr%mesh%np))

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      call zhamiltonian_apply(hm, gr%der, st%zpsi(:, :, ist, ik), hzpsi(:,:), ist, ik, time)

      SAFE_ALLOCATE(xzpsi    (1:gr%mesh%np, 1:st%d%dim, 1:3))
      SAFE_ALLOCATE(vnl_xzpsi(1:gr%mesh%np, 1:st%d%dim))
      xzpsi = M_z0
      do k = 1, gr%mesh%np
        do j = 1, gr%mesh%sb%dim
          xzpsi(k, 1:st%d%dim, j) = gr%mesh%x(k, j)*st%zpsi(k, 1:st%d%dim, ist, ik)
        end do
      end do

      do j = 1, gr%mesh%sb%dim
        call zhamiltonian_apply(hm, gr%der, xzpsi(:, :, j), vnl_xzpsi, ist, ik, time, terms = TERM_NON_LOCAL_POTENTIAL)

        do idim = 1, st%d%dim
          x(j) = x(j) - 2*st%occ(ist, ik)*zmf_dotp(gr%mesh, hzpsi(1:gr%mesh%np, idim), vnl_xzpsi(:, idim) )
        end do
      end do

      xzpsi = M_z0
      do k = 1, gr%mesh%np
        do j = 1, gr%mesh%sb%dim
          xzpsi(k, 1:st%d%dim, j) = gr%mesh%x(k, j)*hzpsi(k, 1:st%d%dim)
        end do
      end do

      do j = 1, gr%mesh%sb%dim
        call zhamiltonian_apply(hm, gr%der, xzpsi(:, :, j), vnl_xzpsi, ist, ik, time, terms = TERM_NON_LOCAL_POTENTIAL)

        do idim = 1, st%d%dim
          x(j) = x(j) + 2*st%occ(ist, ik)* &
            zmf_dotp(gr%mesh, st%zpsi(1:gr%mesh%np, idim, ist, ik), vnl_xzpsi(:, idim) )
        end do
      end do
      SAFE_DEALLOCATE_A(xzpsi)
      SAFE_DEALLOCATE_A(vnl_xzpsi)

    end do
  end do
  SAFE_DEALLOCATE_A(hzpsi)
  SAFE_DEALLOCATE_A(hhzpsi)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(x(1), y(1), gr%mesh%sb%dim, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    x = y
  end if
#endif
  acc = acc + x

  POP_SUB(td_calc_tacc)
end subroutine td_calc_tacc

! ---------------------------------------------------------
! Electronic velocity (to calculate harmonic spectrum...)
! It is calculated as:
!
! d<x>/dt = <p>
! ---------------------------------------------------------
subroutine td_calc_tvel(gr, geo, st, hm, vel, time)
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(inout) :: geo
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  FLOAT,               intent(in)    :: time
  FLOAT,               intent(out)   :: vel(MAX_DIM)

  FLOAT, allocatable :: momentum(:,:,:)
  
  PUSH_SUB(td_calc_tvel)

  SAFE_ALLOCATE(momentum(1:gr%mesh%sb%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
  call states_calc_momentum(st, gr%der, momentum)

  momentum(1:gr%mesh%sb%dim, st%st_start:st%st_end, 1) = & 
    sum(momentum(1:gr%mesh%sb%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end), 3)
  momentum(1:gr%mesh%sb%dim, 1, 1) = & 
    sum(momentum(1:gr%mesh%sb%dim, st%st_start:st%st_end, 1), 2)
  vel = momentum(1:gr%mesh%sb%dim, 1, 1)

  SAFE_DEALLOCATE_A(momentum)
  POP_SUB(td_calc_tvel)
end subroutine td_calc_tvel

end module td_calc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
