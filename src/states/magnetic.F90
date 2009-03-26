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
!! $Id: states.F90 2515 2006-10-24 17:13:30Z acastro $

#include "global.h"

module magnetic_m
  use derivatives_m
  use geometry_m
  use global_m
  use grid_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use profiling_m
  use states_m
  use states_dim_m
  use poisson_m

  implicit none

  private
  public ::                    &
    magnetic_density,          &
    magnetic_moment,           &
    magnetic_local_moments,    &
    calc_physical_current,     &
    magnetic_induced

contains

  ! ---------------------------------------------------------
  subroutine magnetic_density(m, st, rho, md)
    type(mesh_t),   intent(in)  :: m
    type(states_t), intent(in)  :: st
    FLOAT,             intent(in)  :: rho(:,:) ! (np, st%d%nspin)
    FLOAT,             intent(out) :: md(:,:)   ! (np, 3)

    call push_sub('states.magnetic_density')

    select case (st%d%ispin)
    case (UNPOLARIZED)
      md = M_ZERO

    case (SPIN_POLARIZED)
      md = M_ZERO
      md(1:m%np, 3) = rho(1:m%np, 1) - rho(1:m%np, 2)

    case (SPINORS)
      md(1:m%np, 1) =  M_TWO*rho(1:m%np, 3)
      md(1:m%np, 2) = -M_TWO*rho(1:m%np, 4)
      md(1:m%np, 3) = rho(1:m%np, 1) - rho(1:m%np, 2)
    end select

    call pop_sub()
  end subroutine magnetic_density


  ! ---------------------------------------------------------
  subroutine magnetic_moment(m, st, rho, mm)
    type(mesh_t),   intent(in)  :: m
    type(states_t), intent(in)  :: st
    FLOAT,          intent(in)  :: rho(:,:) ! (m%np_part, st%d%nspin)
    FLOAT,          intent(out) :: mm(3)

    FLOAT, allocatable :: md(:,:)

    call push_sub('states.states_magnetic_moment')

    ALLOCATE(md(m%np, MAX_DIM), m%np*MAX_DIM)
    call magnetic_density(m, st, rho, md)

    mm(1) = dmf_integrate(m, md(:, 1))
    mm(2) = dmf_integrate(m, md(:, 2))
    mm(3) = dmf_integrate(m, md(:, 3))

    deallocate(md)

    call pop_sub()
  end subroutine magnetic_moment


  ! ---------------------------------------------------------
  subroutine magnetic_local_moments(m, st, geo, rho, r, lmm)
    type(mesh_t),     intent(in)  :: m
    type(states_t),   intent(in)  :: st
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(in)  :: rho(:,:) ! (m%np_part, st%d%nspin)
    FLOAT,            intent(in)  :: r
    FLOAT,            intent(out) :: lmm(MAX_DIM, geo%natoms)

    integer :: ia, i
    FLOAT :: ri
    FLOAT, allocatable :: md(:, :), aux(:, :)

    call push_sub('magnetic.magnetic_local_moments')

    ALLOCATE(md (m%np, MAX_DIM), m%np*MAX_DIM)
    ALLOCATE(aux(m%np, MAX_DIM), m%np*MAX_DIM)

    call magnetic_density(m, st, rho, md)
    lmm = M_ZERO
    do ia = 1, geo%natoms
      aux = M_ZERO
      do i = 1, m%np
        call mesh_r(m, i, ri, a=geo%atom(ia)%x)
        if (ri > r) cycle
        aux(i, 1:MAX_DIM) = md(i, 1:MAX_DIM)
      end do
      lmm(1, ia) = dmf_integrate(m, aux(:, 1))
      lmm(2, ia) = dmf_integrate(m, aux(:, 2))
      lmm(3, ia) = dmf_integrate(m, aux(:, 3))
    end do
    deallocate(md, aux)

    call pop_sub()
  end subroutine magnetic_local_moments


  ! ---------------------------------------------------------
  subroutine calc_physical_current(gr, st, j)
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: st
    FLOAT,            intent(out)   :: j(:,:,:)   ! j(gr%mesh%np, gr%mesh%sb%dim, st%d%nspin)

    call push_sub('magnetic.calc_physical_current')

    ! Paramagnetic contribution to the physical current
    call states_calc_tau_jp_gn(gr, st, jp = j)

    ! TODO
    ! Diamagnetic contribution to the physical current

    call pop_sub()
  end subroutine calc_physical_current


  ! ---------------------------------------------------------
  ! This soubroutine receives as input a current, and produces
  ! as an output the vector potential that it induces.
  ! WARNING: There is probably a problem for 2D. For 1D none of this makes sense?
  subroutine magnetic_induced(gr, st, a_ind, b_ind)
    type(grid_t), intent(inout) :: gr
    type(states_t), intent(inout) :: st
    FLOAT, intent(out) :: a_ind(:, :) ! a(gr%mesh%np_part, gr%mesh%sb%dim)
    FLOAT, intent(out) :: b_ind(:, :) ! b(gr%mesh%np_part, gr%mesh%sb%dim) if gr%mesh%sb%dim=3, b(gr%mesh%np_part, 1) if gr%mesh%sb%dim=2

    integer :: i
    FLOAT, allocatable :: j(:, :, :)

    call push_sub('magnetic.magnetic_induced')

    ! If the states are real, we should never have reached here, but
    ! just in case we return zero.
    if(st%wfs_type .eq. M_REAL) then
      a_ind = M_ZERO
      b_ind = M_ZERO
      call pop_sub(); return
    end if

    ALLOCATE(j(gr%mesh%np_part, gr%mesh%sb%dim, st%d%nspin), gr%mesh%np_part*gr%mesh%sb%dim*st%d%nspin)
    call states_calc_tau_jp_gn(gr, st, jp=j)

    a_ind = M_ZERO
    do i = 1, gr%mesh%sb%dim
      call dpoisson_solve(gr, a_ind(:, i), j(:, i, 1))
    end do
    ! This minus sign is introduced here because the current that has been used
    ! before is the "number current density", and not the "charge current density",
    ! and therefore there is a minus sign missing (electrons are negative charges...)
    a_ind = - a_ind / P_C

    call dderivatives_curl(gr%der, a_ind, b_ind)

    deallocate(j)
    call pop_sub()
  end subroutine magnetic_induced


end module magnetic_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
