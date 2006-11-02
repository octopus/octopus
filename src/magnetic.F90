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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: states.F90 2515 2006-10-24 17:13:30Z acastro $

#include "global.h"

module magnetic_m
  use global_m
  use messages_m
  use mesh_function_m
  use functions_m
  use mesh_m
  use states_m
  use geometry_m
  use grid_m

  implicit none

  private
  public ::                    &
    magnetic_density,          &
    magnetic_moment,           &
    magnetic_local_moments,    &
    calc_physical_current,     &
    calc_paramagnetic_current

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
  ! This routine (obviously) assumes complex wave-functions
  subroutine calc_paramagnetic_current(gr, st, jp)
    type(grid_t),   intent(inout) :: gr
    type(states_t), intent(inout) :: st
    FLOAT,          intent(out)   :: jp(:,:,:)  ! (NP, NDIM, st%d%nspin)

    integer :: ik, p, sp, k
    CMPLX, allocatable :: grad(:,:)
#ifdef  HAVE_MPI
    FLOAT, allocatable :: red(:,:,:)
#endif

    call push_sub('magnetic.calc_paramagnetic_current')

    ASSERT(st%d%wfs_type == M_CMPLX)

    if(st%d%ispin == SPIN_POLARIZED) then
      sp = 2
    else
      sp = 1
    end if

    jp = M_ZERO
    ALLOCATE(grad(NP_PART, NDIM), NP_PART*NDIM)

    do ik = 1, st%d%nik, sp
      do p  = st%st_start, st%st_end
        call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, 1, p, ik), grad)

        ! spin-up density
        do k = 1, NDIM
          jp(1:NP, k, 1) = jp(1:NP, k, 1) + st%d%kweights(ik)*st%occ(p, ik)       &
            * aimag(conjg(st%zpsi(1:NP, 1, p, ik)) * grad(1:NP, k))
        end do

        ! spin-down density
        if(st%d%ispin == SPIN_POLARIZED) then
          call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, 1, p, ik+1), grad)

          do k = 1, NDIM
            jp(1:NP, k, 2) = jp(1:NP, k, 2) + st%d%kweights(ik+1)*st%occ(p, ik+1) &
              * aimag(conjg(st%zpsi(1:NP, 1, p, ik+1)) * grad(1:NP, k))
          end do

          ! WARNING: the next lines DO NOT work properly
        else if(st%d%ispin == SPINORS) then ! off-diagonal densities
          call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, 2, p, ik), grad)

          do k = 1, NDIM
            jp(1:NP, k, 2) = jp(1:NP, k, 2) + st%d%kweights(ik)*st%occ(p, ik)     &
              * aimag(conjg(st%zpsi(1:NP, 2, p, ik)) * grad(1:NP, k))
          end do
        end if

      end do
    end do
    deallocate(grad)

#if defined(HAVE_MPI)
    ALLOCATE(red(NP_PART, NDIM, st%d%nspin), NP_PART*NDIM*st%d%nspin)
    call MPI_Allreduce(jp(1, 1, 1), red(1, 1, 1), NP*NDIM*st%d%nspin,       &
      MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    jp = red
    deallocate(red)
#endif

    call pop_sub()
  end subroutine calc_paramagnetic_current


  ! ---------------------------------------------------------
  subroutine calc_physical_current(gr, st, j)
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: st
    FLOAT,            intent(out)   :: j(:,:,:)   ! j(NP, NDIM, st%d%nspin)

    call push_sub('magnetic.calc_physical_current')

    ! Paramagnetic contribution to the physical current
    call calc_paramagnetic_current(gr, st, j)

    ! TODO
    ! Diamagnetic contribution to the physical current

    call pop_sub()
  end subroutine calc_physical_current


end module magnetic_m
