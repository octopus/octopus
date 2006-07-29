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
!! $Id$

#include "global.h"

! This module contains routines necessary to the split operator
! methods defined in td_exp
module td_exp_split_m
  use global_m
  use messages_m
  use fft_m
  use lib_basic_alg_m
  use cube_function_m
  use mesh_function_m
  use functions_m
  use simul_box_m
  use mesh_m
  use grid_m
  use states_m
  use hamiltonian_m
  use external_pot_m

  implicit none

contains

  ! ---------------------------------------------------------
  ! Calculates psi = exp{factor*T} psi
  ! where T is the kinetic energy operator
  subroutine zexp_kinetic (gr, h, psi, cf, factor)
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_t), intent(in)    :: h
    CMPLX,               intent(inout) :: psi(:,:) ! (NP_PART, dim)
    type(zcf_t),         intent(inout) :: cf
    CMPLX,               intent(in)    :: factor

    integer :: ix, iy, iz, k(MAX_DIM), idim
    FLOAT :: cutoff, temp(MAX_DIM), g2

    call push_sub('td_exp_split.exp_kinetic')

    if(simul_box_is_periodic(gr%sb)) then
      message(1) = 'Internal error in exp_kinetic'
      call write_fatal(1)
    end if

    if(h%cutoff > M_ZERO) then
      cutoff = h%cutoff
    else
      cutoff = CNST(1e10)
    end if

    temp = M_ZERO
    temp(1:NDIM) = (M_TWO*M_Pi)/(cf%n(1:NDIM)*gr%m%h(1:NDIM))

    call zcf_alloc_RS(cf)
    call zcf_alloc_FS(cf)

    do idim = 1, h%d%dim
      call zmf2cf(gr%m, psi(:, idim), cf)
      call zcf_RS2FS(cf)

      do iz = 1, cf%n(3)
        k(3) = pad_feq(iz, cf%n(3), .true.)
        do iy = 1, cf%n(2)
          k(2) = pad_feq(iy, cf%n(2), .true.)
          do ix = 1, cf%n(1)
            k(1) = pad_feq(ix, cf%n(1), .true.)

            g2 = min(cutoff, sum((temp(1:NDIM)*k(1:NDIM))**2))
            cf%FS(ix, iy, iz) = exp(factor*g2/M_TWO)*cf%FS(ix, iy, iz)
          end do
        end do
      end do

      call zcf_FS2RS(cf)
      call zcf2mf(gr%m, cf, psi(:, idim))
    end do

    call zcf_free_RS(cf)
    call zcf_free_FS(cf)

    call pop_sub()
  end subroutine zexp_kinetic


  ! ---------------------------------------------------------
  ! Calculates psi = exp{factor*V_KS(t)} psi
  ! where V_KS is the Kohn-Sham potential
  subroutine zexp_vlpsi(gr, h, psi, ik, t, factor)
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_t), intent(in)    :: h
    CMPLX,               intent(inout) :: psi(:,:) ! (NP_PART, NDIM)
    integer,             intent(in)    :: ik
    FLOAT,               intent(in)    :: t
    CMPLX,               intent(in)    :: factor

    integer :: k

    call push_sub('td_exp_split.vlpsi')

    ! WARNING: spinors not yet supported.
    select case(h%d%ispin)
    case(UNPOLARIZED)
      psi(1:NP, 1) = exp(factor*(h%ep%vpsl(1:NP)+h%vhxc(1:NP, 1)))*psi(1:NP, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        psi(1:NP, 1) = exp(factor*(h%ep%vpsl(1:NP)+h%vhxc(1:NP, 1)))*psi(1:NP, 1)
      else
        psi(1:NP, 1) = exp(factor*(h%ep%vpsl(1:NP)+h%vhxc(1:NP, 2)))*psi(1:NP, 1)
      end if
    case(SPINORS)
      message(1) = 'Internal error in exp_vlpsi'
      call write_fatal(1)
    end select

    if(h%ep%no_lasers > 0) then
      if(h%gauge /= 1) then  ! only length gauge is supported
        message(1) = "Only the length gauge is supported in exp_vlpsi"
        call write_fatal(1)
      end if

      do k = 1, h%d%dim
        psi(1:NP, k) = exp( factor*epot_laser_scalar_pot(gr%m%np, gr, h%ep, t) ) * psi(1:NP, k)
      end do
    end if

    call pop_sub()
  end subroutine zexp_vlpsi


  ! ---------------------------------------------------------
  ! calculates psi = exp{factor V_nlpp} psi
  ! where V_nlpp is the non-local part of the pseudpotential
  subroutine zexp_vnlpsi (m, h, psi, factor_, order_)
    type(mesh_t),        intent(in) :: m
    type(hamiltonian_t), intent(in) :: h
    CMPLX,            intent(inout) :: psi(m%np, h%d%dim)
    CMPLX,               intent(in) :: factor_
    logical,             intent(in) :: order_

!!$    integer :: idim, ikbc, jkbc, &
!!$         ivnl_start, ivnl_end, step, kbc_start, kbc_end, ivnl
!!$    CMPLX :: uvpsi, p2, ctemp
!!$    CMPLX, allocatable :: lpsi(:), lHpsi(:)
    logical :: order
    CMPLX   :: factor
    CMPLX, allocatable :: initzpsi(:, :)

    call push_sub('td_exp_split.vnlpsi')

    message(1) = 'Error: zexp_vnlpsi is currently broken.'
    call write_fatal(1)

    ALLOCATE(initzpsi(m%np, 1:h%d%dim), m%np*h%d%dim)
!   just to avoid compiler warnings due to unused variables
    factor   = factor_
    order    = order_
    initzpsi = psi

!!$
!!$    dimension_loop: do idim = 1, h%d%dim
!!$
!!$    if(order) then
!!$      step = 1;  ivnl_start = 1; ivnl_end = h%ep%nvnl
!!$    else
!!$      step = -1; ivnl_start = h%ep%nvnl; ivnl_end  = 1
!!$    end if
!!$
!!$    do ivnl = ivnl_start, ivnl_end, step
!!$!      /*
!!$       ALLOCATE( lpsi(h%ep%vnl(ivnl)%n), h%ep%vnl(ivnl)%n)
!!$       ALLOCATE(lhpsi(h%ep%vnl(ivnl)%n), h%ep%vnl(ivnl)%n)
!!$!      */
!!$       lpsi(:) = initzpsi(h%ep%vnl(ivnl)%jxyz(:), idim)
!!$       lhpsi(:) = M_z0
!!$
!!$       if(order) then
!!$          kbc_start = 1; kbc_end = h%ep%vnl(ivnl)%c
!!$       else
!!$          kbc_start = h%ep%vnl(ivnl)%c; kbc_end = 1
!!$       end if
!!$
!!$       do ikbc = kbc_start, kbc_end, step
!!$          do jkbc = kbc_start, kbc_end, step
!!$            stop 'does not work because of vol_pp'
!!$!             p2 = sum(h%ep%vnl(ivnl)%uv(:, ikbc)*h%ep%vnl(ivnl)%uv(:, ikbc))*m%vol_pp
!!$             ctemp = h%ep%vnl(ivnl)%uvu(ikbc, jkbc)*p2*factor
!!$!             uvpsi = sum(h%ep%vnl(ivnl)%uv(:, ikbc)*lpsi(:))*m%vol_pp*(exp(ctemp)-M_z1)/p2
!!$             lhpsi(:) = lhpsi(:) + uvpsi*h%ep%vnl(ivnl)%uv(:, jkbc)
!!$          end do
!!$       end do
!!$
!!$       deallocate(lpsi, lhpsi)
!!$    end do
!!$
!!$    end do dimension_loop
!!$
    deallocate(initzpsi)

    call pop_sub()
  end subroutine zexp_vnlpsi

end module td_exp_split_m
