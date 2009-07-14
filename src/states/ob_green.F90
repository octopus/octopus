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
!! $Id: h.F90 4037 2008-04-03 13:30:00Z xavier $

! Lead-related routines for open boundaries.

#include "global.h"

module ob_green_m
  use global_m
  use grid_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_parser_m
  use math_m
  use messages_m
  use nl_operator_m
  use profiling_m
  use ob_interface_m
  use simul_box_m
  use string_m

  implicit none
  private

  public ::         &
    lead_green

contains

  ! ---------------------------------------------------------
  ! Calculate head of the semi-infinite surface Green`s function with the
  ! algorithm from the paper
  ! Highly convergent schemes for the calculation of bulk and surface Green functions
  ! M. P. Lopez Sanco, J. M. Sancho, and J. Rubio (1984)
  ! J. Phys. F: Met. Phys. 15 (1985) 851-858
  subroutine lead_green(energy, diag, offdiag, np, green, h_is_real)
    FLOAT,   intent(in)  :: energy        ! Energy to calculate Green`s function for.
    CMPLX,   intent(in)  :: diag(:, :)    ! Diagonal block of lead Hamiltonian.
    CMPLX,   intent(in)  :: offdiag(:, :) ! Off-diagonal block of lead Hamiltonian.
    integer, intent(in)  :: np            ! Number of interface points.
    CMPLX,   intent(out) :: green(:, :)   ! The calculated Green`s function.
    logical, intent(in)  :: h_is_real     ! Is the Hamiltonian real? (no vector potential)

    CMPLX, allocatable :: e(:, :), es(:, :), a(:, :), b(:, :), inv(:, :)
    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
    CMPLX              :: cmplx_energy
    integer            :: i, j
    FLOAT              :: det, old_norm, norm, threshold, res, eta

    call push_sub('ob_lead.lead_green')

    SAFE_ALLOCATE(   e(1:np, 1:np))
    SAFE_ALLOCATE(  es(1:np, 1:np))
    SAFE_ALLOCATE(   a(1:np, 1:np))
    SAFE_ALLOCATE(   b(1:np, 1:np))
    SAFE_ALLOCATE( inv(1:np, 1:np))
    SAFE_ALLOCATE(tmp1(1:np, 1:np))
    SAFE_ALLOCATE(tmp2(1:np, 1:np))
    SAFE_ALLOCATE(tmp3(1:np, 1:np))

    tmp1 = M_ZERO
    tmp2 = M_ZERO
    tmp3 = M_ZERO

    eta = CNST(1e-7) ! FIXME: read from input.
    threshold = CNST(1e-14) ! FIXME: read from input.

    ! Fill with start values.
    e(1:np, 1:np) = diag(1:np, 1:np)
    a(1:np, 1:np) = offdiag(1:np, 1:np)
    forall (j = 1:np) b(j, :) = offdiag(:, j)
    es(1:np, 1:np) = diag(1:np, 1:np)
    old_norm = M_ZERO
    cmplx_energy = energy + eta*M_zI

    do i = 1, 1000 ! FIXME: read from input. 2^1000 efective layers
      inv(1:np, 1:np) = - e(1:np, 1:np)

      forall (j = 1:np) inv(j, j) = inv(j, j) + cmplx_energy
      det = lalg_inverter(np, inv, invert=.true.)

      call lalg_gemm(np, np, np, M_z1, a, inv, M_z0, tmp2)
      call lalg_gemm(np, np, np, M_z1, tmp2, b, M_z0, tmp1)

      es(1:np, 1:np) = es(1:np, 1:np) + tmp1(1:np, 1:np)

      norm = infinity_norm(es)
      res  = abs(old_norm/norm-M_ONE)

      if(res.lt.threshold) exit

      old_norm = norm

      e(1:np, 1:np) = e(1:np, 1:np) + tmp1(1:np, 1:np)
      call lalg_gemm(np, np, np, M_z1, b, inv, M_Z0, tmp1)
      call lalg_gemm(np, np, np, M_z1, tmp1, a, M_z1, e)
      tmp3(1:np, 1:np) = a(1:np, 1:np)
      call lalg_gemm(np, np, np, M_z1, tmp2, tmp3, M_z0, a)
      tmp3(1:np, 1:np) = b(1:np, 1:np)
      call lalg_gemm(np, np, np, M_z1, tmp1, tmp3, m_z0, b)
    end do

    green(1:np, 1:np) = - es(1:np, 1:np)

    forall (j = 1:np) green(j, j) = green(j, j) + cmplx_energy
    det = lalg_inverter(np, green, invert = .true.)
    if (h_is_real) then ! the Green`s function is complex symmetric
      call matrix_symmetric_average(green, np)
    end if ! otherwise it is general complex

    ! calculate DOS=-Tr(imag(green)) and check if negative
    det = aimag(green(1, 1))
    do j = 2, np
      det = det + aimag(green(j, j))
    end do
    if(det.gt.M_ZERO) then
      green = conjg(green)
    end if

    if(in_debug_mode) then ! write some info
      ! 1. calculate the residual InfNorm(Inverse(energy-h-offdiag*green*offdiag^T)-green)
      inv(1:np, 1:np) = -diag(1:np, 1:np)
      forall (j = 1:np) inv(j, j) = inv(j, j) + energy
      a(1:np, 1:np) = offdiag(1:np, 1:np)
      forall (j = 1:np) b(j, :) = offdiag(:, j)
      call lalg_gemm(np, np, np, M_z1, a, green, M_z0, tmp1)
      call lalg_gemm(np, np, np, M_z1, tmp1, b, M_z0, tmp2)
      inv(1:np, 1:np) = inv(1:np, 1:np) - tmp2(1:np, 1:np)
      det = lalg_inverter(np, inv, invert = .true.)
      inv(1:np, 1:np) = inv(1:np, 1:np)-green(1:np, 1:np)
      norm = infinity_norm(inv)
      write(*,*) 'Iterations:    ', i-1
      write(*,*) 'InfNorm(green):', infinity_norm(green)
      write(*,*) 'Green-Residual:', norm
    end if

    SAFE_DEALLOCATE_A(e)
    SAFE_DEALLOCATE_A(es)
    SAFE_DEALLOCATE_A(a)
    SAFE_DEALLOCATE_A(b)
    SAFE_DEALLOCATE_A(inv)
    SAFE_DEALLOCATE_A(tmp1)
    SAFE_DEALLOCATE_A(tmp2)
    SAFE_DEALLOCATE_A(tmp3)
    call pop_sub()
  end subroutine lead_green


end module ob_green_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
