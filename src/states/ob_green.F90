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

! Lead related routines for open boundaries.

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
  ! Calculate head of the semi-infinite surface green function with the
  ! algorithm from the paper
  ! Highly convergent schemes for the calculation of bulk and surface Green functions
  ! M. P. Lopez Sanco, J. M. Sancho, and J. Rubio (1984)
  ! J. Phys. F: Met. Phys. 15 (1985) 851-858
  subroutine lead_green(energy, diag, offdiag, np, green, dx)
    FLOAT,   intent(in)  :: energy        ! Energy to calculate Green function for.
    CMPLX,   intent(in)  :: diag(:, :)    ! Diagonal block of lead Hamiltonian.
    CMPLX,   intent(in)  :: offdiag(:, :) ! Off-diagonal block of lead Hamiltonian.
    integer, intent(in)  :: np            ! Number of interface points.
    CMPLX,   intent(out) :: green(:, :)   ! The calculated Green function.
    FLOAT,   intent(in)  :: dx            ! Spacing in transport direction.

    CMPLX, allocatable :: e(:, :), es(:, :), a(:, :), b(:, :), inv(:, :)
    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
    integer            :: i, j
    FLOAT              :: det, old_norm, norm, threshold, res

    call push_sub('ob_lead.lead_green')

    ALLOCATE(e(np, np), np**2)
    ALLOCATE(es(np, np), np**2)
    ALLOCATE(a(np, np), np**2)
    ALLOCATE(b(np, np), np**2)
    ALLOCATE(inv(np, np), np**2)
    ALLOCATE(tmp1(np, np), np**2)
    ALLOCATE(tmp2(np, np), np**2)
    ALLOCATE(tmp3(np, np), np**2)

    threshold = CNST(1e-12) ! FIXME: read from input.

    ! Fill with start values.
    call lalg_copy(np**2, diag(:, 1), e(:, 1))
    call lalg_copy(np**2, offdiag(:, 1), a(:, 1))
    do i = 1, np
      b(i, :) = offdiag(:, i)
    end do
    call lalg_copy(np**2, diag(:, 1), es(:, 1))
    old_norm = M_ZERO

    do i = 1, 1000 ! FIXME: read from input. 2^1000 efective layers
      ! inv <- -e
      inv = M_z0
      call lalg_axpy(np**2, -M_z1, e(:, 1), inv(:, 1))

      do j = 1, np
        inv(j, j) = inv(j, j) + energy + threshold*M_zI
      end do
      det = lalg_inverter(np, inv, invert=.true.)

      call lalg_gemm(np, np, np, M_z1, a, inv, M_z0, tmp2)
      call lalg_gemm(np, np, np, M_z1, tmp2, b, M_z0, tmp1)

      call lalg_axpy(np**2, M_z1, tmp1(:, 1), es(:, 1)) ! es <- es + tmp1

      norm = infinity_norm(es)*dx**2
      res  = abs(norm-old_norm)

      if(res.lt.threshold) then
        exit
      end if

      old_norm = norm

      call lalg_axpy(np**2, M_z1, tmp1(:, 1), e(:, 1)) ! e <- e + tmp1

      call lalg_gemm(np, np, np, M_z1, b, inv, M_Z0, tmp1)
      call lalg_gemm(np, np, np, M_z1, tmp1, a, M_z1, e)
      call lalg_copy(np**2, a(:, 1), tmp3(:, 1))
      call lalg_gemm(np, np, np, M_z1, tmp2, tmp3, M_z0, a)
      call lalg_copy(np**2, b(:, 1), tmp3(:, 1))
      call lalg_gemm(np, np, np, M_z1, tmp1, tmp3, m_z0, b)
    end do

    ! green <- -es
    green = M_z0
    call lalg_axpy(np**2, -M_z1, es(:, 1), green(:, 1))

    do j = 1, np
      green(j, j) = green(j, j) + energy + threshold*M_zI
    end do
    det = lalg_inverter(np, green, invert = .true.)
    call matrix_symmetric_average(green, np)
    det = aimag(green(1, 1))
    do i = 2, np
      det = det + aimag(green(i, i))
    end do
    if(det.gt.M_ZERO) then
      green = conjg(green)
      det   = -det
    end if

    deallocate(e, es, a, b, inv, tmp1, tmp2, tmp3)
    call pop_sub()
  end subroutine lead_green



end module ob_green_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
