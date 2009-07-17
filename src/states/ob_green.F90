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

  subroutine lead_green(energy, diag, offdiag, np, green, h_is_real, il, offdiag_invertible)
    FLOAT,   intent(in)  :: energy        ! Energy to calculate Green`s function for.
    CMPLX,   intent(in)  :: diag(:, :)    ! Diagonal block of lead Hamiltonian.
    CMPLX,   intent(in)  :: offdiag(:, :) ! Off-diagonal block of lead Hamiltonian.
    integer, intent(in)  :: np            ! Number of interface points.
    CMPLX,   intent(out) :: green(:, :)   ! The calculated Green`s function.
    logical, intent(in)  :: h_is_real   ! Is the hamiltonian real? (no vector potential)
    integer, intent(in)  :: il            ! which lead
    logical, intent(in)  :: offdiag_invertible ! Is the offdiagonal invertible

    FLOAT              :: threshold, eta, residual, residual2, eps
    CMPLX              :: c_energy
    CMPLX, allocatable :: green2(:, :)

    call push_sub('ob_lead.lead_green_new')

    eta = CNST(1e-7) ! FIXME: read from input.
    threshold = CNST(1e-14) ! FIXME: read from input.
    eps = CNST(1e-5) ! FIXME: read from input.
    c_energy = energy + eta*M_zI

    ! if the offdiagonal matrix is not invertible we have only one way of calculating
    ! the Green`s function: with the Sancho method
    if(.not.offdiag_invertible) then
      call lead_green_sancho(c_energy, diag, offdiag, np, green, threshold, h_is_real)
      residual = calc_residual_green(energy, green, diag, offdiag, np, il)
      if(in_debug_mode) then ! write info
        write(message(1), '(a,e10.3)') 'Sancho-Residual = ', residual
        call write_info(1)
      end if
    else
      ! 1. calculate with the fastest method
      call lead_green_umerski(c_energy, diag, offdiag, np, green, h_is_real, il)
      ! 2. check if converged
      residual = calc_residual_green(energy, green, diag, offdiag, np, il)
      if(in_debug_mode) then ! write info
        write(message(1), '(a,e10.3)') 'Umerski-Residual = ', residual
        call write_info(1)
      end if
      if(residual.gt.eps) then
        ! 3. if not converged try the other method while saving the old
        SAFE_ALLOCATE(green2(1:np, 1:np))
        call lead_green_sancho(c_energy, diag, offdiag, np, green2, threshold, h_is_real)
       ! call lead_green_umerski(c_energy, diag, offdiag, np, green2, h_is_real, il)
        residual2 = calc_residual_green(energy, green2, diag, offdiag, np, il)
        if(in_debug_mode) then ! write info
          write(message(1), '(a,e10.3)') 'Sancho-Residual = ', residual2
          call write_info(1)
        end if
        ! 4. if also not converged choose the best one and give warning
        if(residual2.lt.eps) then ! finally converged
          green = green2
        else ! write warning
          if(residual2.lt.residual) then
            green = green2
          end if
          message(1) = 'The surface green function did not converge properly'
          message(2) = 'with neither the decimation technique nor the closed form.'
          write(message(3), '(a,e10.3)') 'Umerski-Residual = ', residual
          write(message(4), '(a,e10.3)') 'Sancho-Residual  = ', residual2
          message(5) = 'Now the better converged version is taken.'
          call write_warning(5)
        end if
        ! cleanup
        SAFE_DEALLOCATE_A(green2)
      end if
    end if

    call pop_sub()
  end subroutine lead_green


  ! check if the green function gives the correct density of states
  ! if not compute its hermitian conjugate
  subroutine fix_green(np, green, dos)
    integer, intent(in)    :: np
    CMPLX,   intent(inout) :: green(:, :)
    FLOAT,   intent(out)   :: dos

    integer  :: j

    call push_sub('ob_lead.fix_green')

    ! calculate DOS=-Tr(imag(green)) and check if negative
    dos = -aimag(green(1, 1))
    do j = 2, np
      dos = dos - aimag(green(j, j))
    end do
    if(dos.lt.M_ZERO) then
      green(:, :) = transpose(conjg(green(:, :)))
    end if

    call pop_sub()
  end subroutine fix_green


  ! calculate the residual of the green function
  FLOAT function calc_residual_green(energy, green, diag, offdiag, np, il) result(residual)
    FLOAT,   intent(in) :: energy
    CMPLX,   intent(in) :: green(:, :)
    CMPLX,   intent(in) :: diag(:, :)
    CMPLX,   intent(in) :: offdiag(:, :)
    integer, intent(in) :: np
    integer, intent(in) :: il

    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
    FLOAT              :: det
    integer            :: j

    call push_sub('ob_lead.fix_green')

    SAFE_ALLOCATE(tmp1(1:np, 1:np))
    SAFE_ALLOCATE(tmp2(1:np, 1:np))

    ! 1. calculate the residual InfNorm(Inverse(energy-h-offdiag*green*offdiag^T)-green)
    tmp1(1:np, 1:np) = -diag(1:np, 1:np)
    forall (j = 1:np) tmp1(j, j) = tmp1(j, j) + energy

    tmp2(1:np, 1:np) = green(1:np, 1:np)
    if(il.eq.LEFT) then
      call lalg_trmm(np, np, 'U', 'N', 'L', M_z1, offdiag, tmp2)
      call lalg_trmm(np, np, 'U', 'T', 'R', M_z1, offdiag, tmp2)
    else
      call lalg_trmm(np, np, 'L', 'N', 'L', M_z1, offdiag, tmp2)
      call lalg_trmm(np, np, 'L', 'T', 'R', M_z1, offdiag, tmp2)
    end if

    tmp1(1:np, 1:np) = tmp1(1:np, 1:np) - tmp2(1:np, 1:np)
    det = lalg_inverter(np, tmp1, invert = .true.)
    tmp1(1:np, 1:np) = tmp1(1:np, 1:np) - green(1:np, 1:np)
    residual = infinity_norm(tmp1)

    SAFE_DEALLOCATE_A(tmp1)
    SAFE_DEALLOCATE_A(tmp2)

    call pop_sub()
  end function calc_residual_green


  ! ---------------------------------------------------------
  ! Calculate head of the semi-infinite surface Green`s function with the
  ! algorithm from the paper
  ! Highly convergent schemes for the calculation of bulk and surface Green`s functions
  ! M. P. Lopez Sanco, J. M. Sancho, and J. Rubio (1984)
  ! J. Phys. F: Met. Phys. 15 (1985) 851-858
  subroutine lead_green_sancho(energy, diag, offdiag, np, green, threshold, h_is_real)
    CMPLX,   intent(in)  :: energy        ! Energy to calculate Green`s function for (already shifted).
    CMPLX,   intent(in)  :: diag(:, :)    ! Diagonal block of lead Hamiltonian.
    CMPLX,   intent(in)  :: offdiag(:, :) ! Off-diagonal block of lead Hamiltonian.
    integer, intent(in)  :: np            ! Number of interface points.
    CMPLX,   intent(out) :: green(:, :)   ! The calculated Green`s function.
    FLOAT,   intent(in)  :: threshold     ! this defines convergence
    logical, intent(in)  :: h_is_real     ! Is the hamiltonian real? (no vector potential)

    CMPLX, allocatable :: e(:, :), es(:, :), a(:, :), b(:, :), inv(:, :)
    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
    integer            :: i, j
    FLOAT              :: dos, old_norm, norm, res

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

    ! Fill with start values.
    e(1:np, 1:np) = diag(1:np, 1:np)
    a(1:np, 1:np) = offdiag(1:np, 1:np)
    b(1:np, 1:np) = transpose(conjg(offdiag(1:np, 1:np)))
    
    es(1:np, 1:np) = diag(1:np, 1:np)
    old_norm = M_ZERO

    do i = 1, 1000 ! FIXME: read from input. 2^1000 efective layers
      inv(1:np, 1:np) = - e(1:np, 1:np)

      forall (j = 1:np) inv(j, j) = inv(j, j) + energy
      dos = lalg_inverter(np, inv, invert=.true.)

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

    forall (j = 1:np) green(j, j) = green(j, j) + energy
    dos = lalg_inverter(np, green, invert = .true.)
    if (h_is_real) then ! the Green`s function is complex symmetric
      call matrix_symmetric_average(green, np)
    end if ! otherwise it is general complex

    call fix_green(np, green, dos)

    if(in_debug_mode) then ! write some info
      write(message(1), '(a,e10.3)') 'Sancho-DOS = ', dos
      call write_info(1)
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
  end subroutine lead_green_sancho


  ! compute the semi-infinite surface green function
  ! Algorithm taken from A. Umerski, Closed-form solutions to surface Green`s functions
  ! http://www.city.ac.uk/sems/dps/mathematics/research/nanostructures/prb55_5266.pdf
  subroutine lead_green_umerski(energy, diag, offdiag, np, green, h_is_real, il)
    CMPLX,     intent(in)  :: energy
    CMPLX,     intent(in)  :: diag(:, :)
    CMPLX,     intent(in)  :: offdiag(:, :)
    integer,   intent(in)  :: np ! number of interface points
    CMPLX,     intent(out) :: green(:, :)
    logical,   intent(in)  :: h_is_real ! Is the hamiltonian real? (no vector potential)
    integer,   intent(in)  :: il

    integer              :: i, np2
    CMPLX, allocatable   :: x(:,:), s(:,:), o2(:,:), o4(:,:), d(:), unsorted_d(:)
    FLOAT, allocatable   :: abs_d(:)
    integer, allocatable :: index(:)
    FLOAT                :: dos

    call push_sub('ob_green.lead_green_umerski')

    np2 = 2*np
    SAFE_ALLOCATE( x(1:np2, 1:np2) )
    SAFE_ALLOCATE( s(1:np2, 1:np2) )
    SAFE_ALLOCATE( o2(1:np, 1:np) )
    SAFE_ALLOCATE( o4(1:np, 1:np) )
    SAFE_ALLOCATE( d(1:np2) )
    SAFE_ALLOCATE( unsorted_d(1:np2) )
    SAFE_ALLOCATE( abs_d(1:np2) )
    SAFE_ALLOCATE( index(1:np2) )

    ! 1. create matrix x ( x = {{0,offdiag^(-1)},{-offdiag^H,(energy-diag)*offdiag^(-1)}} )
    ! 2. diagonalize, calculate s, x = s*d*s^(-1)
    ! 3. extract submatrices ( o2 and 04 with s = {{o1,o2},{o3,o4}})
    ! 4. calculate green = o2*o4^(-1)

    ! 1. create matrix x ( x = {{0,offdiag^(-1)},{-offdiag^H,(energy-diag)*offdiag^(-1)}} )
    x(1:np, 1:np) = M_z0
    ! use o2 as tmp variable
    o2(1:np, 1:np) = offdiag(1:np, 1:np)
    if (il.eq.LEFT) then
      call lalg_invert_upper_triangular(np, o2)
    else
      call lalg_invert_lower_triangular(np, o2)
    end if
    x(1:np, np+1:np2) = o2(1:np, 1:np)
    x(np+1:np2, 1:np) = -transpose(conjg(offdiag(1:np, 1:np)))
    ! use o4 as tmp variable
    o4(1:np, 1:np) = -diag(1:np, 1:np)
    forall (i = 1:np) o4(i, i) = o4(i, i) + energy
    if (il.eq.LEFT) then
      call lalg_trmm(np, np, 'U', 'N', 'R', M_z1, o2, o4)
    else
      call lalg_trmm(np, np, 'L', 'N', 'R', M_z1, o2, o4)
    end if
    x(np+1:np2, np+1:np2) = o4(1:np, 1:np)

    ! 2. compute diagonalization matrix s, s^(-1)*x*s = d
    call lalg_eigensolve_nonh(np2, x, unsorted_d)
    ! the eigenvalues and the corresponding eigenvectors have to be sorted in the following way:
    ! let d = diag(d(1),d(2),d(3),...,d(2n))
    ! with |d(1)| < |d(2)| < ... < |d(2n)|
    ! sort the eigenvalues and eigenvectors
    abs_d(:) = abs(unsorted_d(:))
    call sort(abs_d, index)
    do i=1, np2
      d(i)    = unsorted_d(index(i))
      s(:, i) = x(:, index(i))
    end do
    ! 3. extract submatrices ( o2 and o4 of S = {{o1,o2},{o3,o4}};)
    o2(1:np, 1:np) = s(1:np, np+1:np2)
    o4(1:np, 1:np) = s(np+1:np2, np+1:np2)
    dos = lalg_inverter(np, o4, invert = .true.)
    call zgemm('N', 'N', np, np, np, M_z1, o2, np, o4, np, M_z0, green, np)

    call fix_green(np, green, dos)
    if(in_debug_mode) then ! write some info
      write(*,*) 'Umerski: DOS', dos
    end if

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(s)
    SAFE_DEALLOCATE_A(o2)
    SAFE_DEALLOCATE_A(o4)
    SAFE_DEALLOCATE_A(d)
    SAFE_DEALLOCATE_A(unsorted_d)
    SAFE_DEALLOCATE_A(abs_d)
    SAFE_DEALLOCATE_A(index)

    call pop_sub()
  end subroutine lead_green_umerski

end module ob_green_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
