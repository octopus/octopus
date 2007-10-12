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
!! $Id: td_transport.F90 3030 2007-06-25 16:45:05Z marques $

! Calculation of the memory coefficients for the modified Crank-Nicholson propgator.

#include "global.h"

module td_trans_mem_m
  use datasets_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use lib_oct_parser_m
  use messages_m
  use nl_operator_m
  use system_m
  use td_trans_intf_m
  use varinfo_m

  implicit none

  private
  public ::             &
    mbytes_memory_term, &
    memory_init,        &
    memory_end

  integer :: mem_iter
  FLOAT   :: mem_tolerance

contains

  ! ---------------------------------------------------------
  ! Calculate amount of (machine) memory required for
  ! the memory term (in MBytes).
  ! Does not consider spin at the moment.
  integer function mbytes_memory_term(iter, np, nl, sys)
    integer,        intent(in) :: iter
    integer,        intent(in) :: np
    integer,        intent(in) :: nl
    type(system_t), intent(in) :: sys

    call push_sub('td_transport.mbytes_memory_term')

    mbytes_memory_term = nl*iter*2*REAL_PRECISION* &
      np*sys%st%nst / 2**20

    call pop_sub()
  end function mbytes_memory_term


  ! ---------------------------------------------------------
  ! Allocate (machine) memory for the memory coefficents and
  ! calculate them.
  subroutine memory_init(intface, delta, max_iter, op, coeff)
    type(intface_t),     intent(in) :: intface
    FLOAT,               intent(in) :: delta
    integer,             intent(in) :: max_iter
    type(nl_operator_t), intent(in) :: op
    CMPLX,               pointer    :: coeff(:, :, :, :)

    integer            :: il
    FLOAT, allocatable :: h_re(:, :), h_im(:, :)

    call push_sub('td_trans_mem.memory_coeff')

    ALLOCATE(coeff(intface%np, intface%np, 0:max_iter, NLEADS), intface%np**2*(max_iter+1)*NLEADS)
    ALLOCATE(h_re(intface%np, intface%np), intface%np**2)
    ALLOCATE(h_im(intface%np, intface%np), intface%np**2)

    ! Get iteration parameters from input file.

    !%Variable TDTransMemTol
    !%Type float
    !%Default 1e-8
    !%Section Transport
    !%Description
    !% Decides when to consider the memory coefficients converged.
    !%End
    call loct_parse_float(check_inp('TDTransMemTol'), CNST(1e-8), mem_tolerance)
    if(mem_tolerance.le.M_ZERO) then
      write(message(1), '(a,f14.6,a)') "Input : '", mem_tolerance, "' is not a valid TDTransMemTol."
      message(2) = '(0 < TDTransMemTol)'
      call write_fatal(2)
    end if

    !%Variable TDTransMemMaxIter
    !%Type integer
    !%Default 25
    !%Section Transport
    !%Description
    !% Sets the maximum iteration number to converge the memory coefficients.
    !%End
    call loct_parse_int(check_inp('TDTransMemMaxIter'), 25, mem_iter)
    if(mem_iter.le.0) then
      write(message(1), '(a,i6,a)') "Input : '", mem_iter, "' is not a valid TDTransMemMaxIter."
      message(2) = '(0 < TDTransMemMaxIter)'
      call write_fatal(2)
    end if

    do il = 1, NLEADS
      ! Get diagonal matrix.
      call calculate_diag(op, intface, il, h_re, h_im)

      ! Get anchor for recursion.
      call approx_coeff0(intface, delta, il, op, h_re, coeff(:, :, 0, il))
      ! Calculate the subsequent coefficients by the recursive relation.
      call calculate_coeffs(il, max_iter, delta, op, intface, h_re, coeff(:, :, :, il))
    end do

    message(1) = 'Info: Coefficients for memory term calculated.'
    call write_info(1)
    
    call write_coeffs()

    deallocate(h_re, h_im)

    call pop_sub()
  end subroutine memory_init


  ! ---------------------------------------------------------
  ! Solve for zeroth memory coefficient by truncating the continued
  ! matrix fraction.
  subroutine approx_coeff0(intface, delta, il, op, h_re, coeff0)
    type(intface_t),     intent(in)  :: intface
    FLOAT,               intent(in)  :: delta
    integer,             intent(in)  :: il
    type(nl_operator_t), intent(in)  :: op
    FLOAT,               intent(in)  :: h_re(:, :)
    CMPLX,               intent(out) :: coeff0(:, :)

    integer            :: i, j
    CMPLX, allocatable :: coeff0_old(:, :), tmp(:, :)

    call push_sub('td_trans_mem.approx_coeff0')

    ALLOCATE(coeff0_old(intface%np, intface%np), intface%np**2)
    ALLOCATE(tmp(intface%np, intface%np), intface%np**2)

    ! Truncating the continued fraction is the same as iterating the
    ! equation
    ! 
    !     coeff0 = V(op) * (1/(1 + i*delta*h(op) + delta^2*coeff0)) * V(op)^T
    ! 
    ! with start value 0:
    coeff0(:, :) = 0

    do i = 1, mem_iter
      coeff0_old = coeff0

      ! Calculate 1 + i*delta*h(op) + delta^2*coeff0_old
      coeff0 = M_zI*delta*h_re + delta**2*coeff0
      do j = 1, intface%np
        coeff0(j, j) = 1 + coeff0(j, j) 
      end do

      ! Invert.
      call lalg_svd_inverse(intface%np, coeff0)

      ! Apply coupling matrices.
      call apply_coupling(coeff0, op, LEFT, intface, il, tmp)
      call apply_coupling(tmp, op, RIGHT, intface, il, coeff0)

      if(abs(infinity_norm(coeff0)-infinity_norm(coeff0_old)).lt.mem_tolerance) then
        exit
      end if
    end do

    if(i.gt.mem_iter) then
      write(message(1), '(a,i6,a)') 'Memory coefficent for time step 0, ' &
        //trim(lead_name(il))//' lead, not converged'
      call write_warning(1)
    end if
    call pop_sub()
  end subroutine approx_coeff0


  ! ---------------------------------------------------------
  ! Caclulate infinity-norm of matrix.
  FLOAT function infinity_norm(matrix)
    CMPLX, intent(in) :: matrix(:, :)

    integer :: m_min, m_max, i
    FLOAT   :: norm_old, norm

    call push_sub('td_trans_mem.inf_norm')

    norm = 0

    m_min = lbound(matrix, 1)
    m_max = ubound(matrix, 1)
    do i = m_min, m_max
      norm_old = norm
      norm = sum(abs(matrix(i, :)))
      norm = max(norm, norm_old)
    end do
    
    infinity_norm = norm

    call pop_sub()
  end function infinity_norm


  ! ---------------------------------------------------------
  ! Calculate the diagonal block matrix h(op), i. e. the Hamiltonian for
  ! entries contained in the interface region (for zero potential only
  ! at the moment, thus, only kinetic energy).
  subroutine calculate_diag(op, intf, il, h_re, h_im)
    type(nl_operator_t), intent(in)  :: op
    type(intface_t),     intent(in)  :: intf
    integer,             intent(in)  :: il
    FLOAT,               intent(out) :: h_re(:, :), h_im(:, :)

    integer :: i, k, n, k_stencil
    FLOAT   :: w_re, w_im

    call push_sub('td_trans_mem.calculate_diag')

    ! For all interface points...
    do i = 1, intf%np
      n = intf%index(i, il)
      ! ... find the points they are coupled to.
      do k = 1, op%n
        k_stencil = op%i(k, n)
        ! If the coupling point is in the interface...
        if(k_stencil.le.op%np.and. &
          member_of_intface(k_stencil, intf, il)) then
          ! ... get the operator coefficients and...
          if(op%cmplx_op) then
            if(op%const_w) then
              w_re = op%w_re(k, 1)
              w_im = op%w_im(k, 1)
            else
              w_re = op%w_re(k, n)
              w_im = op%w_im(k, n)
            end if
          else
            if(op%const_w) then
              w_re = op%w_re(k, 1)
            else
              w_re = op%w_re(k, n)
            end if
            w_im = w_re
          end if

          ! Calculation if the kinetic term: -1/2 prefactor.
          w_im = -w_im/2
          w_re = -w_re/2

          ! ... write them into the right entry of the diagonal block.
          h_re(i, k_stencil-intf%index_range(1, il)+1) = w_re
          h_im(i, k_stencil-intf%index_range(1, il)+1) = w_im
        end if
      end do
    end do

    call pop_sub()
  end subroutine calculate_diag


  ! ---------------------------------------------------------
  ! Multiplies the coupling block matrix V(op) for interface il
  ! onto matrix.
  ! If il = RIGHT, then the coupling is done to the right,
  ! for il = LEFT, it is to the left.
  ! If lr = RIGHT, then matrix*V(op) is calculated, for lr = LEFT,
  ! V(op)*matrix is calculated
  subroutine apply_coupling(matrix, op, lr, intface, il, v)
    CMPLX,               intent(in)  :: matrix(:, :)
    type(nl_operator_t), intent(in)  :: op
    integer,             intent(in)  :: lr
    type(intface_t),     intent(in)  :: intface
    integer,             intent(in)  :: il
    CMPLX,               intent(out) :: v(:, :)

    integer :: p_n(MAX_DIM), p_k(MAX_DIM), p_matr(MAX_DIM)
    integer :: i, j, k, k_stencil
    integer :: n, n_matr
    integer :: dir
    integer :: x_shift
    FLOAT   :: w_re, w_im

    call push_sub('td_trans_mem.apply_coupling')

    ! Coupling direction.
    select case(il)
    case(LEFT)
      dir = 1
    case(RIGHT)
      dir = -1
    end select

    do i = 1, intface%np
      do j = 1, intface%np
        n = intface%index(i, il)
        if(lr.eq.RIGHT) then
          v(i, j) = 0
        else
          v(j, i) = 0
        end if
        ! Iterate over all stencil points.
        do k = 1, op%n
          k_stencil = op%i(k, n)
          if(k_stencil.le.op%np.and. & ! Is member of simulation region.
            .not.member_of_intface(k_stencil, intface, il)) then
            ! Get coordinates of current point in interface and coupling point.
            p_n = op%m%lxyz(n, :)
            p_k = op%m%lxyz(k_stencil, :)
            ! Get the distance in transport direction.
            x_shift = dir*(p_k(TRANS_DIR)-p_n(TRANS_DIR))
            ! Calculate coordinates of the corresponding point of the matrix.
            p_matr            = p_n
            p_matr(TRANS_DIR) = p_matr(TRANS_DIR) - dir*(intface%extent-x_shift)
            ! Finally, get its point number in the eigenstate array.
            n_matr = intface_index(op%m%lxyz_inv(p_matr(1), p_matr(2), p_matr(3)), intface, il)
            ! Multiply by the coefficient of the operator and sum up.
            if(op%cmplx_op) then
              if(op%const_w) then
                w_re = op%w_re(k, 1)
                w_im = op%w_im(k, 1)
              else
                w_re = op%w_re(k, n)
                w_im = op%w_im(k, n)
              end if
            else
              if(op%const_w) then
                w_re = op%w_re(k, 1)
              else
                w_re = op%w_re(k, n)
              end if
              w_im = w_re
            end if

            ! Calculation of the kinetic term: -1/2 prefactor.
            w_im = -w_im/2
            w_re = -w_re/2

            ! Sum up.
            if(lr.eq.RIGHT) then
              v(i, j) = v(i, j) + &
                cmplx(w_re*real(matrix(n_matr, j)), w_im*aimag(matrix(n_matr, j)))
            else
              v(j, i) = v(j, i) + &
                cmplx(w_re*real(matrix(n_matr, j)), w_im*aimag(matrix(n_matr, j)))
            end if
          end if
        end do
      end do
    end do

    call pop_sub()
  end subroutine apply_coupling


  ! ---------------------------------------------------------
  ! coeffs(:, :, :, 0) given, calculate the subsequent ones by
  ! the recursive relation.
  subroutine calculate_coeffs(il, iter, delta, op, intf, h_re, coeffs)
    integer,             intent(in)    :: il
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: delta
    type(nl_operator_t), intent(in)    :: op
    type(intface_t),     intent(in)    :: intf
    FLOAT,               intent(in)    :: h_re(:, :)
    CMPLX,               intent(inout) :: coeffs(intf%np, intf%np, 0:iter)

    integer            :: i,j, k, np
    CMPLX, allocatable :: coeff_p(:, :, :), p_prev(:, :), tmp(:, :)
    CMPLX, allocatable :: prefactor_plus(:, :), prefactor_minus(:, :)

    call push_sub('td_trans_mem.calculate_coeffs')

    ALLOCATE(coeff_p(intf%np, intf%np, 0:iter), intf%np**2*(iter+1))
    ALLOCATE(p_prev(intf%np, intf%np), intf%np**2)
    ALLOCATE(tmp(intf%np, intf%np), intf%np**2)
    ALLOCATE(prefactor_plus(intf%np, intf%np), intf%np**2)
    ALLOCATE(prefactor_minus(intf%np, intf%np), intf%np**2)

    np = intf%np

    coeff_p = 0

    prefactor_plus  = M_zI*delta*h_re
    prefactor_minus = M_zI*delta*h_re
    do i = 1, intf%np
      prefactor_plus(i, i)  = 1 + prefactor_plus(i, i)
      prefactor_minus(i, i) = 1 - prefactor_minus(i, i)
    end do

    ! Calculate p_\alpha = 1/(1 + i*delta*h_\alpha + delta^2*q_\alpha)
    coeff_p(:, :, 0) = prefactor_plus + delta**2 * coeffs(:, :, 0)
    call lalg_svd_inverse(intf%np, coeff_p(:, :, 0))
    call lalg_svd_inverse(intf%np, prefactor_plus)

    do i = 1, iter
      do j = 1, mem_iter
        p_prev = coeff_p(:, :, i)

        ! (3)
        call apply_coupling(p_prev, op, LEFT, intf, il, tmp)
        call apply_coupling(tmp, op, RIGHT, intf, il, coeff_p(:, :, i))
        coeff_p(:, :, i) = coeff_p(:, :, i) + 2*coeffs(:, :, i-1)
        if(i.gt.1) then
          coeff_p(:, :, i) = coeff_p(:, :, i) + coeffs(:, :, i-2)
        end if

        call lalg_gemm(np, np, np, M_z1, coeff_p(:, :, i), coeff_p(:, :, 0), &
          M_z0, tmp)
        coeff_p(:, :, i) = tmp

        ! (1)
        call lalg_gemm(np, np, np, M_z1, coeffs(:, :, 0), p_prev, M_z1, coeff_p(:, :, i))

        ! (2)
        do k = 1, i-1
          tmp = coeffs(:, :, k) + 2*coeffs(:, :, k-1)
          if(k.gt.1) then 
            tmp = tmp + coeffs(:, :, k-2)
          end if
          call lalg_gemm(np, np, np, M_z1, tmp, coeff_p(:, :, i-k), M_z1, coeff_p(:, :, i))
        end do

        coeff_p(:, :, i) = -delta**2 * coeff_p(:, :, i)
        call lalg_gemm(np, np, np, M_z1, prefactor_minus, coeff_p(:, :, i-1), &
          M_z1, coeff_p(:, :, i))

        call lalg_gemm(np, np, np, M_z1, prefactor_plus, coeff_p(:, :, i), &
          M_z0, tmp)
        coeff_p(:, :, i) = tmp

        if(abs(infinity_norm(coeff_p(:, :, i))-infinity_norm(p_prev)).lt.mem_tolerance) then
          exit
        end if
      end do

      ! Write a warning if a coefficient is not converged.
      if(j.gt.mem_iter) then
        write(message(1), '(a,i6,a)') 'Memory coefficent for time step ', i, &
          ', '//trim(lead_name(il))//' lead, not converged'
        call write_warning(1)
      end if

      call apply_coupling(coeff_p(:, :, i), op, LEFT, intf, il, tmp)
      call apply_coupling(tmp, op, RIGHT, intf, il, coeffs(:, :, i))
    end do
    
    deallocate(coeff_p, p_prev, tmp)
    deallocate(prefactor_plus, prefactor_minus)

    call pop_sub()
  end subroutine calculate_coeffs


  ! ---------------------------------------------------------
  ! Write memory coefficients to file.
  subroutine write_coeffs()
    call push_sub('td_trans_mem.write_coeffs')

    call pop_sub()
  end subroutine write_coeffs


  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine memory_end(coeff)
    CMPLX, pointer :: coeff(:, :, :, :)

    call push_sub('td_trans_mem.memory_end')

    if(associated(coeff)) then
      deallocate(coeff)
      nullify(coeff)
    end if

    call pop_sub()
  end subroutine memory_end
end module td_trans_mem_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
