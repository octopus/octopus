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

! Calculation of the memory coefficients for the modified
! Crank-Nicholson propagator.

#include "global.h"

module ob_mem_m
  use datasets_m
  use global_m
  use hamiltonian_m
  use io_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use loct_parser_m
  use math_m
  use messages_m
  use mpi_m
  use nl_operator_m
  use ob_interface_m
  use ob_lead_m
  use ob_terms_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_dim_m
  use system_m
  use varinfo_m

  implicit none

  integer, parameter, public :: &
    SAVE_CPU_TIME  = 1,         &
    SAVE_RAM_USAGE = 2

  private
  public ::             &
    mbytes_memory_term, &
    ob_mem_init,        &
    ob_mem_end,         &
    make_full_matrix

  integer :: mem_iter
  FLOAT   :: mem_tolerance

contains

  ! ---------------------------------------------------------
  ! Allocate (machine) memory for the memory coefficents and
  ! calculate them.
  subroutine ob_mem_init(intface, hm, ob, delta, max_iter, op, &
    spacing, order, mpi_grp)
    type(interface_t),   intent(in)    :: intface(NLEADS)
    type(hamiltonian_t), intent(in)    :: hm
    type(ob_terms_t),    intent(inout) :: ob
    FLOAT,               intent(in)    :: delta
    integer,             intent(in)    :: max_iter
    type(nl_operator_t), intent(in)    :: op
    FLOAT,               intent(in)    :: spacing
    integer,             intent(in)    :: order
    type(mpi_grp_t),     intent(in)    :: mpi_grp

    integer :: il, np, saved_iter, ii

    call push_sub('ob_mem.ob_mem_init')

    np = intface(LEFT)%np

    ! Allocate arrays depending on the type of coefficients requested.
    select case(ob%mem_type)
    case(SAVE_CPU_TIME)
      ALLOCATE(ob%mem_coeff(np, np, 0:max_iter, NLEADS), np**2*(max_iter+1)*NLEADS)
      allocate(ob%mem_sp_coeff(0, 0, NLEADS), ob%mem_s(0, 0, 2, NLEADS), ob%sp2full_map(0))
      ob%mem_coeff = M_z0
    case(SAVE_RAM_USAGE) ! FIXME: only 2D.
      ASSERT(calc_dim.eq.2)
      allocate(ob%mem_coeff(0, 0, 0, NLEADS))
      ! sp_coeff has the size ny*nz*do**2, (np=ny*nz*do).
      ALLOCATE(ob%mem_sp_coeff(np*order, 0:max_iter, NLEADS), order*np*(max_iter+1)*NLEADS)
      ALLOCATE(ob%mem_s(np, np, 2, NLEADS), np**2*2*NLEADS)
      ALLOCATE(ob%sp2full_map(np*order), order*np)
      ! Fill sparse to full mapping matrix.
      do ii = 0, np-1
        do il = 1, order
          ob%sp2full_map(ii*order+il) = mod((il-1)*(np/order) + ii , np) + 1
        end do
      end do
    end select

    ! Get iteration parameters from input file.

    !%Variable MemoryTol
    !%Type float
    !%Default 1e-12
    !%Section Time Dependent::Open Boundaries
    !%Description
    !% Decides when to consider the memory coefficients converged.
    !%End
    call loct_parse_float(datasets_check('MemoryTol'), CNST(1e-12), mem_tolerance)
    if(mem_tolerance.le.M_ZERO) then
      write(message(1), '(a,f14.6,a)') "Input : '", mem_tolerance, "' is not a valid MemoryTol."
      message(2) = '(0 < TDTransMemTol)'
      call write_fatal(2)
    end if

    !%Variable MemoryMaxIter
    !%Type integer
    !%Default 500
    !%Section Time Dependent::Open Boundaries
    !%Description
    !% Sets the maximum iteration number to converge the memory coefficients.
    !%End
    call loct_parse_int(datasets_check('MemoryMaxIter'), 500, mem_iter)
    if(mem_iter.le.0) then
      write(message(1), '(a,i6,a)') "Input : '", mem_iter, "' is not a valid MemoryMaxIter."
      message(2) = '(0 <= TDTransMemMaxIter)'
      call write_fatal(2)
    end if

    do il = 1, NLEADS
      ! Try to read the coefficients from file
      call read_coeffs(trim(restart_dir)//'open_boundaries/', saved_iter, ob%mem_coeff(:, :, :, il),     &
        ob%mem_sp_coeff(:, :, il), ob%mem_s(:, :, :, il), calc_dim, max_iter, np, spacing, delta, op%stencil%size, &
        ob%mem_type, np*order, il)

      if (saved_iter.lt.max_iter) then ! Calculate missing coefficients.
        if (saved_iter.gt.0) then
          write(message(1),'(a,i5,a)') 'Info: Successfully loaded the first', saved_iter, &
            ' memory coefficients of '//trim(lead_name(il))//' lead.'
          call write_info(1)
        end if
        message(1) = 'Info: Calculating missing coefficients for memory term of '// &
          trim(lead_name(il))//' lead.'
        call write_info(1)

        ! Initialize progress bar.
        call loct_progress_bar(-1, max_iter+1)

        ! FIXME: the spinor index of hm%lead_h_diag is ignored here.
        ASSERT(hm%d%ispin.ne.SPINORS)
        if(saved_iter.eq.0) then 
          ! Get anchor for recursion.
          call approx_coeff0(intface(il), delta, il, hm%lead_h_diag(:, :, 1, il),              &
            hm%lead_h_offdiag(:, :, il), ob%mem_coeff(:, :, 0, il), ob%mem_sp_coeff(:, 0, il), &
            ob%mem_s(:,:,:,il), order, ob%mem_type, ob%sp2full_map, spacing)
          call loct_progress_bar(1, max_iter+1)
        end if
        ! Calculate the subsequent coefficients by the recursive relation.
        select case(ob%mem_type)
        case(SAVE_CPU_TIME)
          if(intface(il)%offdiag_invertible) then
            call calculate_coeffs(il, saved_iter+1, max_iter, delta, intface(il), hm%lead_h_diag(:, :, 1, il), &
              hm%lead_h_offdiag(:, :, il), ob%mem_coeff(:, :, :, il), spacing)
          else
            call calculate_coeffs_ni(il, saved_iter+1, max_iter, delta, intface(il), hm%lead_h_diag(:, :, 1, il), &
              hm%lead_h_offdiag(:, :, il), ob%mem_coeff(:, :, :, il))
          end if
        case(SAVE_RAM_USAGE) ! FIXME: only 2D.
          ASSERT(calc_dim.eq.2)
          call calculate_sp_coeffs(il, saved_iter+1, max_iter, delta, intface(il), hm%lead_h_diag(:, :, 1, il), &
            hm%lead_h_offdiag(:, :, il), ob%mem_sp_coeff(:, :, il), ob%mem_s(:, :, :, il), np*order,            &
            order, calc_dim, ob%sp2full_map, spacing)
        end select

        if(saved_iter.lt.max_iter) then
          message(1) = ''
          message(2) = 'Info: Writing memory coefficients of '// &
            trim(lead_name(il))//' lead.'
          call write_info(2)
          if(mpi_grp_is_root(mpi_grp)) then
            call write_coeffs(trim(restart_dir)//'open_boundaries/', ob%mem_coeff(:, :, :, il),          &
              ob%mem_sp_coeff(:, :, il), ob%mem_s(:, :, :, il), calc_dim, max_iter, intface(il)%np, spacing, &
              delta, op%stencil%size, ob%mem_type, np*order, il)
          end if
        end if
      else
        message(1) = 'Info: Successfully loaded memory coefficients from '// &
          trim(lead_name(il))//' lead.'
        call write_info(1)
      end if

    end do

    call pop_sub()
  end subroutine ob_mem_init


  ! ---------------------------------------------------------
  ! Solve for zeroth memory coefficient by truncating the continued
  ! matrix fraction. Since the coefficient must be symmetric,
  ! a symmetric inversion is used.
  subroutine approx_coeff0(intface, delta, il, diag, offdiag, coeff0, sp_coeff0, &
    mem_s, order, mem_type, mapping, spacing)
    type(interface_t), intent(in)  :: intface
    FLOAT,             intent(in)  :: delta
    integer,           intent(in)  :: il
    CMPLX,             intent(in)  :: diag(:, :)
    CMPLX,             intent(in)  :: offdiag(:, :)
    CMPLX,             intent(out) :: coeff0(:, :)
    CMPLX,             intent(out) :: sp_coeff0(:)   ! 0th coefficient in packed storage.
    CMPLX,             intent(out) :: mem_s(:, :, :) ! S & S^(-1).
    integer,           intent(in)  :: order
    integer,           intent(in)  :: mem_type
    integer,           intent(in)  :: mapping(:)     ! Mapping.
    FLOAT,             intent(in)  :: spacing

    integer            :: i, j, np
    CMPLX, allocatable :: q0(:, :)
    FLOAT              :: norm, old_norm, sp2
    CMPLX              :: h

    call push_sub('ob_mem.approx_coeff0')

    np = intface%np

    ! If we are in 1D and have only a number we can solve the equation explicitly.
    ! So check this first for faster calculation.
    if(np.eq.1) then
      h            = M_z1 + M_zI*delta*diag(1,1)
      coeff0(1, 1) = (-h + sqrt(h**2 + (M_TWO*delta*offdiag(1,1))**2)) / (M_TWO*delta**2)
    else ! We have the general case of a matrix, so solve the equation by iteration.
      ! Truncating the continued fraction is the same as iterating the equation
      !
      !     q0 = V^T * (1/(1 + i*delta*h + delta^2*q0)) * V
      !
      ! with start value 0:
      sp2 = spacing**2
      ALLOCATE(q0(np, np), np**2)

      q0(:, :) = M_z0
      old_norm = M_ZERO
      do i = 1, mem_iter
        ! Calculate 1 + i*delta*h + delta^2*coeff0_old
        q0 = M_zI*delta*diag + delta**2*q0
        do j = 1, np
          q0(j, j) = 1 + q0(j, j)
        end do

        ! Invert.
        call lalg_sym_inverter('U', np, q0)
        call matrix_symmetrize(q0, np)

        ! Apply coupling matrices.
        call apply_coupling(q0, offdiag, q0, np, il)
        call matrix_symmetric_average(q0, np)
        norm = infinity_norm(q0)
        if((abs(norm-old_norm)*sp2).lt.mem_tolerance) then
          exit
        end if
        old_norm = norm
      end do
      if(i.gt.mem_iter) then
        write(message(1), '(a,i6,a)') 'Memory coefficent for time step 0, ' &
          //trim(lead_name(il))//' lead, not converged'
        call write_warning(1)
      end if

      select case(mem_type)
      case(SAVE_CPU_TIME)
        coeff0 = q0
      case(SAVE_RAM_USAGE)
        ! Diagonalization procedure.
        mem_s(:, :, 1) = q0(:, :)
        call lalg_eigensolve_nonh(np, mem_s(:, :, 1), mem_s(:, 1, 2))
        mem_s(:, :, 2) = mem_s(:, :, 1)
        norm = lalg_inverter(np, mem_s(:, :, 2), invert=.true.)
        call make_sparse_matrix(np, order, 2, q0, mem_s, sp_coeff0, mapping)
      end select
      deallocate(q0)
    end if

    call pop_sub()
  end subroutine approx_coeff0


  ! ---------------------------------------------------------
  ! coeffs(:, :, :, 0) given, calculate the subsequent ones by
  ! the recursive relation. Since the 0th coefficiant is symmetric
  ! all subsequent are also, therefore symmetric matrix multiplications
  ! can be used.
  subroutine calculate_coeffs(il, start_iter, iter, delta, intf, diag, offdiag, coeffs, spacing)
    integer,           intent(in)    :: il
    integer,           intent(in)    :: start_iter
    integer,           intent(in)    :: iter
    FLOAT,             intent(in)    :: delta
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: diag(intf%np, intf%np)
    CMPLX,             intent(in)    :: offdiag(intf%np, intf%np)
    CMPLX,             intent(inout) :: coeffs(intf%np, intf%np, 0:iter)
    FLOAT,             intent(in)    :: spacing

    FLOAT              :: old_norm, norm, sp2
    integer            :: i, j, k, np
    CMPLX, allocatable :: tmp(:, :), tmp2(:, :), inv_offdiag(:, :)
    CMPLX, allocatable :: prefactor_plus(:, :), prefactor_minus(:, :), fr(:, :)
    CMPLX              :: tp, d2
    character          :: uplo

    call push_sub('ob_mem.calculate_coeffs')

    ASSERT(intf%offdiag_invertible)
    np  = intf%np
    sp2 = spacing**2
    d2  = TOCMPLX(delta**2, M_ZERO)

    ALLOCATE(tmp(np, np), np**2)
    ALLOCATE(tmp2(np, np), np**2)
    ALLOCATE(inv_offdiag(np, np), np**2)
    ALLOCATE(prefactor_plus(np, np), np**2)
    ALLOCATE(prefactor_minus(np, np), np**2)
    ALLOCATE(fr(np, np), np**2)

    prefactor_plus(:, :)  = M_zI*delta*diag(:, :)
    prefactor_minus(:, :) = -M_zI*delta*diag(:, :)

    do i = 1, np
      prefactor_plus(i, i)  = M_ONE + prefactor_plus(i, i)
      prefactor_minus(i, i) = M_ONE + prefactor_minus(i, i)
    end do

    ! If we are in 1D and have only a number we can calculate the coefficients explicitly.
    ! So check this first for faster calculation.
    if(np.eq.1) then
      tp = M_ONE/(prefactor_plus(1, 1) + M_TWO*d2*coeffs(1, 1, 0))
      do i = start_iter, iter ! i = start_iter to m.
        if(i.eq.1) then
          coeffs(1, 1, 1) = coeffs(1,1,0)*(prefactor_minus(1, 1) - M_TWO*d2*coeffs(1, 1, 0))*tp
        else ! i >= 2.
          ! j = 1.
          coeffs(1, 1, i) = (coeffs(1, 1, 1)+M_TWO*coeffs(1, 1, 0))*coeffs(1, 1, i-1) + &
            coeffs(1, 1, 0)*coeffs(1, 1, i-2)
          do j = 2, i-1
            coeffs(1, 1, i) = coeffs(1, 1, i) + (coeffs(1, 1, j) + &
              M_TWO*coeffs(1, 1, j-1)+coeffs(1, 1, j-2))*coeffs(1, 1, i-j)
          end do
          coeffs(1, 1, i) = coeffs(1, 1, 1)*coeffs(1, 1, i-1)/coeffs(1, 1, 0) - d2*tp*coeffs(1, 1, i)
        end if
        call loct_progress_bar(i+1, iter+1)
      end do
    else ! We have the general case of a matrix, so solve the equation by iteration.
      inv_offdiag(:, :) = offdiag(:, :)
      if(il.eq.LEFT) then
        call lalg_invert_upper_triangular(np, inv_offdiag)
      else
        call lalg_invert_lower_triangular(np, inv_offdiag)
      end if

      prefactor_plus(:, :) = prefactor_plus(:, :) + d2*coeffs(:, :, 0)
      call lalg_sym_inverter('U', np, prefactor_plus)
      call matrix_symmetrize(prefactor_plus, np)
      fr(:, :) = coeffs(:, :, 0)
      if (il.eq.LEFT) then
        uplo = 'U'
      else
        uplo = 'L'
      end if
      call lalg_trmm(np, np, uplo, 'N', 'L', M_z1, offdiag, prefactor_plus)
      call lalg_trmm(np, np, uplo, 'N', 'R', M_z1, inv_offdiag, prefactor_minus)
      call lalg_trmm(np, np, uplo, 'N', 'L', d2, inv_offdiag, fr)

      do i = start_iter, iter ! i = start_iter to m.
        tmp2(:, :) = M_z0
        ! First calculate the part without q(m) and store it in tmp2.
        ! k = 1 to k = m, but without q(m).
        do k = 1, i
          tmp(:, :) = M_TWO*coeffs(:, :, k-1)
          if(k.ne.i) tmp(:, :) = tmp(:, :) + coeffs(:, :, k)
          if(k.gt.1) tmp(:, :) = tmp(:, :) + coeffs(:, :, k-2)
          if (il.eq.LEFT) then
            call lalg_trmm(np, np,'U','N','R', M_z1, inv_offdiag, tmp)
          else
            call lalg_trmm(np, np, 'L', 'N', 'R', M_z1, inv_offdiag, tmp)
          end if
          call lalg_symm(np, np, 'R', M_z1, coeffs(:, :, i-k), tmp, M_z1, tmp2)
        end do
        call lalg_symm(np, np, 'R', M_z1, coeffs(:, :, i-1), prefactor_minus, -d2, tmp2)

        ! Now solve the equation via iteration (for now).
        ! Save some time, so set the first iteration explicitly.
        coeffs(:, :, i) = M_z0
        call lalg_gemm(np, np, np, M_z1, prefactor_plus, tmp2, M_z0, coeffs(:, :, i))
        call matrix_symmetric_average(coeffs(:, :, i), np)
        old_norm = infinity_norm(coeffs(:, :, i))

        do j = 2, mem_iter ! maxiter to converge matrix equation for q(i).
          tmp(:, :) = tmp2(:, :)
          call lalg_symm(np, np, 'L', -M_z1, coeffs(:, :, i), fr, M_z1, tmp)
          call lalg_gemm(np, np, np, M_z1, prefactor_plus, tmp, M_z0, coeffs(:, :, i))
          ! Use for numerical stability.
          call matrix_symmetric_average(coeffs(:, :, i), np)
          norm = infinity_norm(coeffs(:, :, i))
          if((abs(norm-old_norm)*sp2).lt.mem_tolerance) then
            exit
          end if
          old_norm = norm
        end do

        ! Write a warning if a coefficient is not converged.
        if(j.gt.mem_iter) then
          write(message(1), '(a,i6,a)') 'Memory coefficent for time step ', i, &
            ', '//trim(lead_name(il))//' lead, not converged.'
          call write_warning(1)
        end if

        call loct_progress_bar(i+1, iter+1)
      end do
    end if

    deallocate(tmp, tmp2, inv_offdiag)
    deallocate(prefactor_plus, prefactor_minus, fr)

    call pop_sub()
  end subroutine calculate_coeffs


  ! ---------------------------------------------------------
  ! A bit faster, uses twice the memory.
  ! But it is able to calculate the memory coefficients in the
  ! case that the offdiagonal block is not invertible also.
  ! coeffs(:, :, :, 0) given, calculate the subsequent ones by
  ! the recursive relation. Since the 0th coefficiant is symmetric
  ! all subsequent are also, therefore symmetric matrix multiplications
  ! can be used.
  subroutine calculate_coeffs_ni(il, start_iter, iter, delta, intf, diag, offdiag, coeffs)
    integer,           intent(in)    :: il
    integer,           intent(in)    :: start_iter
    integer,           intent(in)    :: iter
    FLOAT,             intent(in)    :: delta
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: diag(:, :)
    CMPLX,             intent(in)    :: offdiag(:, :)
    CMPLX,             intent(inout) :: coeffs(intf%np, intf%np, 0:iter)

    integer            :: i,j, k, np
    CMPLX, allocatable :: coeff_p(:, :, :), p_prev(:, :), tmp(:, :), tmp2(:, :)
    CMPLX, allocatable :: prefactor_plus(:, :), prefactor_minus(:, :)
    FLOAT              :: norm, old_norm

    call push_sub('ob_mem.calculate_coeffs_ni')

    np = intf%np

    ALLOCATE(coeff_p(np, np, 0:iter), np**2*(iter+1))
    ALLOCATE(p_prev(np, np), np**2)
    ALLOCATE(tmp(np, np), np**2)
    ALLOCATE(tmp2(np, np), np**2)
    ALLOCATE(prefactor_plus(np, np), np**2)
    ALLOCATE(prefactor_minus(np, np), np**2)

    coeff_p         = M_z0
    prefactor_plus  = M_zI*delta*diag(:, :)
    prefactor_minus = -M_zI*delta*diag(:, :)

    do i = 1, np
      prefactor_plus(i, i)  = M_ONE + prefactor_plus(i, i)
      prefactor_minus(i, i) = M_ONE + prefactor_minus(i, i)
    end do

    ! Calculate p_\alpha = 1/(1 + i*delta*h_\alpha + delta^2*q_\alpha)
    coeff_p(:, :, 0) = prefactor_plus + delta**2 * coeffs(:, :, 0)
    call lalg_sym_inverter('U', np, coeff_p(:, :, 0))
    call matrix_symmetrize(coeff_p(:, :, 0), np)

    ! We only need the inverse if (1 + i*delta*h_\alpha) in the
    ! following.
    call lalg_sym_inverter('U', np, prefactor_plus)
    call matrix_symmetrize(prefactor_plus, np)

    ! FIXME: this routine is not able to restart because the coeff_p
    ! matrices are required. They have to be written to a file also
    ! to make restarts work.
    ASSERT(start_iter.eq.1)

    do i = start_iter, iter
      old_norm = M_ZERO

      ! The part of the sum independent of coeff_p(:, :, i),
      ! accumulated in tmp2.
      tmp2 = M_z0
      do k = 1, i-1
        call lalg_copy(np**2, coeffs(1, :, k), tmp(:, 1))
        call lalg_axpy(np**2, M_z2, coeffs(:, 1, k-1), tmp(:, 1))
        if(k.gt.1) then
          call lalg_axpy(np**2, M_z1, coeffs(:, 1, k-2), tmp(:, 1))
        end if
        call lalg_symm(np, np, 'L', M_z1, tmp, coeff_p(:, :, i-k), M_z1, tmp2)
      end do

      ! The convergence loop.
      do j = 1, mem_iter
        call lalg_copy(np**2, coeff_p(:, 1, i), p_prev(:, 1))

        ! k = m.
        call apply_coupling(p_prev, offdiag, coeff_p(:, :, i), np, il)

        call lalg_axpy(np**2, M_z2, coeffs(:, 1, i-1), coeff_p(:, 1, i))

        if(i.gt.1) then
          call lalg_axpy(np**2, M_z1, coeffs(:, 1, i-2), coeff_p(:, 1, i))
        end if

        call lalg_symm(np, np, 'L', M_z1, coeff_p(:, :, i), coeff_p(:, :, 0), M_z0, tmp)
        call lalg_copy(np**2, tmp(:, 1), coeff_p(:, 1, i))

        ! k = 0.
        call lalg_symm(np, np, 'L', M_z1, coeffs(:, :, 0), p_prev, M_z1, coeff_p(:, :, i))

        ! Add the constant part from above, and multiply with prefactors.
        call lalg_axpy(np**2, M_z1, tmp2(:, 1), coeff_p(:, 1, i))
        call lalg_scal(np**2, TOCMPLX(-delta**2, 0), coeff_p(:, 1, i))
        call lalg_symm(np, np, 'L', M_z1, prefactor_minus, coeff_p(:, :, i-1), M_z1, coeff_p(:, :, i))
        call lalg_symm(np, np, 'L', M_z1, prefactor_plus, coeff_p(:, :, i), M_z0, tmp)
        call lalg_copy(np**2, tmp(:, 1), coeff_p(:, 1, i))

        norm = infinity_norm(coeff_p(:, :, i))
        if(abs(norm-old_norm).lt.mem_tolerance) then
          exit
        end if
        old_norm = norm
      end do

      ! Write a warning if a coefficient is not converged.
      if(j.gt.mem_iter) then
        write(message(1), '(a,i6,a)') 'Memory coefficent for time step ', i, &
          ', '//trim(lead_name(il))//' lead, not converged.'
        call write_warning(1)
      end if

      call apply_coupling(coeff_p(:, :, i), offdiag, coeffs(:, :, i), np, il)

      call loct_progress_bar(i+1, iter+1)
    end do

    message(1) = ''
    call write_info(1)

    deallocate(coeff_p, p_prev, tmp, tmp2)
    deallocate(prefactor_plus, prefactor_minus)

    call pop_sub()
  end subroutine calculate_coeffs_ni


  ! ---------------------------------------------------------
  ! Same as calculate_coeffs but with the sparse matrix
  ! sp_coeffs(:, :, 0) given, calculate the subsequent ones by
  ! the recursive relation. We can only use the (sparse) mem_q, so use the
  ! mem_q only recursive relation.
  subroutine calculate_sp_coeffs(il, start_iter, iter, delta, intf, diag, offdiag, &
    sp_coeffs, mem_s, length, order, dim, mapping, spacing)
    integer,           intent(in)    :: il
    integer,           intent(in)    :: start_iter
    integer,           intent(in)    :: iter
    FLOAT,             intent(in)    :: delta
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: diag(:, :)
    CMPLX,             intent(in)    :: offdiag(:, :)
    CMPLX,             intent(inout) :: sp_coeffs(1:length, 0:iter)
    CMPLX,             intent(in)    :: mem_s(intf%np, intf%np, 2)
    integer,           intent(in)    :: length
    integer,           intent(in)    :: order
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: mapping(:)   ! the mapping
    FLOAT,             intent(in)    :: spacing

    integer            :: i,j, k, np
    character          :: uplo
    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :), inv_offdiag(:, :)
    CMPLX, allocatable :: prefactor_plus(:, :), prefactor_minus(:, :)
    CMPLX, allocatable :: sp_tmp(:)
    FLOAT              :: old_norm, norm, sp2

    call push_sub('ob_mem.calculate_sp_coeffs')
    ASSERT(intf%offdiag_invertible)

    np  = intf%np
    sp2 = spacing**2

    ALLOCATE(tmp1(np, np), np**2)
    ALLOCATE(tmp2(np, np), np**2)
    ALLOCATE(tmp3(np, np), np**2)
    ALLOCATE(inv_offdiag(np, np), np**2)
    ALLOCATE(prefactor_plus(np, np), np**2)
    ALLOCATE(prefactor_minus(np, np), np**2)
    ALLOCATE(sp_tmp(length), length)

    prefactor_plus  = M_zI*delta*diag(:, :)
    prefactor_minus = -M_zI*delta*diag(:, :)

    do i = 1, np
      prefactor_plus(i, i)  = M_ONE + prefactor_plus(i, i)
      prefactor_minus(i, i) = M_ONE + prefactor_minus(i, i)
    end do
    inv_offdiag(:, :) = offdiag(:, :)
    if (il.eq.LEFT) then
      call lalg_invert_upper_triangular(np, inv_offdiag)
    else
      call lalg_invert_lower_triangular(np, inv_offdiag)
    end if

    call lalg_sym_inverter('U', np, prefactor_plus)
    call matrix_symmetrize(prefactor_plus, np)
    select case(il)
    case(LEFT)
      uplo = 'U'
    case(RIGHT)
      uplo = 'L'
    end select
    call lalg_trmm(np, np, uplo, 'N', 'L', M_z1, offdiag, prefactor_plus)
    call lalg_trmm(np, np, uplo, 'N', 'R', M_z1, inv_offdiag, prefactor_minus)

    do i = start_iter, iter
      sp_coeffs(:, i) = M_z0
      old_norm = M_ZERO

      do j = 1, mem_iter
        tmp2(:, :) = M_z0
        ! k = 0 to k = m.
        do k = 0, i
          tmp1(:, :) = M_z0
          sp_tmp = sp_coeffs(:, k)
          if(k.gt.0) then
            sp_tmp = sp_tmp + M_TWO*sp_coeffs(:, k-1)
          end if
          if(k.gt.1) then
            sp_tmp = sp_tmp + sp_coeffs(:, k-2)
          end if
          call make_full_matrix(np, order, dim, sp_tmp, mem_s, tmp1, mapping)
          if (il.eq.LEFT) then
            call lalg_trmm(np, np, 'U', 'N', 'R', M_z1, inv_offdiag, tmp1)
          else
            call lalg_trmm(np, np, 'L', 'N', 'R', M_z1, inv_offdiag, tmp1)
          end if
          call make_full_matrix(np, order, dim, sp_coeffs(:, i-k), mem_s, tmp3, mapping)
          call lalg_symm(np, np, 'R', M_z1, tmp3, tmp1, M_z1, tmp2)
        end do
        call make_full_matrix(np, order, dim, sp_coeffs(:, i-1), mem_s, tmp3, mapping)
        call lalg_symm(np, np, 'R', M_z1, tmp3, prefactor_minus, TOCMPLX(-delta**2, M_ZERO), tmp2)
        call lalg_gemm(np, np, np, M_z1, prefactor_plus, tmp2, M_z0, tmp1)
        ! Use for numerical stability.
        call matrix_symmetrize(tmp1, np)
        call make_sparse_matrix(np, order, dim, tmp1, mem_s, sp_coeffs(:, i), mapping)
        norm = infinity_norm(tmp1)
        if((abs(norm-old_norm)*sp2).lt.mem_tolerance) then
          exit
        end if
        old_norm = norm
      end do
      ! Write a warning if a coefficient is not converged.
      if(j.gt.mem_iter) then
        write(message(1), '(a,i6,a)') 'Memory coefficent for time step ', i, &
          ', '//trim(lead_name(il))//' lead, not converged.'
        call write_warning(1)
      end if

      call loct_progress_bar(i+1, iter+1)
    end do

    message(1) = ''
    call write_info(1)

    deallocate(tmp1, tmp2, tmp3, inv_offdiag)
    deallocate(prefactor_plus, prefactor_minus)
    deallocate(sp_tmp)

    call pop_sub()
  end subroutine calculate_sp_coeffs


  ! ---------------------------------------------------------
  ! Write memory coefficients to file.
  ! FIXME: this routine writes compiler and/or system dependent binary files.
  ! It should be chanegd to a platform independent format.
  subroutine write_coeffs(dir, coeffs, sp_coeffs, mem_s, dim, iter, np, &
    spacing, delta, op_n, mem_type, length, il)
    character(len=*), intent(in) :: dir
    CMPLX,            intent(in) :: coeffs(np, np, 0:iter)    ! Saved coefficients.
    CMPLX,            intent(in) :: sp_coeffs(length, 0:iter) ! saved coefficients.
    CMPLX,            intent(in) :: mem_s(np, np, 2)          ! Diagonalization matrix s.
    integer,          intent(in) :: dim                       ! Dimension.
    integer,          intent(in) :: iter                      ! Number of coefficients.
    integer,          intent(in) :: np                        ! Number of points in the interface.
    FLOAT,            intent(in) :: spacing                   ! Grid spacing.
    FLOAT,            intent(in) :: delta                     ! Timestep.
    integer,          intent(in) :: op_n                      ! Number of operator points.
    integer,          intent(in) :: mem_type                  ! SAVE_CPU_TIME or SAVE_RAM_USAGE.
    integer,          intent(in) :: length                    ! Length of the packed array.
    integer,          intent(in) :: il                        ! Which lead.

    integer :: ntime, j, iunit

    call push_sub('ob_mem.write_coeffs')

    call io_mkdir(dir, is_tmp=.true.)
    iunit = io_open(trim(dir)//trim(lead_name(il)), action='write', form='unformatted', is_tmp=.true.)
    if(iunit.lt.0) then
      message(1) = 'Cannot write memory coefficents to file.'
      call write_warning(1)
      call io_close(iunit)
      call pop_sub(); return
    end if

    ! Write numerical parameters.
    write(iunit) dim
    write(iunit) iter
    write(iunit) np
    write(iunit) spacing
    write(iunit) delta
    write(iunit) op_n
    write(iunit) mem_type

    ! Write matrices.
    select case(mem_type)
    case(SAVE_CPU_TIME)
      do ntime = 0, iter
        do j = 1, np
          write(iunit) coeffs(j, j:np, ntime)
        end do
      end do
    case(SAVE_RAM_USAGE) ! FIXME: only 2D.
      ASSERT(dim.eq.2)
      write(iunit) mem_s(:, :, 1)
      do ntime = 0, iter
        write(iunit) sp_coeffs(1:length, ntime)
      end do
    end select

    call io_close(iunit)

    call pop_sub()
  end subroutine write_coeffs


  ! ---------------------------------------------------------
  ! Read memory coefficients from file.
  ! FIXME: this routine reads compiler and/or system dependent binary files.
  ! It should be chanegd to a platform independent format.
  subroutine read_coeffs(dir, s_iter, coeffs, sp_coeffs, mem_s, dim, iter, &
    np, spacing, delta, op_n, mem_type, length, il)
    character(len=*), intent(in)    :: dir
    integer,          intent(out)   :: s_iter                    ! Number of saved coefficients.
    CMPLX,            intent(inout) :: coeffs(np, np, 0:iter)    ! Saved coefficients.
    CMPLX,            intent(inout) :: sp_coeffs(length, 0:iter) ! Saved coefficients.
    CMPLX,            intent(out)   :: mem_s(np, np, 2)          ! Diagonalization matrices.
    integer,          intent(in)    :: dim                       ! Dimension of the problem.
    integer,          intent(in)    :: iter                      ! Number of coefficients.
    integer,          intent(in)    :: np                        ! Number of points in the interface.
    FLOAT,            intent(in)    :: spacing                   ! Spacing.
    FLOAT,            intent(in)    :: delta                     ! Timestep.
    integer,          intent(in)    :: op_n                      ! Number of operator points.
    integer,          intent(in)    :: mem_type                  ! Which type of coefficient.
    integer,          intent(in)    :: length                    ! Length of the packed array.
    integer,          intent(in)    :: il                        ! Which lead.

    integer :: ntime, j, iunit, s_dim, s_np, s_op_n, s_mem_type
    FLOAT   :: s_spacing, s_delta, det

    call push_sub('ob_mem.read_coeffs')

    s_iter = 0

    ! Try to open file.
    iunit = io_open(trim(dir)//trim(lead_name(il)), action='read', &
      status='old', die=.false., is_tmp=.true., form='unformatted')
    if(iunit.lt.0) then
      call io_close(iunit)
      call pop_sub(); return
    end if

    ! Now read the data.
    read(iunit) s_dim
    read(iunit) s_iter
    read(iunit) s_np
    read(iunit) s_spacing
    read(iunit) s_delta
    read(iunit) s_op_n
    read(iunit) s_mem_type

    ! Check if numerical parameters of saved coefficients match
    ! current parameter set.
    if((s_dim.eq.dim) .and. (s_np.eq.np) .and. (s_op_n.eq.op_n) &
      .and. (s_spacing.eq.spacing) .and. (s_delta.eq.delta)     &
      .and. (s_mem_type.eq.mem_type)) then
      ! Read the coefficients.
      if (mem_type.eq.SAVE_CPU_TIME) then ! Full (upper half) matrices.
        do ntime = 0, min(iter, s_iter)
          do j = 1, np
            read(iunit) coeffs(j, j:np, ntime)
            coeffs(j:np, j, ntime) = coeffs(j, j:np, ntime)
          end do
        end do
      else ! Packed matrices (FIXME: yet only 2D).
        ASSERT(dim.eq.2)
        read(iunit) mem_s(:, :, 1)
        mem_s(:, :, 2) = mem_s(1:np, 1:np, 1)
        det = lalg_inverter(np, mem_s(:, :, 2), invert=.true.)
        do ntime = 0, min(iter, s_iter)
          read(iunit) sp_coeffs(1:length, ntime)
        end do
      end if
    else
      s_iter = 0
    end if

    call io_close(iunit)

    call pop_sub()
  end subroutine read_coeffs


  ! ---------------------------------------------------------
  ! Calculate amount of (machine) memory required for
  ! the memory term (in MBytes).
  ! Does not consider spin at the moment.
  FLOAT function mbytes_memory_term(iter, np, nl, st, mem_type, order)
    integer,        intent(in) :: iter
    integer,        intent(in) :: np(nl)
    integer,        intent(in) :: nl
    type(states_t), intent(in) :: st
    integer,        intent(in) :: mem_type
    integer,        intent(in) :: order

    integer :: il

    call push_sub('ob_mem.mbytes_memory_term')

    mbytes_memory_term = M_ZERO
    do il = 1, nl
      select case(mem_type)
      case(SAVE_CPU_TIME) ! Lots of memory, but faster.
        mbytes_memory_term = mbytes_memory_term + iter*M_TWO*REAL_PRECISION* &
          (st%nst*np(il)+np(il)**2) / M_TWO**20
      case(SAVE_RAM_USAGE) ! Less memory, but slower.
        mbytes_memory_term = mbytes_memory_term + iter*M_TWO*REAL_PRECISION* &
          (st%nst*np(il)+np(il)*order) / M_TWO**20
      end select
    end do

    call pop_sub()
  end function mbytes_memory_term


  ! ---------------------------------------------------------
  ! Create sparse matrix out of the full matrix.
  ! Only 2D case.
  ! sp = s^(-1).m.s
  subroutine make_sparse_matrix(np, order, dim, m, s, sp, mapping)
    integer, intent(in)  :: np          ! Matrix size.
    integer, intent(in)  :: order       ! Derivative order.
    integer, intent(in)  :: dim         ! Dimension
    CMPLX,   intent(in)  :: m(np, np)   ! Full matrix.
    CMPLX,   intent(in)  :: s(np, np, 2)! Diagonalization matrices.
    CMPLX,   intent(out) :: sp(:)       ! Sparse matrix.
    integer, intent(in)  :: mapping(:)  ! Mapping.

    integer            :: id, i
    CMPLX, allocatable :: tmp(:, :), tmp2(:, :)

    call push_sub('ob_mem.make_sparse_matrix')

    ASSERT(dim.eq.2)

    ALLOCATE(tmp(np, np), np**2)
    ALLOCATE(tmp2(np, np), np**2)

    ! Now calculate [S^(-1)*sp*S] to get the sparse matrix.
    ! FIXME: s, inv_s and sp are symmetric.
    call lalg_gemm(np, np, np, M_z1, s(:, :, 2), m, M_z0, tmp)
    call lalg_gemm(np, np, np, M_z1, tmp, s(:, :, 1), M_z0, tmp2)

    select case(dim)
    case(2)
      ! For each row.
      do i = 0, np-1
        do id = 1, order
          sp(i*order+id) = tmp2(i+1, mapping(i*order+id))
        end do
      end do
    case default
      message(1) = 'Sparse matrices for memory coefficients not supported'
      message(2) = 'in 1D or 3D.'
      call write_fatal(2)
    end select

    deallocate(tmp, tmp2)

    call pop_sub()
  end subroutine make_sparse_matrix


  ! ---------------------------------------------------------
  ! Fast multiplication of the sparse matrix with a dense matrix
  ! res = sp*m.
  subroutine mul_sp(np, order, dim, sp, m, mapping, res)
    integer, intent(in)  :: np          ! matrix size (np**2)
    integer, intent(in)  :: order       ! derivative order
    integer, intent(in)  :: dim         ! dimension
    CMPLX,   intent(in)  :: sp(:)       ! sparse matrix (np*order)
    CMPLX,   intent(in)  :: m(np, np)   ! the matrix to multipliy
    integer, intent(in)  :: mapping(:)  ! the mapping
    CMPLX,   intent(out) :: res(np, np) ! the full matrix

    integer :: ii,ij,ik, ijo
    CMPLX   :: tmp

    call push_sub('ob_mem.mul_sp')

    do ij=1, np
      do ii=1, np
        tmp = M_z0
        ijo = (ij-1)*order
        do ik=1,order
          tmp = tmp + sp(ijo+ik)*m(mapping(ijo+ik),ii)
        end do
        res(ij,ii) = tmp
      end do
    end do

    call pop_sub()
  end subroutine mul_sp


  ! ---------------------------------------------------------
  ! Create the original matrix out of the sparse matrix.
  ! m = s.sp.s^(-1)
  subroutine make_full_matrix(np, order, dim, sp, s, m, mapping)
    integer, intent(in)  :: np          ! matrix size (np**2)
    integer, intent(in)  :: order       ! derivative order
    integer, intent(in)  :: dim         ! dimension
    CMPLX,   intent(in)  :: sp(:)       ! sparse matrix (np*order)
    CMPLX,   intent(in)  :: s(np, np, 2)! diagonalization matrix
    CMPLX,   intent(out) :: m(np, np)   ! the full matrix
    integer, intent(in)  :: mapping(:)  ! the mapping

    CMPLX, allocatable :: tmp(:, :)
    integer            :: id, i

    call push_sub('td_transport.make_full_matrix')

    ALLOCATE(tmp(np,np),np**2)

    ! make a full matrix of the packed storage form
    ! the sparse matrix has a (symmetric) banded form
    m = M_z0
    select case(dim)
    case(2)
      do i=0, np-1
        do id=1, order
          m(i+1,mapping(i*order+id)) = sp(i*order+id)
        end do
      end do
    case(3)
      ! FIXME (not yet done)
    end select
    ! now calculate [S*sp*S^(-1)] to get the dense matrix
    call mul_sp(np, order, dim, sp, s(:,:,2), mapping, tmp)
!    call lalg_gemm(np,np,np,M_z1,m,s(:,:,2),M_z0,tmp)
    call lalg_gemm(np,np,np,M_z1,s(:,:,1),tmp,M_z0,m)
    call matrix_symmetric_average(m, np)

    deallocate(tmp)
    call pop_sub()
  end subroutine make_full_matrix


  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine ob_mem_end(ob)
    type(ob_terms_t), intent(inout) :: ob

    call push_sub('ob_mem.ob_mem_end')

    DEALLOCATE(ob%mem_coeff)
    DEALLOCATE(ob%mem_sp_coeff)
    DEALLOCATE(ob%sp2full_map)
    DEALLOCATE(ob%mem_s)

    call pop_sub()
  end subroutine ob_mem_end
end module ob_mem_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
