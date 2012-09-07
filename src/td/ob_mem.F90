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

!> Calculation of the memory coefficients for the modified
!! Crank-Nicholson propagator.

#include "global.h"

module ob_mem_m
  use datasets_m
  use global_m
  use hamiltonian_m
  use io_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use math_m
  use messages_m
  use mpi_m
  use nl_operator_m
  use ob_interface_m
  use ob_lead_m
  use ob_terms_m
  use parser_m
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
  !> Allocate (machine) memory for the memory coefficents and
  !! calculate them.
  subroutine ob_mem_init(intf, hm, ob, delta, max_iter, op, spacing, order, mpi_grp)
    type(interface_t),   intent(in)    :: intf(:)
    type(hamiltonian_t), intent(in)    :: hm
    type(ob_terms_t),    intent(inout) :: ob
    FLOAT,               intent(in)    :: delta
    integer,             intent(in)    :: max_iter
    type(nl_operator_t), intent(in)    :: op
    FLOAT,               intent(in)    :: spacing
    integer,             intent(in)    :: order
    type(mpi_grp_t),     intent(in)    :: mpi_grp

    integer :: il, np, saved_iter, ii, ij

    PUSH_SUB(ob_mem_init)

      do il=1, NLEADS
        ASSERT(intf(il)%np_intf .eq. intf(il)%np_uc)
      end do

    ! Allocate arrays depending on the type of coefficients requested.
    select case(ob%mem_type)
    case(SAVE_CPU_TIME)
      do il = 1, NLEADS
        np = intf(il)%np_intf
        SAFE_ALLOCATE(ob%lead(il)%q(1:np, 1:np, 0:max_iter))
        nullify(ob%lead(il)%q_sp, ob%lead(il)%q_s, ob%lead(il)%sp2full_map)
        ob%lead(il)%q = M_z0
      end do
    case(SAVE_RAM_USAGE) ! FIXME: only 2D.
      ASSERT(op%mesh%sb%dim .eq. 2)
      ! sp_coeff has the size ny*nz*do**2, (np=ny*nz*do).
      do il = 1, NLEADS
        np = intf(il)%np_intf
        SAFE_ALLOCATE(ob%lead(il)%q_sp(1:np*order, 0:max_iter))
        SAFE_ALLOCATE(ob%lead(il)%q_s(1:np, 1:np, 1:2))
        SAFE_ALLOCATE(ob%lead(il)%sp2full_map(1:np*order))
        ob%lead(il)%q_sp = M_z0
        ob%lead(il)%q_s  = M_z0
        ob%lead(il)%sp2full_map = M_z0
        nullify(ob%lead(il)%q)
        ! Fill sparse to full mapping matrix.
        do ii = 0, np-1
          do ij = 1, order
            ob%lead(il)%sp2full_map(ii*order+ij) = mod((ij-1)*(np/order) + ii , np) + 1
          end do
        end do
      end do
    end select

    ! Get iteration parameters from input file.

    !%Variable MemoryTol
    !%Type float
    !%Default 1e-12
    !%Section Time-Dependent::Open Boundaries
    !%Description
    !% Decides when to consider the memory coefficients converged.
    !%End
    call parse_float(datasets_check('MemoryTol'), CNST(1e-12), mem_tolerance)
    if(mem_tolerance .le. M_ZERO) then
      write(message(1), '(a,f14.6,a)') "Input : '", mem_tolerance, "' is not a valid MemoryTol."
      message(2) = '(0 < TDTransMemTol)'
      call messages_fatal(2)
    end if

    !%Variable MemoryMaxIter
    !%Type integer
    !%Default 500
    !%Section Time-Dependent::Open Boundaries
    !%Description
    !% Sets the maximum iteration number to converge the memory coefficients.
    !%End
    call parse_integer(datasets_check('MemoryMaxIter'), 500, mem_iter)
    if(mem_iter .le. 0) then
      write(message(1), '(a,i6,a)') "Input : '", mem_iter, "' is not a valid MemoryMaxIter."
      message(2) = '(0 <= TDTransMemMaxIter)'
      call messages_fatal(2)
    end if

    do il = 1, NLEADS
      np = intf(il)%np_intf
      ! Try to read the coefficients from file
      call read_coeffs(trim(restart_dir)//'open_boundaries/', saved_iter, ob, hm, intf(il), &
                       op%mesh%sb%dim, max_iter, spacing, delta, op%stencil%size, order)

      if (saved_iter .lt. max_iter) then ! Calculate missing coefficients.
        if (saved_iter .gt. 0) then
          write(message(1),'(a,i5,a)') 'Info: Successfully loaded the first', saved_iter, &
            ' memory coefficients of '//trim(lead_name(il))//' lead.'
          call messages_info(1)
        end if
        message(1) = 'Info: Calculating missing coefficients for memory term of '// &
          trim(lead_name(il))//' lead.'
        call messages_info(1)

        ! Initialize progress bar.
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, max_iter+1)

        ! FIXME: the spinor index of hm%lead(il)%h_diag is ignored here.
        ASSERT(hm%d%ispin .ne. SPINORS)
        if(saved_iter .eq. 0) then 
          ! Get anchor for recursion.
          select case(ob%mem_type)
          case(SAVE_CPU_TIME)
            call approx_coeff0(intf(il), delta, hm%lead(il)%h_diag(:, :, 1),   &
              hm%lead(il)%h_offdiag(:, :), ob%lead(il)%q(:, :, 0))
          case(SAVE_RAM_USAGE) ! FIXME: only 2D.
            ASSERT(op%mesh%sb%dim .eq. 2)
            call approx_sp_coeff0(intf(il), delta, hm%lead(il)%h_diag(:, :, 1), &
                 hm%lead(il)%h_offdiag(:, :), ob%lead(il)%q_sp(:, 0), &
                 ob%lead(il)%q_s(:,:,:), order, ob%lead(il)%sp2full_map)
          end select
          if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(1, max_iter+1)
        end if
        ! Calculate the subsequent coefficients by the recursive relation.
        select case(ob%mem_type)
        case(SAVE_CPU_TIME)
          call calculate_coeffs_ni(saved_iter+1, max_iter, delta, intf(il), hm%lead(il)%h_diag(:, :, 1), &
            hm%lead(il)%h_offdiag(:, :), ob%lead(il)%q(:, :, :))
        case(SAVE_RAM_USAGE) ! FIXME: only 2D.
          ASSERT(op%mesh%sb%dim .eq. 2)
          call apply_coupling(ob%lead(il)%q(:, :, 0), hm%lead(il)%h_offdiag(:, :),&
              ob%lead(il)%q(:, :, 0), np, il)
          call calculate_sp_coeffs(saved_iter+1, max_iter, delta, intf(il), hm%lead(il)%h_diag(:, :, 1), &
            hm%lead(il)%h_offdiag(:, :), ob%lead(il)%q_sp(:, :), ob%lead(il)%q_s(:, :, :), np*order,     &
            order, op%mesh%sb%dim, ob%lead(il)%sp2full_map, spacing)
        end select

        if(saved_iter .lt. max_iter) then
          message(1) = ''
          message(2) = 'Info: Writing memory coefficients of '// &
            trim(lead_name(il))//' lead.'
          call messages_info(2)
          if(mpi_grp_is_root(mpi_grp)) then
            call write_coeffs(trim(restart_dir)//'open_boundaries/', ob, hm, intf(il),     &
              op%mesh%sb%dim, max_iter, spacing, delta, op%stencil%size, order)
          end if
        end if
      else
        message(1) = 'Info: Successfully loaded memory coefficients from '// &
          trim(lead_name(il))//' lead.'
        call messages_info(1)
      end if

      if(ob%mem_type .eq. SAVE_CPU_TIME) then
        do ii=0, max_iter
          call apply_coupling(ob%lead(il)%q(:, :, ii), hm%lead(il)%h_offdiag(:, :),&
              ob%lead(il)%q(:, :, ii), np, il)
        end do
      end if

    end do

    POP_SUB(ob_mem_init)
  end subroutine ob_mem_init


  ! ---------------------------------------------------------
  !> Solve for zeroth memory coefficient by truncating the continued
  !! matrix fraction. Since the coefficient must be symmetric (non-magnetic),
  !! a symmetric inversion is used.
  subroutine approx_coeff0(intf, delta, diag, offdiag, coeff0)
    type(interface_t), intent(in)  :: intf
    FLOAT,             intent(in)  :: delta
    CMPLX,             intent(in)  :: diag(:, :)
    CMPLX,             intent(in)  :: offdiag(:, :)
    CMPLX,             intent(out) :: coeff0(:, :)

    integer            :: iter, ip, np
    CMPLX, allocatable :: q0(:, :)
    FLOAT              :: norm, old_norm, d2
    CMPLX              :: hh

    PUSH_SUB(approx_coeff0)

    np = intf%np_intf
    d2 = delta**2

    ! If we are in 1D and have only a number we can solve the equation explicitly.
    ! So check this first for faster calculation.
    if(np .eq. 1) then
      hh = M_z1 + M_zI*delta*diag(1,1)
      if(infinity_norm(offdiag) > M_ZERO) then
        coeff0(1, 1) = (-hh + sqrt(hh**2 + d2*(M_TWO*offdiag(1,1))**2)) / (M_TWO*d2*offdiag(1,1)**2)
      else
        message(1) = 'Error in approx_coeff0:'
        message(2) = 'Off-diagonal term of Hamiltonian must not be the zero matrix!'
        call messages_fatal(2)
      end if
    else ! We have the general case of a matrix, so solve the equation by iteration.
      ! Truncating the continued fraction is the same as iterating the equation
      !
      !     p0 = (1/(1 + i*delta*h + delta^2*q0))
      !     q0 = V^T * p0 * V
      !
      ! with start value 0:
      SAFE_ALLOCATE(q0(1:np, 1:np))

      q0(:, :) = M_z0

      old_norm = M_ZERO
      do iter = 1, mem_iter
        ! Calculate 1 + i*delta*h + delta^2*coeff0_old
        coeff0(:, :) = M_zI*delta*diag(:, :) + d2*q0(:, :)
        do ip = 1, np
          coeff0(ip, ip) = M_ONE + coeff0(ip, ip)
        end do

        ! Invert.
        call lalg_sym_inverter('U', np, coeff0)
        call matrix_symmetrize(coeff0, np)
        norm = infinity_norm(coeff0)
        if((abs(old_norm/norm-M_ONE)) .lt. mem_tolerance) then
          exit
        end if
        ! Apply coupling matrices.
        call apply_coupling(coeff0, offdiag, q0, np, intf%il)
        call matrix_symmetric_average(q0, np)
        old_norm = norm
      end do
      if(iter .gt. mem_iter) then
        write(message(1), '(a,i6,a)') 'Memory coefficent for time step 0, ' &
          //trim(lead_name(intf%il))//' lead, not converged'
        call messages_warning(1)
      end if

      SAFE_DEALLOCATE_A(q0)
    end if

    POP_SUB(approx_coeff0)
  end subroutine approx_coeff0


  ! ---------------------------------------------------------
  !> Solve for zeroth memory coefficient by truncating the continued
  !! matrix fraction. Since the coefficient must be symmetric (non-magnetic),
  !! a symmetric inversion is used. Sparse version
  subroutine approx_sp_coeff0(intf, delta, diag, offdiag, sp_coeff0, mem_s, order, mapping)
    type(interface_t), intent(in)  :: intf
    FLOAT,             intent(in)  :: delta
    CMPLX,             intent(in)  :: diag(:, :)
    CMPLX,             intent(in)  :: offdiag(:, :)
    CMPLX,             intent(out) :: sp_coeff0(:)   ! 0th coefficient in packed storage.
    CMPLX,             intent(out) :: mem_s(:, :, :) ! S & S^(-1).
    integer,           intent(in)  :: order
    integer,           intent(in)  :: mapping(:)     ! Mapping.

    integer            :: iter, ip, np
    CMPLX, allocatable :: q0(:, :)
    FLOAT              :: norm, old_norm, d2

    PUSH_SUB(approx_sp_coeff0)

    np = intf%np_intf
    d2 = delta**2

    ! Truncating the continued fraction is the same as iterating the equation
    !
    !     p0 = (1/(1 + i*delta*h + delta^2*q0))
    !     q0 = V^T * p0 * V
    !
    ! with start value 0:
    SAFE_ALLOCATE(q0(1:np, 1:np))

    q0(:, :) = M_z0

    old_norm = M_ZERO
    do iter = 1, mem_iter
      ! Calculate 1 + i*delta*h + delta^2*coeff0_old
      q0(:, :) = M_zI*delta*diag(:, :) + d2*q0(:, :)
      do ip = 1, np
        q0(ip, ip) = M_ONE + q0(ip, ip)
      end do

      ! Invert.
      call lalg_sym_inverter('U', np, q0)
      call matrix_symmetrize(q0, np)
      norm = infinity_norm(q0)
      if((abs(old_norm/norm-M_ONE)) .lt. mem_tolerance) then
        exit
      end if
      ! Apply coupling matrices.
      call apply_coupling(q0, offdiag, q0, np, intf%il)
      call matrix_symmetric_average(q0, np)
      old_norm = norm
    end do
    if(iter .gt. mem_iter) then
      write(message(1), '(a,i6,a)') 'Memory coefficent for time step 0, ' &
        //trim(lead_name(intf%il))//' lead, not converged'
      call messages_warning(1)
    end if

    ! Diagonalization procedure.
    mem_s(:, :, 1) = q0(:, :)
    call lalg_eigensolve_nonh(np, mem_s(:, :, 1), mem_s(:, 1, 2))
    mem_s(:, :, 2) = mem_s(:, :, 1)
    norm = lalg_inverter(np, mem_s(:, :, 2), invert=.true.)
    call make_sparse_matrix(np, order, 2, q0, mem_s, sp_coeff0, mapping)

    SAFE_DEALLOCATE_A(q0)

    POP_SUB(approx_sp_coeff0)
  end subroutine approx_sp_coeff0


  ! ---------------------------------------------------------
  !> New version of calculating he memory coefficients.
  !!
  !! It now calculates and stores the p-coefficients.
  !! This method is more general and can even continue if
  !! the matrices are not invertable. It uses only the memory
  !! which is needed by the q-matrices (defined by \f$q=v^T.p.v\f$).
  !! These coefficients have to be calculated after this routine
  !! is finished and the coefficients are stored on the disk.
  !! This version cannot handle magnetic fields in the leads, 
  !! as the assumption of symmetric matrices for multiplications
  !! is not valid anymore.
  subroutine calculate_coeffs_ni(start_iter, iter, delta, intf, diag, offdiag, coeffs)
    integer,           intent(in)    :: start_iter
    integer,           intent(in)    :: iter
    FLOAT,             intent(in)    :: delta
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: diag(:, :)
    CMPLX,             intent(in)    :: offdiag(:, :)
    CMPLX,             intent(inout) :: coeffs(intf%np_intf, intf%np_intf, 0:iter)

    integer            :: ip, ii, ji, ki, np
    CMPLX, allocatable :: tmp(:, :), tmp2(:, :), m0(:, :), m_l(:, :), m_r(:, :)
    FLOAT              :: norm, old_norm

    PUSH_SUB(calculate_coeffs_ni)

    np = intf%np_intf

    SAFE_ALLOCATE( tmp(1:np, 1:np))
    SAFE_ALLOCATE(tmp2(1:np, 1:np))
    SAFE_ALLOCATE(  m0(1:np, 1:np))
    SAFE_ALLOCATE( m_l(1:np, 1:np))
    SAFE_ALLOCATE( m_r(1:np, 1:np))

    m0 = M_z0
    tmp(:, :) = -M_zI*delta*diag(:, :)
    do ip = 1, np
      tmp(ip, ip) = M_ONE + tmp(ip, ip)
    end do

    call lalg_symm(np, np, 'L', M_z1, coeffs(:, :, 0), tmp, M_z0, m0)

    m_l(:, :) = delta*coeffs(:, :, 0)
    m_r(:, :) = m_l(:, :)
    if(mod(intf%il+1,2)+1 .eq. 1) then
      call lalg_trmm(np, np, 'U', 'N', 'R', M_z1, offdiag, m_l)
      call lalg_trmm(np, np, 'U', 'T', 'L', M_z1, offdiag, m_r)
    else
      call lalg_trmm(np, np, 'L', 'N', 'R', M_z1, offdiag, m_l)
      call lalg_trmm(np, np, 'L', 'T', 'L', M_z1, offdiag, m_r)
    end if

    do ii = start_iter, iter
      ! The part of the sum independent of coeff_p(:, :, ii),
      ! accumulated in tmp2.
      tmp2 = M_z0
      if(ii .gt. 1) then
        tmp(:, :) = M_TWO*coeffs(:, :, ii-1) + coeffs(:, :, ii-2)
      else
        tmp(:, :) = M_TWO*coeffs(:, :, ii-1)
      end if
      call lalg_symm(np, np, 'L', M_z1, tmp, m_r, M_z0, tmp2)
      
      do ki = 1, ii-1       
        if(ki .gt. 1) then
          tmp(:, :) = coeffs(:, :, ki) + M_TWO*coeffs(:, :, ki-1) + coeffs(:, :, ki-2)
        else
          tmp(:, :) = coeffs(:, :, ki) + M_TWO*coeffs(:, :, ki-1)
        end if
        if(mod(intf%il+1,2)+1.eq.1) then
          call lalg_trmm(np, np, 'U', 'T', 'R', M_z1, offdiag, tmp)
        else
          call lalg_trmm(np, np, 'L', 'T', 'R', M_z1, offdiag, tmp)
        end if
        call lalg_symm(np, np, 'R', TOCMPLX(delta, M_ZERO), coeffs(:, :, ii-ki), tmp, M_z1, tmp2)
      end do

      call lalg_symm(np, np, 'R', M_z1, coeffs(:, :, ii-1), m0, M_z0, tmp)
      call lalg_gemm(np, np, np, -M_z1, m_l, tmp2, M_z1, tmp)
      call matrix_symmetric_average(tmp, np)

      
      ! initialize coefficient
      coeffs(:, :, ii) = tmp(:, :)
      old_norm = infinity_norm(coeffs(:, :, ii))
      ! The convergence loop.
      do ji = 1, mem_iter
        call lalg_symm(np, np, 'L', M_z1, coeffs(:, :, ii), m_r, M_z0, tmp2)
!        call lalg_gemm(np, np, np, M_z1, coeffs(:, :, ii), m_r, M_z0, tmp2)
        call lalg_gemm(np, np, np, M_z1, m_l, tmp2, M_z0, coeffs(:, :, ii))
        coeffs(:, :, ii) = tmp(:, :) - coeffs(:, :, ii)
        call matrix_symmetric_average(coeffs(:, :, ii), np)

        norm = infinity_norm(coeffs(:, :, ii))
        if(abs(old_norm/norm-M_ONE) .lt. mem_tolerance) then
          exit
        end if
        old_norm = norm
      end do

      ! Write a warning if a coefficient is not converged.
      if(ji .gt. mem_iter) then
        write(message(1), '(a,i6,a)') 'Memory coefficent for time step ', ii, &
          ', '//trim(lead_name(intf%il))//' lead, not converged.'
        call messages_warning(1)
      end if

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ii+1, iter+1)
    end do

    message(1) = ''
    call messages_info(1)

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp2)
    SAFE_DEALLOCATE_A(m0)
    SAFE_DEALLOCATE_A(m_l)
    SAFE_DEALLOCATE_A(m_r)

    POP_SUB(calculate_coeffs_ni)
  end subroutine calculate_coeffs_ni


  ! ---------------------------------------------------------
  !> Same as calculate_coeffs but with the sparse matrix
  !! sp_coeffs(:, :, 0) given, calculate the subsequent ones by
  !! the recursive relation. We can only use the (sparse) mem_q, so use the
  !! mem_q only recursive relation.
  subroutine calculate_sp_coeffs(start_iter, iter, delta, intf, diag, offdiag, &
    sp_coeffs, mem_s, length, order, dim, mapping, spacing)
    integer,           intent(in)    :: start_iter
    integer,           intent(in)    :: iter
    FLOAT,             intent(in)    :: delta
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: diag(:, :)
    CMPLX,             intent(in)    :: offdiag(:, :)
    CMPLX,             intent(inout) :: sp_coeffs(1:length, 0:iter)
    CMPLX,             intent(in)    :: mem_s(intf%np_intf, intf%np_intf, 2)
    integer,           intent(in)    :: length
    integer,           intent(in)    :: order
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: mapping(:)   ! the mapping
    FLOAT,             intent(in)    :: spacing

    integer            :: ip, ii, ji, ki, np
    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :), inv_offdiag(:, :)
    CMPLX, allocatable :: prefactor_plus(:, :), prefactor_minus(:, :)
    CMPLX, allocatable :: sp_tmp(:)
    FLOAT              :: old_norm, norm, sp2

    PUSH_SUB(calculate_sp_coeffs)
    ASSERT(intf%np_intf.eq.intf%np_uc)

    np  = intf%np_intf
    sp2 = spacing**2

    SAFE_ALLOCATE(           tmp1(1:np, 1:np))
    SAFE_ALLOCATE(           tmp2(1:np, 1:np))
    SAFE_ALLOCATE(           tmp3(1:np, 1:np))
    SAFE_ALLOCATE(    inv_offdiag(1:np, 1:np))
    SAFE_ALLOCATE( prefactor_plus(1:np, 1:np))
    SAFE_ALLOCATE(prefactor_minus(1:np, 1:np))
    SAFE_ALLOCATE(sp_tmp(1:length))

    prefactor_plus  = M_zI*delta*diag(:, :)
    prefactor_minus = -M_zI*delta*diag(:, :)

    do ip = 1, np
      prefactor_plus(ip, ip)  = M_ONE + prefactor_plus(ip, ip)
      prefactor_minus(ip, ip) = M_ONE + prefactor_minus(ip, ip)
    end do
    inv_offdiag(:, :) = offdiag(:, :)
    if (mod(intf%il+1,2)+1 .eq. 1) then
      call lalg_invert_upper_triangular(np, inv_offdiag)
    else
      call lalg_invert_lower_triangular(np, inv_offdiag)
    end if

    call lalg_sym_inverter('U', np, prefactor_plus)
    call matrix_symmetrize(prefactor_plus, np)
    if (mod(intf%il+1,2)+1.eq.1) then
      call lalg_trmm(np, np, 'U', 'N', 'L', M_z1, offdiag, prefactor_plus)
      call lalg_trmm(np, np, 'U', 'N', 'R', M_z1, inv_offdiag, prefactor_minus)
    else
      call lalg_trmm(np, np, 'L', 'N', 'L', M_z1, offdiag, prefactor_plus)
      call lalg_trmm(np, np, 'L', 'N', 'R', M_z1, inv_offdiag, prefactor_minus)
    end if

    do ii = start_iter, iter
      sp_coeffs(:, ii) = M_z0
      old_norm = M_ZERO

      do ji = 1, mem_iter
        tmp2(:, :) = M_z0
        ! k = 0 to k = m.
        do ki = 0, ii
          tmp1(:, :) = M_z0
          sp_tmp = sp_coeffs(:, ki)
          if(ki.gt.0) then
            sp_tmp = sp_tmp + M_TWO*sp_coeffs(:, ki-1)
          end if
          if(ki.gt.1) then
            sp_tmp = sp_tmp + sp_coeffs(:, ki-2)
          end if
          call make_full_matrix(np, order, dim, sp_tmp, mem_s, tmp1, mapping)
          if (intf%il.eq.LEFT) then
            call lalg_trmm(np, np, 'U', 'N', 'R', M_z1, inv_offdiag, tmp1)
          else
            call lalg_trmm(np, np, 'L', 'N', 'R', M_z1, inv_offdiag, tmp1)
          end if
          call make_full_matrix(np, order, dim, sp_coeffs(:, ii-ki), mem_s, tmp3, mapping)
          call lalg_symm(np, np, 'R', M_z1, tmp3, tmp1, M_z1, tmp2)
        end do
        call make_full_matrix(np, order, dim, sp_coeffs(:, ii-1), mem_s, tmp3, mapping)
        call lalg_symm(np, np, 'R', M_z1, tmp3, prefactor_minus, TOCMPLX(-delta**2, M_ZERO), tmp2)
        call lalg_gemm(np, np, np, M_z1, prefactor_plus, tmp2, M_z0, tmp1)
        ! Use for numerical stability.
        call matrix_symmetrize(tmp1, np)
        call make_sparse_matrix(np, order, dim, tmp1, mem_s, sp_coeffs(:, ii), mapping)
        norm = infinity_norm(tmp1)
        if((abs(norm-old_norm)*sp2) .lt. mem_tolerance) then
          exit
        end if
        old_norm = norm
      end do
      ! Write a warning if a coefficient is not converged.
      if(ji .gt. mem_iter) then
        write(message(1), '(a,i6,a)') 'Memory coefficent for time step ', ii, &
          ', '//trim(lead_name(intf%il))//' lead, not converged.'
        call messages_warning(1)
      end if

      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ii+1, iter+1)
    end do

    message(1) = ''
    call messages_info(1)

    SAFE_DEALLOCATE_A(tmp1)
    SAFE_DEALLOCATE_A(tmp2)
    SAFE_DEALLOCATE_A(tmp3)
    SAFE_DEALLOCATE_A(inv_offdiag)
    SAFE_DEALLOCATE_A(prefactor_plus)
    SAFE_DEALLOCATE_A(prefactor_minus)
    SAFE_DEALLOCATE_A(sp_tmp)

    POP_SUB(calculate_sp_coeffs)
  end subroutine calculate_sp_coeffs


  ! ---------------------------------------------------------
  !> Write memory coefficients to file.
  !! FIXME: this routine writes compiler- and/or system-dependent binary files.
  !! It should be changed to a platform-independent format.
  subroutine write_coeffs(dir, ob, hm, intf, dim, iter, spacing, delta, op_n, order)
    character(len=*),    intent(in)    :: dir
    type(ob_terms_t),    intent(inout) :: ob
    type(hamiltonian_t), intent(in)    :: hm
    type(interface_t),   intent(in)    :: intf
    integer,             intent(in)    :: dim             !< Dimension.
    integer,             intent(in)    :: iter            !< Number of coefficients.
    FLOAT,               intent(in)    :: spacing         !< Grid spacing.
    FLOAT,               intent(in)    :: delta           !< Timestep.
    integer,             intent(in)    :: op_n            !< Number of operator points.
    integer,             intent(in)    :: order           !< Discretization order.

    integer :: ntime, ip, iunit, np

    PUSH_SUB(write_coeffs)

    call io_mkdir(dir, is_tmp=.true.)
    iunit = io_open(trim(dir)//trim(lead_name(intf%il)), &
                    action='write', form='unformatted', is_tmp=.true.)
    if(iunit.lt.0) then
      message(1) = 'Cannot write memory coefficients to file.'
      call messages_warning(1)
      call io_close(iunit)
      POP_SUB(write_coeffs)
      return
    end if

    np = intf%np_intf

    ! Write numerical parameters.
    write(iunit) dim
    write(iunit) iter
    write(iunit) np
    write(iunit) spacing
    write(iunit) delta
    write(iunit) op_n
    write(iunit) ob%mem_type
    write(iunit) intf%reducible
    write(iunit) hm%lead(intf%il)%vks(1:np, 1)

    ! Write matrices.
    select case(ob%mem_type)
    case(SAVE_CPU_TIME)
      do ntime = 0, iter
        do ip = 1, np
          write(iunit) ob%lead(intf%il)%q(ip, ip:np, ntime)
        end do
      end do
    case(SAVE_RAM_USAGE) ! FIXME: only 2D.
      ASSERT(dim.eq.2)
      write(iunit) ob%lead(intf%il)%q_s(:, :, 1)
      do ntime = 0, iter
        write(iunit) ob%lead(intf%il)%q_sp(1:np*order, ntime)
      end do
    end select

    call io_close(iunit)

    POP_SUB(write_coeffs)
  end subroutine write_coeffs


  ! ---------------------------------------------------------
  !> Read memory coefficients from file.
  !! FIXME: this routine reads compiler- and/or system-dependent binary files.
  !! It should be chanegd to a platform-independent format.
  subroutine read_coeffs(dir, s_iter, ob, hm, intf, dim, iter, spacing, delta, op_n, order)
    character(len=*),    intent(in)    :: dir
    integer,             intent(out)   :: s_iter       !< Number of saved coefficients.
    type(ob_terms_t),    intent(inout) :: ob
    type(hamiltonian_t), intent(in)    :: hm
    type(interface_t),   intent(in)    :: intf
    integer,             intent(in)    :: dim          !< Dimension of the problem.
    integer,             intent(in)    :: iter         !< Number of coefficients.
    FLOAT,               intent(in)    :: spacing      !< Spacing.
    FLOAT,               intent(in)    :: delta        !< Timestep.
    integer,             intent(in)    :: op_n         !< Number of operator points.
    integer,             intent(in)    :: order        !< Discretization order.

    integer :: ntime, ip, iunit, s_dim, s_np, s_op_n, s_mem_type, np, il
    FLOAT   :: s_spacing, s_delta, det
    FLOAT, allocatable :: s_vks(:)
    logical :: s_reducible

    PUSH_SUB(read_coeffs)

    s_iter = 0
    np = intf%np_intf
    il = intf%il

    ! Try to open file.
    iunit = io_open(trim(dir)//trim(lead_name(il)), action='read', &
      status='old', die=.false., is_tmp=.true., form='unformatted')
    if(iunit .lt. 0) then ! no file found
      POP_SUB(read_coeffs)
      return
    end if

    SAFE_ALLOCATE(s_vks(1:np))

    ! Now read the data.
    read(iunit) s_dim
    read(iunit) s_iter
    read(iunit) s_np
    read(iunit) s_spacing
    read(iunit) s_delta
    read(iunit) s_op_n
    read(iunit) s_mem_type
    read(iunit) s_reducible

    ! Check if numerical parameters of saved coefficients match
    ! current parameter set.
    if((s_dim.eq.dim) .and. (s_np.eq.np) .and. (s_op_n.eq.op_n) &
      .and. (s_spacing.eq.spacing) .and. (s_delta.eq.delta) .and. (s_mem_type.eq.ob%mem_type) &
      .and. (s_reducible.eqv.intf%reducible) ) then
      ! read the potential
      read(iunit) s_vks(1:np)
      if(hm%lead(il)%vks(1:np, 1).app.s_vks(1:np)) then
        ! Read the coefficients.
        if (ob%mem_type.eq.SAVE_CPU_TIME) then ! Full (upper half) matrices.
          do ntime = 0, min(iter, s_iter)
            do ip = 1, np
              read(iunit) ob%lead(il)%q(ip, ip:np, ntime)
              ob%lead(il)%q(ip:np, ip, ntime) = ob%lead(il)%q(ip, ip:np, ntime)
            end do
          end do
        else ! Packed matrices (FIXME: yet only 2D).
          ASSERT(dim.eq.2)
          read(iunit) ob%lead(il)%q_s(:, :, 1)
          ob%lead(il)%q_s(:, :, 2) = ob%lead(il)%q_s(1:np, 1:np, 1)
          det = lalg_inverter(np, ob%lead(il)%q_s(:, :, 2), invert=.true.)
          do ntime = 0, min(iter, s_iter)
            read(iunit) ob%lead(il)%q_sp(1:np*order, ntime)
          end do
        end if
      else
        s_iter = 0
      end if
    else
      s_iter = 0
    end if

    call io_close(iunit)

    SAFE_DEALLOCATE_A(s_vks)

    POP_SUB(read_coeffs)
  end subroutine read_coeffs


  ! ---------------------------------------------------------
  !> Calculate amount of (machine) memory required for
  !! the memory term (in MBytes).
  !! Does not consider spin at the moment.
  FLOAT function mbytes_memory_term(iter, np, st, mem_type, order)
    integer,        intent(in) :: iter
    integer,        intent(in) :: np(:)
    type(states_t), intent(in) :: st
    integer,        intent(in) :: mem_type
    integer,        intent(in) :: order

    integer :: il

    PUSH_SUB(mbytes_memory_term)

    mbytes_memory_term = M_ZERO
    do il = 1, NLEADS
      select case(mem_type)
      case(SAVE_CPU_TIME) ! Lots of memory, but faster.
        mbytes_memory_term = mbytes_memory_term + iter*M_TWO*REAL_PRECISION* &
          (st%d%nik*np(il)+np(il)**2) / M_TWO**20
      case(SAVE_RAM_USAGE) ! Less memory, but slower.
        mbytes_memory_term = mbytes_memory_term + iter*M_TWO*REAL_PRECISION* &
          (st%d%nik*np(il)+np(il)*order) / M_TWO**20
      end select
    end do

    POP_SUB(mbytes_memory_term)
  end function mbytes_memory_term


  ! ---------------------------------------------------------
  !> Create sparse matrix out of the full matrix.
  !! Only 2D case.
  !! sp = ss^(-1).mm.ss
  subroutine make_sparse_matrix(np, order, dim, mm, ss, sp, mapping)
    integer, intent(in)  :: np            !< Matrix size.
    integer, intent(in)  :: order         !< Derivative order.
    integer, intent(in)  :: dim           !< Dimension
    CMPLX,   intent(in)  :: mm(np, np)    !< Full matrix.
    CMPLX,   intent(in)  :: ss(np, np, 2) !< Diagonalization matrices.
    CMPLX,   intent(out) :: sp(:)         !< Sparse matrix.
    integer, intent(in)  :: mapping(:)    !< Mapping.

    integer            :: id, ip
    CMPLX, allocatable :: tmp(:, :), tmp2(:, :)

    PUSH_SUB(make_sparse_matrix)

    ASSERT(dim .eq. 2)

    SAFE_ALLOCATE( tmp(1:np, 1:np))
    SAFE_ALLOCATE(tmp2(1:np, 1:np))

    ! Now calculate [S^(-1)*sp*S] to get the sparse matrix.
    ! FIXME: s, inv_s and sp are symmetric.
    call lalg_gemm(np, np, np, M_z1, ss(:, :, 2), mm, M_z0, tmp)
    call lalg_gemm(np, np, np, M_z1, tmp, ss(:, :, 1), M_z0, tmp2)

    select case(dim)
    case(2)
      ! For each row.
      do ip = 0, np-1
        do id = 1, order
          sp(ip*order+id) = tmp2(ip+1, mapping(ip*order+id))
        end do
      end do
    case default
      message(1) = 'Sparse matrices for memory coefficients not supported'
      message(2) = 'in 1D or 3D.'
      call messages_fatal(2)
    end select

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp2)

    POP_SUB(make_sparse_matrix)
  end subroutine make_sparse_matrix


  ! ---------------------------------------------------------
  !> Fast multiplication of the sparse matrix with a dense matrix
  !! res = sp*mm.
  subroutine mul_sp(np, order, sp, mm, mapping, res)
    integer, intent(in)  :: np           !< matrix size (np**2)
    integer, intent(in)  :: order        !< derivative order
    CMPLX,   intent(in)  :: sp(:)        !< sparse matrix (np*order)
    CMPLX,   intent(in)  :: mm(np, np)   !< the matrix to multipliy
    integer, intent(in)  :: mapping(:)   !< the mapping
    CMPLX,   intent(out) :: res(np, np)  ! <the full matrix

    integer :: ii,ij,ik, ijo
    CMPLX   :: tmp

    PUSH_SUB(mul_sp)

    do ij = 1, np
      do ii = 1, np
        tmp = M_z0
        ijo = (ij-1)*order
        do ik = 1, order
          tmp = tmp + sp(ijo+ik)*mm(mapping(ijo+ik),ii)
        end do
        res(ij,ii) = tmp
      end do
    end do

    POP_SUB(mul_sp)
  end subroutine mul_sp


  ! ---------------------------------------------------------
  !> Create the original matrix out of the sparse matrix.
  !! mm = ss.sp.ss^(-1)
  subroutine make_full_matrix(np, order, dim, sp, ss, mm, mapping)
    integer, intent(in)  :: np            !< matrix size (np**2)
    integer, intent(in)  :: order         !< derivative order
    integer, intent(in)  :: dim           !< dimension
    CMPLX,   intent(in)  :: sp(:)         !< sparse matrix (np*order)
    CMPLX,   intent(in)  :: ss(np, np, 2) !< diagonalization matrix
    CMPLX,   intent(out) :: mm(np, np)    !< the full matrix
    integer, intent(in)  :: mapping(:)    !< the mapping

    CMPLX, allocatable :: tmp(:, :)
    integer            :: id, ip

    PUSH_SUB(make_full_matrix)

    SAFE_ALLOCATE(tmp(1:np, 1:np))

    ! make a full matrix of the packed storage form
    ! the sparse matrix has a (symmetric) banded form
    mm = M_z0
    select case(dim)
    case(2)
      do ip = 0, np-1
        do id = 1, order
          mm(ip+1,mapping(ip*order+id)) = sp(ip*order+id)
        end do
      end do
    case(3)
      ! FIXME (not yet done)
    end select
    ! now calculate [ss*sp*ss^(-1)] to get the dense matrix
    call mul_sp(np, order, sp, ss(:,:,2), mapping, tmp)
!    call lalg_gemm(np, np, np, M_z1, mm, ss(:,:,2), M_z0, tmp)
    call lalg_gemm(np, np, np, M_z1, ss(:,:,1), tmp, M_z0, mm)
    call matrix_symmetric_average(mm, np)

    SAFE_DEALLOCATE_A(tmp)
    POP_SUB(make_full_matrix)
  end subroutine make_full_matrix


  ! ---------------------------------------------------------
  !> Free arrays.
  subroutine ob_mem_end(ob)
    type(ob_terms_t), intent(inout) :: ob
      
    integer :: il

    PUSH_SUB(ob_mem_end)

    do il=1, NLEADS
      SAFE_DEALLOCATE_P(ob%lead(il)%q)
      SAFE_DEALLOCATE_P(ob%lead(il)%q_sp)
      SAFE_DEALLOCATE_P(ob%lead(il)%sp2full_map)
      SAFE_DEALLOCATE_P(ob%lead(il)%q_s)
    end do

    POP_SUB(ob_mem_end)
  end subroutine ob_mem_end
end module ob_mem_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
