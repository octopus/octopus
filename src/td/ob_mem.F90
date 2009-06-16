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
    type(interface_t),   intent(in)    :: intface(1:NLEADS)
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
      SAFE_ALLOCATE(ob%mem_coeff(1:np, 1:np, 0:max_iter, 1:NLEADS))
      nullify(ob%mem_sp_coeff, ob%mem_s, ob%sp2full_map)
      ob%mem_coeff = M_z0
    case(SAVE_RAM_USAGE) ! FIXME: only 2D.
      ASSERT(calc_dim.eq.2)
      ! sp_coeff has the size ny*nz*do**2, (np=ny*nz*do).
      SAFE_ALLOCATE(ob%mem_sp_coeff(1:np*order, 0:max_iter, 1:NLEADS))
      SAFE_ALLOCATE(ob%mem_s(1:np, 1:np, 1:2, 1:NLEADS))
      SAFE_ALLOCATE(ob%sp2full_map(1:np*order))
      nullify(ob%mem_coeff)
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
          call calculate_coeffs_ni(il, saved_iter+1, max_iter, delta, intface(il), hm%lead_h_diag(:, :, 1, il), &
            hm%lead_h_offdiag(:, :, il), ob%mem_coeff(:, :, :, il))
        case(SAVE_RAM_USAGE) ! FIXME: only 2D.
          ASSERT(calc_dim.eq.2)
          call apply_coupling(ob%mem_coeff(:, :, 0, il), hm%lead_h_offdiag(:, :, il),&
              ob%mem_coeff(:, :, 0, il), np, il)
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

      if(ob%mem_type.eq.SAVE_CPU_TIME) then
        do ii=0, max_iter
          call apply_coupling(ob%mem_coeff(:, :, ii, il), hm%lead_h_offdiag(:, :, il),&
              ob%mem_coeff(:, :, ii, il), np, il)
        end do
      end if

    end do

    call pop_sub()
  end subroutine ob_mem_init


  ! ---------------------------------------------------------
  ! Solve for zeroth memory coefficient by truncating the continued
  ! matrix fraction. Since the coefficient must be symmetric (non-magnetic),
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
    FLOAT              :: norm, old_norm, d2
    CMPLX              :: h

    call push_sub('ob_mem.approx_coeff0_new')

    np = intface%np
    d2 = delta**2

    ! If we are in 1D and have only a number we can solve the equation explicitly.
    ! So check this first for faster calculation.
    if(np.eq.1) then
      h            = M_z1 + M_zI*delta*diag(1,1)
      coeff0(1, 1) = (-h + sqrt(h**2 + d2*(M_TWO*offdiag(1,1))**2)) / (M_TWO*d2*offdiag(1,1)**2)
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
      do i = 1, mem_iter
        ! Calculate 1 + i*delta*h + delta^2*coeff0_old
        coeff0(:, :) = M_zI*delta*diag(:, :) + d2*q0(:, :)
        do j = 1, np
          coeff0(j, j) = M_ONE + coeff0(j, j)
        end do

        ! Invert.
        call lalg_sym_inverter('U', np, coeff0)
        call matrix_symmetrize(coeff0, np)
        norm = infinity_norm(coeff0)
        if((abs(old_norm/norm-M_ONE)).lt.mem_tolerance) then
          exit
        end if
        ! Apply coupling matrices.
        call apply_coupling(coeff0, offdiag, q0, np, il)
        call matrix_symmetric_average(q0, np)
        old_norm = norm
      end do
      if(i.gt.mem_iter) then
        write(message(1), '(a,i6,a)') 'Memory coefficent for time step 0, ' &
          //trim(lead_name(il))//' lead, not converged'
        call write_warning(1)
      end if

      if(mem_type.eq.SAVE_RAM_USAGE) then
        ! Diagonalization procedure.
        mem_s(:, :, 1) = q0(:, :)
        call lalg_eigensolve_nonh(np, mem_s(:, :, 1), mem_s(:, 1, 2))
        mem_s(:, :, 2) = mem_s(:, :, 1)
        norm = lalg_inverter(np, mem_s(:, :, 2), invert=.true.)
        call make_sparse_matrix(np, order, 2, q0, mem_s, sp_coeff0, mapping)
      end if
      SAFE_DEALLOCATE_A(q0)
    end if

    call pop_sub()
  end subroutine approx_coeff0


  ! ---------------------------------------------------------
  ! New version of calculating he memory coefficients.
  ! It now calculates and stores the p-coefficients.
  ! This method is more general and can even continue if
  ! the matrices are not invertable. It uses only the memory
  ! which is needed by the q-matrices (defined by q=v^T.p.v).
  ! These coefficients have to be calculated after this routine
  ! is finished and the coefficients are stored on the disk.
  ! This version cannot handle magnetic fields in the leads, 
  ! as the assumption of symmetric matrices for multiplications
  ! is not valid anymore.
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
    CMPLX, allocatable :: tmp(:, :), tmp2(:, :), m0(:, :), m_l(:, :), m_r(:, :)
    FLOAT              :: norm, old_norm

    call push_sub('ob_mem.calculate_coeffs_ni_new')

    np = intf%np

    SAFE_ALLOCATE( tmp(1:np, 1:np))
    SAFE_ALLOCATE(tmp2(1:np, 1:np))
    SAFE_ALLOCATE(  m0(1:np, 1:np))
    SAFE_ALLOCATE( m_l(1:np, 1:np))
    SAFE_ALLOCATE( m_r(1:np, 1:np))

    m0 = M_z0
    tmp(:, :) = -M_zI*delta*diag(:, :)
    do i = 1, np
      tmp(i, i) = M_ONE + tmp(i, i)
    end do

    call lalg_symm(np, np, 'L', M_z1, coeffs(:, :, 0), tmp, M_z0, m0)

    m_l(:, :) = delta*coeffs(:, :, 0)
    m_r(:, :) = m_l(:, :)
    if(il.eq.LEFT) then
      call lalg_trmm(np, np, 'U', 'N', 'R', M_z1, offdiag, m_l)
      call lalg_trmm(np, np, 'U', 'T', 'L', M_z1, offdiag, m_r)
    else
      call lalg_trmm(np, np, 'L', 'N', 'R', M_z1, offdiag, m_l)
      call lalg_trmm(np, np, 'L', 'T', 'L', M_z1, offdiag, m_r)
    end if

    do i = start_iter, iter
      ! The part of the sum independent of coeff_p(:, :, i),
      ! accumulated in tmp2.
      tmp2 = M_z0
      if(i.gt.1) then
        tmp(:, :) = M_TWO*coeffs(:, :, i-1) + coeffs(:, :, i-2)
      else
        tmp(:, :) = M_TWO*coeffs(:, :, i-1)
      end if
      call lalg_symm(np, np, 'L', M_z1, tmp, m_r, M_z0, tmp2)
      
      do k = 1, i-1       
        if(k.gt.1) then
          tmp(:, :) = coeffs(:, :, k) + M_TWO*coeffs(:, :, k-1) + coeffs(:, :, k-2)
        else
          tmp(:, :) = coeffs(:, :, k) + M_TWO*coeffs(:, :, k-1)
        end if
        if(il.eq.LEFT) then
          call lalg_trmm(np, np, 'U', 'T', 'R', M_z1, offdiag, tmp)
        else
          call lalg_trmm(np, np, 'L', 'T', 'R', M_z1, offdiag, tmp)
        end if
        call lalg_symm(np, np, 'R', TOCMPLX(delta, M_ZERO), coeffs(:, :, i-k), tmp, M_z1, tmp2)
      end do

      call lalg_symm(np, np, 'R', M_z1, coeffs(:, :, i-1), m0, M_z0, tmp)
      call lalg_gemm(np, np, np, -M_z1, m_l, tmp2, M_z1, tmp)
      call matrix_symmetric_average(tmp, np)

      
      ! initialize coefficient
      coeffs(:, :, i) = tmp(:, :)
      old_norm = infinity_norm(coeffs(:, :, i))
      ! The convergence loop.
      do j = 1, mem_iter
        call lalg_symm(np, np, 'L', M_z1, coeffs(:, :, i), m_r, M_z0, tmp2)
!        call lalg_gemm(np, np, np, M_z1, coeffs(:, :, i), m_r, M_z0, tmp2)
        call lalg_gemm(np, np, np, M_z1, m_l, tmp2, M_z0, coeffs(:, :, i))
        coeffs(:, :, i) = tmp(:, :) - coeffs(:, :, i)
        call matrix_symmetric_average(coeffs(:, :, i), np)

        norm = infinity_norm(coeffs(:, :, i))
        if(abs(old_norm/norm-M_ONE).lt.mem_tolerance) then
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

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp2)
    SAFE_DEALLOCATE_A(m0)
    SAFE_DEALLOCATE_A(m_l)
    SAFE_DEALLOCATE_A(m_r)

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

    SAFE_ALLOCATE(           tmp1(1:np, 1:np))
    SAFE_ALLOCATE(           tmp2(1:np, 1:np))
    SAFE_ALLOCATE(           tmp3(1:np, 1:np))
    SAFE_ALLOCATE(    inv_offdiag(1:np, 1:np))
    SAFE_ALLOCATE( prefactor_plus(1:np, 1:np))
    SAFE_ALLOCATE(prefactor_minus(1:np, 1:np))
    SAFE_ALLOCATE(sp_tmp(1:length))

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

    SAFE_DEALLOCATE_A(tmp1)
    SAFE_DEALLOCATE_A(tmp2)
    SAFE_DEALLOCATE_A(tmp3)
    SAFE_DEALLOCATE_A(inv_offdiag)
    SAFE_DEALLOCATE_A(prefactor_plus)
    SAFE_DEALLOCATE_A(prefactor_minus)
    SAFE_DEALLOCATE_A(sp_tmp)

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
    if(iunit.lt.0) then ! no file found
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
    integer,        intent(in) :: nl
    integer,        intent(in) :: np(nl)
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

    SAFE_ALLOCATE( tmp(1:np, 1:np))
    SAFE_ALLOCATE(tmp2(1:np, 1:np))

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

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp2)

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

    SAFE_ALLOCATE(tmp(1:np, 1:np))

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

    SAFE_DEALLOCATE_A(tmp)
    call pop_sub()
  end subroutine make_full_matrix


  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine ob_mem_end(ob)
    type(ob_terms_t), intent(inout) :: ob

    call push_sub('ob_mem.ob_mem_end')

    SAFE_DEALLOCATE_P(ob%mem_coeff)
    SAFE_DEALLOCATE_P(ob%mem_sp_coeff)
    SAFE_DEALLOCATE_P(ob%sp2full_map)
    SAFE_DEALLOCATE_P(ob%mem_s)

    call pop_sub()
  end subroutine ob_mem_end
end module ob_mem_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
