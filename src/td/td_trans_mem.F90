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

! Calculation of the memory coefficients for the modified Crank-Nicholson propagator.

#include "global.h"

module td_trans_mem_m
  use datasets_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_parser_m
  use loct_m
  use math_m
  use messages_m
  use mpi_m
  use nl_operator_m
  use states_m
  use system_m
  use td_trans_intf_m
  use varinfo_m
  use io_m

  implicit none

  private
  public ::                 &
    mbytes_memory_term,     &
    memory_init,            &
    memory_end,             &
    make_full_matrix,       &
    make_symmetric,         &
    make_symmetric_average, &
    apply_coupling
    

  integer :: mem_iter
  FLOAT   :: mem_tolerance

contains

  ! ---------------------------------------------------------
  ! Calculate amount of (machine) memory required for
  ! the memory term (in MBytes).
  ! Does not consider spin at the moment.
  FLOAT function mbytes_memory_term(iter, np, nl, st, mem_type, order)
    integer,        intent(in) :: iter
    integer,        intent(in) :: np
    integer,        intent(in) :: nl
    type(states_t), intent(in) :: st
    integer,        intent(in) :: mem_type
    integer,        intent(in) :: order

    call push_sub('td_transport.mbytes_memory_term')
    mbytes_memory_term = M_ZERO
    select case(mem_type)
    case(1) ! lots of memory, but fast
      mbytes_memory_term = nl*iter*M_TWO*REAL_PRECISION* &
        (st%nst*np+np**2) / M_TWO**20
    case(2) ! less memory, but slower
      mbytes_memory_term = nl*iter*M_TWO*REAL_PRECISION* &
        (st%nst*np+np*order) / M_TWO**20
    end select

    call pop_sub()
  end function mbytes_memory_term

  ! ---------------------------------------------------------
  ! create sparse matrix out of the full matrix
  ! only 2D case yet
  ! sp = s^(-1).m.s
  subroutine make_sparse_matrix(np, order, dim, m, s, sp, mapping)
    integer,        intent(in)  :: np          ! matrix size (np**2)
    integer,        intent(in)  :: order       ! derivative order
    integer,        intent(in)  :: dim         ! dimension
    CMPLX,          intent(in)  :: m(np, np)   ! full matrix
    CMPLX,          intent(in)  :: s(np, np, 2)! diagonalization matrices
    CMPLX,          intent(out) :: sp(:)       ! the sparse matrix
    integer,        intent(in)  :: mapping(:)  ! the mapping

    CMPLX, allocatable :: tmp(:, :), tmp2(:, :)
    integer :: id, i
    call push_sub('td_transport.make_sparse_matrix')

    ALLOCATE(tmp(np,np),np**2)
    ALLOCATE(tmp2(np,np),np**2)

    ! now calculate [S^(-1)*sp*S] to get the sparse matrix
    ! FIXME: s, inv_s and sp are symmetric
    call lalg_gemm(np, np, np, M_z1, s(:,:,2), m, M_z0, tmp)
    call lalg_gemm(np, np, np, M_z1, tmp, s(:,:,1), M_z0, tmp2)

    select case(dim)
    case(2)
      ! for each row
      do i=0, np-1
        do id=1, order
          sp(i*order+id) = tmp2(i+1,mapping(i*order+id))
        end do
      end do

    case(3)
      ! FIXME (not yet done)
      ! tribanded offdiagonals
    end select
!write(*,*) count, order*np

    deallocate(tmp, tmp2)
    call pop_sub()
  end subroutine make_sparse_matrix

  ! ---------------------------------------------------------
  ! fast multiplication of the sparse matrix with a dense matrix
  ! res = sp*m
  subroutine mul_sp(np, order, dim, sp, m, mapping, res)
    integer,        intent(in)  :: np          ! matrix size (np**2)
    integer,        intent(in)  :: order       ! derivative order
    integer,        intent(in)  :: dim         ! dimension
    CMPLX,          intent(in)  :: sp(:)       ! sparse matrix (np*order)
    CMPLX,          intent(in)  :: m(np, np)   ! the matrix to multipliy
    integer,        intent(in)  :: mapping(:)  ! the mapping
    CMPLX,          intent(out) :: res(np, np) ! the full matrix

    integer :: ii,ij,ik, ijo
    CMPLX   :: tmp

    call push_sub('td_transport.mul_sp')
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
  ! create the original matrix out of the sparse matrix
  ! m = s.sp.s^(-1)
  subroutine make_full_matrix(np, order, dim, sp, s, m, mapping)
    integer,        intent(in)  :: np          ! matrix size (np**2)
    integer,        intent(in)  :: order       ! derivative order
    integer,        intent(in)  :: dim         ! dimension
    CMPLX,          intent(in)  :: sp(:)       ! sparse matrix (np*order)
    CMPLX,          intent(in)  :: s(np, np, 2)! diagonalization matrix
    CMPLX,          intent(out) :: m(np, np)   ! the full matrix
    integer,        intent(in)  :: mapping(:)  ! the mapping

    CMPLX, allocatable :: tmp(:, :)
    integer :: id, i
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
    call make_symmetric_average(m, np)

    deallocate(tmp)
    call pop_sub()
  end subroutine make_full_matrix

  ! ---------------------------------------------------------
  ! Allocate (machine) memory for the memory coefficents and
  ! calculate them.
  subroutine memory_init(intface, delta, max_iter, op, coeff, sp_coeff, mem_s, &
    diag, offdiag, dim, spacing, mem_type, order, sp2full_map, mpi_grp)
    type(intface_t),     intent(in) :: intface
    FLOAT,               intent(in) :: delta
    integer,             intent(in) :: max_iter
    type(nl_operator_t), intent(in) :: op
    CMPLX,               pointer    :: coeff(:, :, :, :)
    CMPLX,               pointer    :: sp_coeff(:, :, :)
    CMPLX,               pointer    :: mem_s(:, :, :, :)
    CMPLX,               pointer    :: diag(:, :, :)
    CMPLX,               pointer    :: offdiag(:, :, :)
    integer,             intent(in) :: dim
    FLOAT,               intent(in) :: spacing
    integer,             intent(in) :: mem_type
    integer,             intent(in) :: order
    integer,             pointer    :: sp2full_map(:)
    type(mpi_grp_t),     intent(in) :: mpi_grp

    integer            :: il, np, saved_iter, ii

    call push_sub('td_trans_mem.memory_coeff')

    np = intface%np
    ALLOCATE(diag(np, np, NLEADS), np**2*NLEADS)
    ALLOCATE(offdiag(np, np, NLEADS), np**2*NLEADS)
    if (mem_type.eq.1) then
      ALLOCATE(coeff(np, np, 0:max_iter, NLEADS), np**2*(max_iter+1)*NLEADS)
      ALLOCATE(sp_coeff(0, 0, 0), 0)
      ALLOCATE(mem_s(0, 0, 0, 0), 0)
      ALLOCATE(sp2full_map(0), 0)
    else ! FIXME: yet only 2D
      ALLOCATE(coeff(0, 0, 0, 0), 0)
      ! sp_coeff has the size ny*nz*do**2, (np=ny*nz*do)
      ALLOCATE(sp_coeff(np*order, 0:max_iter, NLEADS), order*np*(max_iter+1)*NLEADS)
      ALLOCATE(mem_s(np, np, 2, NLEADS), np**2*2*NLEADS)
      ALLOCATE(sp2full_map(np*order), order*np)
      ! fill sparse to full mapping matrix
      do ii = 0, np-1
        do il = 1, order
          sp2full_map(ii*order+il) = mod((il-1)*(np/order) + ii , np) + 1
        end do
      end do
    end if

    coeff = M_z0
    ! Get iteration parameters from input file.

    !%Variable TDTransMemTol
    !%Type float
    !%Default 1e-12
    !%Section Transport
    !%Description
    !% Decides when to consider the memory coefficients converged.
    !%End
    call loct_parse_float(check_inp('TDTransMemTol'), CNST(1e-12), mem_tolerance)
    if(mem_tolerance.le.M_ZERO) then
      write(message(1), '(a,f14.6,a)') "Input : '", mem_tolerance, "' is not a valid TDTransMemTol."
      message(2) = '(0 < TDTransMemTol)'
      call write_fatal(2)
    end if

    !%Variable TDTransMemMaxIter
    !%Type integer
    !%Default 500
    !%Section Transport
    !%Description
    !% Sets the maximum iteration number to converge the memory coefficients.
    !%End
    call loct_parse_int(check_inp('TDTransMemMaxIter'), 500, mem_iter)
    if(mem_iter.lt.0) then
      write(message(1), '(a,i6,a)') "Input : '", mem_iter, "' is not a valid TDTransMemMaxIter."
      message(2) = '(0 <= TDTransMemMaxIter)'
      call write_fatal(2)
    end if
    do il = 1, NLEADS
      ! Get diagonal matrix.
      call calculate_diag(op, intface, il, diag(:, :, il))
      ! Get offdiagonal matrix.
      call calculate_offdiag(op, intface, il, offdiag(:, :, il))
!write(message(1),'(a, i3)') 'Info: Number of threads ', mc%nthreads

      ! Try to read the coefficients from file
      call read_coeffs(io_workpath("transport/"), saved_iter, coeff(:, :, :, il), sp_coeff(:, :, il), mem_s(:,:,:,il), &
                       dim, max_iter, intface%np, spacing, delta, op%n, mem_type, np*order, il)

      if (saved_iter.lt.max_iter) then ! calculate coefficients
        if (saved_iter.gt.0) then
          write(message(1),'(a,i5,a)') 'Info: Successfully loaded the first',saved_iter,' memory coefficients from '// &
                       trim(lead_name(il))//' lead.'
          call write_info(1)
        end if
        message(1) = 'Info: Calculating missing coefficients for memory term of '// &
                      trim(lead_name(il))//' lead.'
        call write_info(1)
        ! initialize progress bar
        call loct_progress_bar(-1, max_iter+1)

        if (saved_iter.eq.0) then 
          ! Get anchor for recursion.
          call approx_coeff0(intface, delta, il, diag(:, :, il), offdiag(:, :, il), coeff(:, :, 0, il), &
                             sp_coeff(:, 0, il), mem_s(:,:,:,il), order, mem_type, sp2full_map, spacing)
          call loct_progress_bar(1, max_iter+1)
        end if
        ! Calculate the subsequent coefficients by the recursive relation.
        if (mem_type.eq.1) then
          call calculate_coeffs(il, saved_iter+1, max_iter, delta, intface, diag(:, :, il), &
                     offdiag(:, :, il), coeff(:, :, :, il), spacing)
        else ! FIXME: yet only 2D
          call calculate_sp_coeffs(il, saved_iter+1, max_iter, delta, intface, diag(:, :, il), &
                     offdiag(:, :, il), sp_coeff(:, :, il), mem_s(:,:,:,il),np*order, &
                      order, dim, sp2full_map, spacing)
        end if

        if (saved_iter.lt.max_iter) then
          message(1) = ''
          message(2) = 'Info: Writing memory coefficients of '// &
            trim(lead_name(il))//' lead.'
          call write_info(2)
          if(mpi_grp_is_root(mpi_grp)) then
            call write_coeffs(io_workpath("transport/"), coeff(:, :, :, il), sp_coeff(:, :, il), mem_s(:,:,:,il), dim, &
              max_iter, intface%np, spacing, delta, op%n, mem_type, np*order, il)
          end if
        end if
      else
        message(1) = 'Info: Successfully loaded memory coefficients from '// &
          trim(lead_name(il))//' lead.'
        call write_info(1)
      end if

    end do

    call pop_sub()
  end subroutine memory_init

  ! ---------------------------------------------------------
  ! copies the upper matrix into the lower part
  subroutine make_symmetric(matrix, np)
    CMPLX,          intent(inout) :: matrix(np, np)
    integer,        intent(in) :: np

    integer   ::  j
    call push_sub('td_transport.make_symmetric')

    do j=1, np-1
      matrix(j+1:np,j) = matrix(j,j+1:np)
    end do

    call pop_sub()
  end subroutine make_symmetric

  ! ---------------------------------------------------------
  ! takes the average of the matrix and its transposed
  subroutine make_symmetric_average(matrix, np)
    CMPLX,          intent(inout) :: matrix(np, np)
    integer,        intent(in) :: np

    integer   ::  j
    CMPLX,  allocatable :: tmp(:,:)
    call push_sub('td_transport.make_symmetric_average')
    ALLOCATE(tmp(np,np),np**2)
   
    if (np.gt.1) then
      do j=1, np
        tmp(:,j) = M_HALF*(matrix(j,:)+matrix(:,j))
      end do
      matrix = tmp
    end if

    deallocate(tmp)
    call pop_sub()
  end subroutine make_symmetric_average

  ! ---------------------------------------------------------
  ! Solve for zeroth memory coefficient by truncating the continued
  ! matrix fraction. Since the coefficient must be symmetric,
  ! a symmetric inversion is used.
  subroutine approx_coeff0(intface, delta, il, diag, offdiag, coeff0, sp_coeff0, &
                            mem_s, order, mem_type, mapping, spacing)
    type(intface_t),     intent(in)  :: intface
    FLOAT,               intent(in)  :: delta
    integer,             intent(in)  :: il
    CMPLX,               intent(in)  :: diag(:, :)
    CMPLX,               intent(in)  :: offdiag(:, :)
    CMPLX,               intent(out) :: coeff0(:, :)
    CMPLX,               intent(out) :: sp_coeff0(:)   ! the 0th coefficient in packed storage
    CMPLX,               intent(out) :: mem_s(:, :, :) ! S & S^(-1)
    integer,             intent(in)  :: order
    integer,             intent(in)  :: mem_type
    integer,             intent(in)  :: mapping(:)   ! the mapping
    FLOAT,               intent(in)  :: spacing

    integer            :: i, j, np
    CMPLX, allocatable :: q0(:, :)
    FLOAT              :: norm, old_norm, sp2
    CMPLX              :: h


    call push_sub('td_trans_mem.approx_coeff0')

    np = intface%np

    ! If we are in 1d mode and have only a number we can solve the equation explicitely.
    ! So check this first for faster calculation.
    if (np.eq.1) then
      h = M_z1 + M_zI*delta*diag(1,1)
      coeff0(1,1) = ( -h + sqrt(h**2 + (M_TWO*delta*offdiag(1,1))**2) ) / (M_TWO*delta**2)
    else ! we have the general case of a matrix, so solve the equation by iteration
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
        call make_symmetric(q0, np)

        ! Apply coupling matrices.
        call apply_coupling(q0, offdiag, q0, np, il)
        call make_symmetric_average(q0, np)
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

      if (mem_type.eq.1) then
        coeff0 = q0

      else
        ! diagonalization procedure
        mem_s(:, :, 1) = q0(:, :)
        call lalg_eigensolve_nonh(np, mem_s(:, :, 1), mem_s(:, 1, 2))
        mem_s(:, :, 2) = mem_s(:, :, 1)
        norm = lalg_inverter(np, mem_s(:, :, 2), invert = .true.)
        call make_sparse_matrix(np, order, 2, q0, mem_s, sp_coeff0, mapping)
      end if
      deallocate(q0)
    end if

    call pop_sub()
  end subroutine approx_coeff0

  ! ---------------------------------------------------------
  ! output of a matrix (for debugging)
  subroutine write_matrix(matr, iter, np)
    CMPLX,     intent(in)  :: matr(np, np)
    integer,   intent(in)  :: iter,np
    integer            :: i, j

    write(*,'(a,i4,a)') "--- Matrix(",iter,") ---"
    do j = 1, np
      write(*,'(i3,a,i3,a,11e20.9)') iter, "m(",j,",i)=",abs(matr(j,:))
    end do
!    do j = 1, np
!      do i = 1, np
!        write(*,'(i3,a,i3,a,i3,a,1e20.9)') iter, "m(",j,",",i,")=",abs(matr(j,i))
!        !write(*,'(i3,a,i3,a,i3,a,2e20.9)') iter, "m(",j,",",i,")=",real(matr(j,i)), aimag(matr(j,i))
!      end do
!    end do
   end subroutine write_matrix


  ! ---------------------------------------------------------
  ! Calculate the diagonal block matrix, i. e. the Hamiltonian for
  ! entries contained in the interface region (for zero potential only
  ! at the moment, thus, only kinetic energy).
  subroutine calculate_diag(op, intf, il, diag)
    type(nl_operator_t), intent(in)  :: op
    type(intface_t),     intent(in)  :: intf
    integer,             intent(in)  :: il
    CMPLX,               intent(out) :: diag(:, :)

    integer :: i, k, n, k_stencil
    FLOAT   :: w_re, w_im

    call push_sub('td_trans_mem.calculate_diag')

    diag = M_z0

    ! For all interface points...
    do i = 1, intf%np
      n = intf%index(i, il)
      ! ... find the points they are coupled to.
      do k = 1, op%n
        k_stencil = op%m%lxyz_inv(op%m%lxyz(n,1)+op%stencil(1,k), &
                                  op%m%lxyz(n,2)+op%stencil(2,k), &
                                  op%m%lxyz(n,3)+op%stencil(3,k))
        ! If the coupling point is in the interface...
        if(k_stencil.le.op%np.and. &
          member_of_intface(k_stencil, intf, il)) then
          ! ... get the operator coefficients.
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
            w_im = M_ZERO
          end if
          ! Calculation if the kinetic term: -1/2 prefactor.
          w_im = -w_im/M_TWO
          w_re = -w_re/M_TWO
          ! ... write them into the right entry of the diagonal block.
          diag(i, k_stencil-intf%index_range(1, il)+1) = TOCMPLX(w_re, w_im)
        end if
      end do
    end do
    call pop_sub()
  end subroutine calculate_diag


  ! ---------------------------------------------------------
  ! Calculate the offdiagonal block matrix, i. e. the Hamiltonian for
  ! entries not contained in the interface region (for zero potential only
  ! at the moment, thus, only kinetic energy), which is the matrix
  ! V^T_{\alpha } or H_{C\alpha}
  subroutine calculate_offdiag(op, intface, il, offdiag)
    type(nl_operator_t), intent(in)  :: op
    type(intface_t),     intent(in)  :: intface
    integer,             intent(in)  :: il
    CMPLX,               intent(out) :: offdiag(:, :)

    integer :: p_n(MAX_DIM), p_k(MAX_DIM), p_matr(MAX_DIM)
    integer :: j, k, k_stencil
    integer :: n, n_matr
    integer :: dir
    integer :: x_shift
    FLOAT   :: w_re, w_im

    call push_sub('td_trans_mem.calculate_offdiag')

    ! Coupling direction.
    select case(il)
    case(LEFT)
      dir = 1
    case(RIGHT)
      dir = -1
    end select

    offdiag(:, :) = M_z0

    ! j iterates over rows of the block matrix.
    do j = 1, intface%np
      n = intface%index(j, il)
      ! k iterates over all stencil points.
      do k = 1, op%n
        ! Get point number of coupling point.
        k_stencil = op%m%lxyz_inv(          &
          op%m%lxyz(n, 1)+op%stencil(1, k), &
          op%m%lxyz(n, 2)+op%stencil(2, k), &
          op%m%lxyz(n, 3)+op%stencil(3, k))

        ! Get coordinates of current interface point n and current stencil point k_stencil.
        p_n = op%m%lxyz(n, :)
        p_k = op%m%lxyz(k_stencil, :)

        ! Now, we shift the stencil by the size of one unit cell (intface%extent)
        ! and check if the coupling point with point number n_matr is in the interface.
        ! The sketch if for a 4x2 unit cell and a 5 point stencil in 2D:
        !
        !    |       ||                     |       ||
        ! ...+---+---++---+---+...       ...+---+---++---+---+...
        !    | o |   ||   |   |             |   |   || o |   |
        ! ...+-|-+---++---+---+...       ...+---+---++-|-+---+...
        !  x---o---o ||   |   |             |   | #----o---o |
        ! ...+-|-+---++---+---+...  ==>  ...+---+---++-|-+---+...
        !    | o |   ||   |   |             |   |   || o |   |
        ! ...+---+---++---+---+...       ...+---+---++---+---+...
        !    |   |   ||   |   |             |   |   ||   |   |
        ! ...+---+---++---+---+...       ...+---+---++---+---+...
        !    |   L   ||   C                 |   L   ||   C
        !      unit
        !      cell
        !
        ! The point p_k is markes with an x, the point p_matr with a #.
        x_shift           = abs(p_k(TRANS_DIR)-p_n(TRANS_DIR))
        p_matr            = p_n
        p_matr(TRANS_DIR) = p_matr(TRANS_DIR) + dir*(intface%extent-x_shift)
        n_matr            = op%m%lxyz_inv(p_matr(1), p_matr(2), p_matr(3))

        if(member_of_intface(n_matr, intface, il)) then
          n_matr = intface_index(n_matr, intface, il)

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
            w_im = M_ZERO
          end if

          ! Calculation of the kinetic term: -1/2 prefactor.
          w_im = -w_im/M_TWO
          w_re = -w_re/M_TWO
          offdiag(j, n_matr) = TOCMPLX(w_re, w_im)
        end if
      end do
    end do

    call pop_sub()
  end subroutine calculate_offdiag


  ! ---------------------------------------------------------
  ! Calculate res <- offdiag^T matrix offdiag. with all matrices np x np.
  ! If matrix is symmetric, so is the result.
  subroutine apply_coupling(matrix, offdiag, res, np, il)
    CMPLX,               intent(inout)  :: matrix(np, np)
    CMPLX,               intent(in)  :: offdiag(np, np)
    integer,             intent(in)  :: np, il
    CMPLX,               intent(out) :: res(np, np)

    call push_sub('td_trans_mem.apply_coupling')

    res(:, :) = matrix(:, :)
    if (il.eq.LEFT) then
      call lalg_trmm(np,np,'U','N','l',M_z1,offdiag,res)
      call lalg_trmm(np,np,'U','T','r',M_z1,offdiag,res)
    else
      call lalg_trmm(np,np,'L','N','l',M_z1,offdiag,res)
      call lalg_trmm(np,np,'L','T','r',M_z1,offdiag,res)
    end if

    call pop_sub()
  end subroutine apply_coupling

  ! ---------------------------------------------------------
  ! coeffs(:, :, :, 0) given, calculate the subsequent ones by
  ! the recursive relation. Since the 0th coefficiant is symmetric
  ! all subsequent are also, therefore symmetric matrix multiplications
  ! can be used. use recursive relation without p (uses only half the memory)
  subroutine calculate_coeffs(il, start_iter, iter, delta, intf, diag, offdiag, coeffs, spacing)
    integer,             intent(in)    :: il
    integer,             intent(in)    :: start_iter
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: delta
    type(intface_t),     intent(in)    :: intf
    CMPLX,               intent(in)    :: diag(intf%np, intf%np)
    CMPLX,               intent(in)    :: offdiag(intf%np, intf%np)
    CMPLX,               intent(inout) :: coeffs(intf%np, intf%np, 0:iter)
    FLOAT,               intent(in)    :: spacing

    FLOAT              :: old_norm, norm, sp2, d2
    integer            :: i,j, k, np
    CMPLX, allocatable :: tmp(:, :), tmp2(:, :), inv_offdiag(:, :)
    CMPLX, allocatable :: prefactor_plus(:, :), prefactor_minus(:, :)
    CMPLX              :: tp

    call push_sub('td_trans_mem.calculate_coeffs')

    np = intf%np
    sp2 = spacing**2
    d2 = delta**2

    ALLOCATE(tmp(np, np), np**2)
    ALLOCATE(tmp2(np, np), np**2)
    ALLOCATE(inv_offdiag(np, np), np**2)
    ALLOCATE(prefactor_plus(np, np), np**2)
    ALLOCATE(prefactor_minus(np, np), np**2)

    prefactor_plus(:, :)  = M_zI*delta*diag(:, :)
    prefactor_minus(:, :) = -M_zI*delta*diag(:, :)

    do i = 1, np
      prefactor_plus(i, i)  = M_one + prefactor_plus(i, i)
      prefactor_minus(i, i) = M_one + prefactor_minus(i, i)
    end do

    ! If we are in 1d mode and have only a number we can calculate the coefficients explicitely.
    ! So check this first for faster calculation.
    if (np.eq.1) then
      tp = M_ONE/(prefactor_plus(1, 1) + M_two*d2*coeffs(1, 1, 0))
      do i = start_iter, iter       ! i=start_iter to m.
        if (i.eq.1) then
          coeffs(1, 1, 1) = coeffs(1,1,0)*(prefactor_minus(1, 1) - M_two*d2*coeffs(1, 1, 0))*tp
        else ! i>=2
          !j = 1
          coeffs(1, 1, i) = (coeffs(1,1,1)+M_two*coeffs(1,1,0))*coeffs(1,1,i-1) + coeffs(1,1,0)*coeffs(1,1,i-2)
          do j = 2, i-1
            coeffs(1,1,i) = coeffs(1,1,i) + (coeffs(1,1,j)+M_two*coeffs(1,1,j-1)+coeffs(1,1,j-2))*coeffs(1,1,i-j)
          end do
          coeffs(1,1,i) = coeffs(1,1,1)*coeffs(1,1,i-1)/coeffs(1,1,0) - d2*tp*coeffs(1,1,i)
        end if
        call loct_progress_bar(i+1, iter+1)
      end do
    else ! we have the general case of a matrix, so solve the equation by iteration
      inv_offdiag(:, :) = offdiag(:, :)
      if (il.eq.LEFT) then
        call lalg_invert_upper_triangular(np, inv_offdiag)
      else
        call lalg_invert_lower_triangular(np, inv_offdiag)
      end if

      call lalg_sym_inverter('U', np, prefactor_plus)
      call make_symmetric(prefactor_plus, np)
      if (il.eq.LEFT) then
        call lalg_trmm(np,np,'U','N','L',M_z1,offdiag,prefactor_plus)
        call lalg_trmm(np,np,'U','N','R',M_z1,inv_offdiag,prefactor_minus)
      else
        call lalg_trmm(np,np,'L','N','L',M_z1,offdiag,prefactor_plus)
        call lalg_trmm(np,np,'L','N','R',M_z1,inv_offdiag,prefactor_minus)
      end if
      do i = start_iter, iter       ! i=start_iter to m.
        coeffs(:, :, i) = M_z0
        old_norm = M_ZERO
        do j = 1, mem_iter ! maxiter to converge matrix equation for q(i)

          tmp2(:,:) = M_z0

          ! k = 0 to k = m.
          do k = 0, i
            tmp(:, :) = coeffs(:, :, k)
            if(k.gt.0) then
              tmp(:, :) = tmp(:, :) + M_TWO*coeffs(:, :, k-1)
            end if
            if(k.gt.1) then
              tmp(:, :) = tmp(:, :) + coeffs(:, :, k-2)
            end if
            if (il.eq.LEFT) then
              call lalg_trmm(np,np,'U','N','R',M_z1,inv_offdiag,tmp)
            else
              call lalg_trmm(np,np,'L','N','R',M_z1,inv_offdiag,tmp)
            end if
            call zsymm('R','U',np,np,M_z1,coeffs(:, :, i-k),np,tmp,np,M_z1,tmp2,np)
          end do
          call zsymm('R','U',np,np,M_z1,coeffs(:, :, i-1),np,prefactor_minus,np,TOCMPLX(-delta**2, M_ZERO),tmp2,np)
          call zgemm('N','N',np,np,np,M_z1,prefactor_plus,np,tmp2,np,M_z0,coeffs(:, :, i),np)
          ! use for numerical stability
          call make_symmetric_average(coeffs(:, :, i), np)

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
    deallocate(prefactor_plus, prefactor_minus)

    call pop_sub()
  end subroutine calculate_coeffs

  ! ---------------------------------------------------------
  ! a bit faster, but uses twice the memory
  ! coeffs(:, :, :, 0) given, calculate the subsequent ones by
  ! the recursive relation. Since the 0th coefficiant is symmetric
  ! all subsequent are also, therefore symmetric matrix multiplications
  ! can be used.
  subroutine calculate_coeffs_old(il, start_iter, iter, delta, intf, diag, offdiag, coeffs)
    integer,             intent(in)    :: il
    integer,             intent(in)    :: start_iter
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: delta
    type(intface_t),     intent(in)    :: intf
    CMPLX,               intent(in)    :: diag(:, :)
    CMPLX,               intent(in)    :: offdiag(:, :)
    CMPLX,               intent(inout) :: coeffs(intf%np, intf%np, 0:iter)

    integer            :: i,j, k, np
    CMPLX, allocatable :: coeff_p(:, :, :), p_prev(:, :), tmp(:, :)
    CMPLX, allocatable :: prefactor_plus(:, :), prefactor_minus(:, :)
    FLOAT              :: norm, old_norm

    call push_sub('td_trans_mem.calculate_coeffs')

    np = intf%np

    ALLOCATE(coeff_p(np, np, 0:iter), np**2*(iter+1))
    ALLOCATE(p_prev(np, np), np**2)
    ALLOCATE(tmp(np, np), np**2)
    ALLOCATE(prefactor_plus(np, np), np**2)
    ALLOCATE(prefactor_minus(np, np), np**2)

    coeff_p = M_z0
    prefactor_plus  = M_zI*delta*diag(:, :)
    prefactor_minus = -M_zI*delta*diag(:, :)

    do i = 1, np
      prefactor_plus(i, i)  = M_one + prefactor_plus(i, i)
      prefactor_minus(i, i) = M_one + prefactor_minus(i, i)
    end do

    ! Calculate p_\alpha = 1/(1 + i*delta*h_\alpha + delta^2*q_\alpha)
    coeff_p(:, :, 0) = prefactor_plus + delta**2 * coeffs(:, :, 0)

    call lalg_sym_inverter('U', np, coeff_p(:, :, 0))
    call make_symmetric(coeff_p(:, :, 0), np)
    call lalg_sym_inverter('U', np, prefactor_plus)
    call make_symmetric(prefactor_plus, np)

    tmp = offdiag
    if (il.eq.LEFT) then
      call lalg_invert_upper_triangular(np, tmp)
    else
      call lalg_invert_lower_triangular(np, tmp)
    end if

    do i = 1, start_iter-1
      call apply_coupling(coeffs(:, :, i), tmp, coeff_p(:, :, i), np, il)
    end do

    do i = start_iter, iter       ! i=1 to m.
      old_norm = M_ZERO
      do j = 1, mem_iter ! maxiter to converge matrix equation for p(i) and q(i)

        p_prev = coeff_p(:, :, i)

        ! k = m.
        call apply_coupling(p_prev, offdiag, coeff_p(:, :, i), np, il)
        coeff_p(:, :, i) = coeff_p(:, :, i) + M_TWO*coeffs(:, :, i-1)

        if(i.gt.1) then
          coeff_p(:, :, i) = coeff_p(:, :, i) + coeffs(:, :, i-2)
        end if

        call zsymm('L','U',np,np,M_z1,coeff_p(:, :, i),np,coeff_p(:, :, 0),np,M_z0,tmp,np)
        coeff_p(:, :, i) = tmp

        ! k = 0.
        call zsymm('L','U',np,np,M_z1,coeffs(:, :, 0),np,p_prev,np,M_z1,coeff_p(:, :, i),np)

        ! k = 1 to k = m-1.
        do k = 1, i-1
          tmp = coeffs(:, :, k) + M_TWO*coeffs(:, :, k-1)
          if(k.gt.1) then
            tmp = tmp + coeffs(:, :, k-2)
          end if
          call zsymm('L','U',np,np,M_z1,tmp,np,coeff_p(:, :, i-k),np,M_z1,coeff_p(:, :, i),np)
        end do

        coeff_p(:, :, i) = -delta**2 * coeff_p(:, :, i)
        call zsymm('L','U',np,np,M_z1,prefactor_minus,np,coeff_p(:, :, i-1),np,M_z1,coeff_p(:, :, i),np)
        call zsymm('L','U',np,np,M_z1,prefactor_plus,np,coeff_p(:, :, i),np,M_z0,tmp,np)
        coeff_p(:, :, i) = tmp
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

    deallocate(coeff_p, p_prev, tmp)
    deallocate(prefactor_plus, prefactor_minus)

    call pop_sub()
  end subroutine calculate_coeffs_old

  ! ---------------------------------------------------------
  ! same as calculate_coeffs but with the sparse matrix
  ! sp_coeffs(:, :, 0) given, calculate the subsequent ones by
  ! the recursive relation. We can only use the (sparse) mem_q, so use the
  ! mem_q only recursive relation.
  subroutine calculate_sp_coeffs(il, start_iter, iter, delta, intf, diag, offdiag, &
                                 sp_coeffs, mem_s, length, order, dim, mapping, spacing)
    integer,             intent(in)    :: il
    integer,             intent(in)    :: start_iter
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: delta
    type(intface_t),     intent(in)    :: intf
    CMPLX,               intent(in)    :: diag(:, :)
    CMPLX,               intent(in)    :: offdiag(:, :)
    CMPLX,               intent(inout) :: sp_coeffs(1:length, 0:iter)
    CMPLX,               intent(in)    :: mem_s(intf%np, intf%np, 2)
    integer,             intent(in)    :: length
    integer,             intent(in)    :: order
    integer,             intent(in)    :: dim
    integer,             intent(in)    :: mapping(:)   ! the mapping
    FLOAT,               intent(in)    :: spacing

    integer            :: i,j, k, np
    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :), inv_offdiag(:, :)
    CMPLX, allocatable :: prefactor_plus(:, :), prefactor_minus(:, :)
    CMPLX, allocatable :: sp_tmp(:)
    FLOAT              :: old_norm, norm, sp2

    call push_sub('td_trans_mem.calculate_sp_coeffs')
    np = intf%np
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
      prefactor_plus(i, i)  = M_one + prefactor_plus(i, i)
      prefactor_minus(i, i) = M_one + prefactor_minus(i, i)
    end do
    inv_offdiag(:, :) = offdiag(:, :)
    if (il.eq.LEFT) then
      call lalg_invert_upper_triangular(np, inv_offdiag)
    else
      call lalg_invert_lower_triangular(np, inv_offdiag)
    end if

    call lalg_sym_inverter('U', np, prefactor_plus)
    call make_symmetric(prefactor_plus, np)
    if (il.eq.LEFT) then
      call lalg_trmm(np,np,'U','N','L',M_z1,offdiag,prefactor_plus)
      call lalg_trmm(np,np,'U','N','R',M_z1,inv_offdiag,prefactor_minus)
    else
      call lalg_trmm(np,np,'L','N','L',M_z1,offdiag,prefactor_plus)
      call lalg_trmm(np,np,'L','N','R',M_z1,inv_offdiag,prefactor_minus)
    end if

    do i = start_iter, iter       ! i=start_iter to m.
      sp_coeffs(:, i) = M_z0
      old_norm = M_ZERO

      do j = 1, mem_iter ! maxiter to converge matrix equation for q(i)
        tmp2(:,:) = M_z0
        ! k = 0 to k = m.
        do k = 0, i
          tmp1(:,:) = M_z0
          sp_tmp = sp_coeffs(:, k)
          if(k.gt.0) then
            sp_tmp = sp_tmp + M_TWO*sp_coeffs(:, k-1)
          end if
          if(k.gt.1) then
            sp_tmp = sp_tmp + sp_coeffs(:, k-2)
          end if
          call make_full_matrix(np, order, dim, sp_tmp, mem_s, tmp1, mapping)
          if (il.eq.LEFT) then
            call lalg_trmm(np,np,'U','N','R',M_z1,inv_offdiag,tmp1)
          else
            call lalg_trmm(np,np,'L','N','R',M_z1,inv_offdiag,tmp1)
          end if
          call make_full_matrix(np, order, dim, sp_coeffs(:, i-k), mem_s, tmp3, mapping)
          call zsymm('R','U',np,np,M_z1,tmp3,np,tmp1,np,M_z1,tmp2,np)
        end do
        call make_full_matrix(np, order, dim, sp_coeffs(:, i-1), mem_s, tmp3, mapping)
        call zsymm('R','U',np,np,M_z1,tmp3,np,prefactor_minus,np,TOCMPLX(-delta**2, M_ZERO),tmp2,np)
        call zgemm('N','N',np,np,np,M_z1,prefactor_plus,np,tmp2,np,M_z0,tmp1,np)
        ! use for numerical stability
        call make_symmetric_average(tmp1, np)
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
  subroutine write_coeffs(dir, coeffs, sp_coeffs, mem_s, dim, iter, np, spacing, delta, op_n, mem_type, length, il)
    character(len=*), intent(in) :: dir
    CMPLX,     intent(in)        :: coeffs(np, np, 0:iter) ! the saved coefficients
    CMPLX,     intent(in)        :: sp_coeffs(length, 0:iter) ! the saved coefficients
    CMPLX,     intent(in)        :: mem_s(np, np, 2)! the diagonalization matrix s
    integer,   intent(in)        :: dim          ! dimension
    integer,   intent(in)        :: iter         ! the number of coefficients
    integer,   intent(in)        :: np           ! number of points in the interface
    FLOAT,     intent(in)        :: spacing      ! grid spacing
    FLOAT,     intent(in)        :: delta        ! timestep
    integer,   intent(in)        :: op_n         ! number of operator points
    integer,   intent(in)        :: mem_type     ! 1: normal  2: packed
    integer,   intent(in)        :: length       ! length of the packed array
    integer,   intent(in)        :: il           ! which lead

    integer     :: ntime, j, iunit

    call push_sub('td_trans_mem.write_coeffs')

    call io_mkdir(dir, is_tmp=.true.)
    iunit = io_open(trim(dir)//trim(lead_name(il)), action='write', form='unformatted')
    if(iunit < 0) then
      write(*,*) "Error, can't write to file!"
      call io_close(iunit)
      return
    end if

    write(iunit) dim
    write(iunit) iter
    write(iunit) np
    write(iunit) spacing
    write(iunit) delta
    write(iunit) op_n
    write(iunit) mem_type

    if (mem_type.eq.1) then
      do ntime=0, iter
        do j=1, np
          write(iunit) coeffs(j,j:np,ntime)
        end do
      end do
    else ! FIXME: yet only 2D
      write(iunit) mem_s(:,:,1)
      do ntime=0, iter
        write(iunit) sp_coeffs(1:length,ntime)
      end do
    end if

    call io_close(iunit)

    call pop_sub()
  end subroutine write_coeffs

  ! ---------------------------------------------------------
  ! Read memory coefficients from file.
  subroutine read_coeffs(dir, s_iter, coeffs, sp_coeffs, mem_s, dim, iter, np, spacing, delta, op_n, mem_type, length, il)
    character(len=*), intent(in)    :: dir
    integer,   intent(out)   :: s_iter   ! the number of saved coefficients
    CMPLX,     intent(inout) :: coeffs(np, np, 0:iter) ! the saved coefficients
    CMPLX,     intent(inout) :: sp_coeffs(length, 0:iter) ! the saved coefficients
    CMPLX,     intent(out)   :: mem_s(np, np, 2) ! the diagonalization matrices
    integer,   intent(in)    :: dim          ! dimension of the problem
    integer,   intent(in)    :: iter         ! the number of coefficients
    integer,   intent(in)    :: np           ! number of points in the interface
    FLOAT,     intent(in)    :: spacing      ! spacing 
    FLOAT,     intent(in)    :: delta        ! timestep
    integer,   intent(in)    :: op_n         ! number of operator points
    integer,   intent(in)    :: mem_type     ! which lead
    integer,   intent(in)    :: length       ! length of the packed array
    integer,   intent(in)    :: il           ! which lead

    integer     :: ntime, j, iunit, s_dim, s_np, s_op_n, s_mem_type
    FLOAT       :: s_spacing, s_delta, det

    call push_sub('td_trans_mem.read_coeffs')
    s_iter = 0
    ! try to read from file
    iunit = io_open(trim(dir)//trim(lead_name(il)), action='read', &
                    status='old', die=.false., is_tmp = .true., form='unformatted')
    if(iunit < 0) then
      call io_close(iunit)
      return
    end if
    ! now read the data
    read(iunit) s_dim
    read(iunit) s_iter
    read(iunit) s_np
    read(iunit) s_spacing
    read(iunit) s_delta
    read(iunit) s_op_n
    read(iunit) s_mem_type

    if ((s_dim.eq.dim) .and. (s_np.eq.np) .and. (s_op_n.eq.op_n) &
          .and. (s_spacing.eq.spacing) .and. (s_delta.eq.delta) &
          .and. (s_mem_type.eq.mem_type)) then
      ! read the coefficients
      if (mem_type.eq.1) then ! full (upper half) matrices
        do ntime=0, min(iter,s_iter)
          do j=1, np
            read(iunit) coeffs(j,j:np,ntime)
            coeffs(j:np,j,ntime) = coeffs(j,j:np,ntime)
          end do
        end do
      else ! packed matrices ! FIXME: yet only 2D
        read(iunit) mem_s(:,:,1)
        mem_s(:, :, 2) = mem_s(1:np, 1:np, 1)
        det = lalg_inverter(np, mem_s(:, :, 2), invert = .true.)
        do ntime=0, min(iter,s_iter)
          read(iunit) sp_coeffs(1:length,ntime)
        end do
      end if

    else
      s_iter = 0
    end if

    call io_close(iunit)

    call pop_sub()
  end subroutine read_coeffs

  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine memory_end(coeff, sp_coeff, diag, offdiag, mem_s, sp2full_map)
    CMPLX, pointer   :: coeff(:, :, :, :)
    CMPLX, pointer   :: sp_coeff(:, :, :)
    CMPLX, pointer   :: diag(:, :, :)
    CMPLX, pointer   :: offdiag(:, :, :)
    CMPLX, pointer   :: mem_s(:, :, :, :)
    integer, pointer :: sp2full_map(:)

    call push_sub('td_trans_mem.memory_end')

    if(associated(coeff)) then
      deallocate(coeff)
      nullify(coeff)
    end if

    if(associated(sp_coeff)) then
      deallocate(sp_coeff)
      nullify(sp_coeff)
    end if

    if(associated(diag)) then
      deallocate(diag)
      nullify(diag)
    end if

    if(associated(offdiag)) then
      deallocate(offdiag)
      nullify(offdiag)
    end if

    if(associated(mem_s)) then
      deallocate(mem_s)
      nullify(mem_s)
    end if

    if(associated(sp2full_map)) then
      deallocate(sp2full_map)
      nullify(sp2full_map)
    end if

    call pop_sub()
  end subroutine memory_end
end module td_trans_mem_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
