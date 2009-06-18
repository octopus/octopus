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

module ob_lead_m
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
    apply_coupling, &
    lead_diag,      &
    lead_offdiag,   &
    lead_td_pot

contains

  ! ---------------------------------------------------------
  ! Calculate the diagonal block matrix, i. e. the Hamiltonian for
  ! entries contained in the interface region.
  subroutine lead_diag(lapl, vks, intf, diag)
    type(nl_operator_t), intent(in)  :: lapl
    FLOAT,               intent(in)  :: vks(:)
    type(interface_t),   intent(in)  :: intf
    CMPLX,               intent(out) :: diag(:, :)

    integer :: i, k, n, k_stencil, intf_index
    FLOAT   :: w_re, w_im

    call push_sub('ob_lead.lead_diag')

    diag = M_z0

    ! For all interface points...
    do i = 1, intf%np
      n = intf%index(i)
      ! ... find the points they are coupled to.
      do k = 1, lapl%stencil%size
        k_stencil = lapl%m%idx%Lxyz_inv(          &
          lapl%m%idx%Lxyz(n,1)+lapl%stencil%points(1,k), &
          lapl%m%idx%Lxyz(n,2)+lapl%stencil%points(2,k), &
          lapl%m%idx%Lxyz(n,3)+lapl%stencil%points(3,k))
        ! If the coupling point is in the interface...
        if(k_stencil.le.lapl%np.and. &
          member_of_interface(k_stencil, intf, intf_index)) then
          ! ... get the operator coefficients.
          if(lapl%cmplx_op) then
            if(lapl%const_w) then
              w_re = lapl%w_re(k, 1)
              w_im = lapl%w_im(k, 1)
            else
              w_re = lapl%w_re(k, n)
              w_im = lapl%w_im(k, n)
            end if
          else
            if(lapl%const_w) then
              w_re = lapl%w_re(k, 1)
            else
              w_re = lapl%w_re(k, n)
            end if
            w_im = M_ZERO
          end if
          ! Calculation if the kinetic term: -1/2 prefactor.
          w_im = -w_im/M_TWO
          w_re = -w_re/M_TWO
          ! ... write them into the right entry of the diagonal block.
          diag(i, intf_index) = TOCMPLX(w_re, w_im)
        end if
      end do
    end do

    ! Add potential.
    do i = 1, intf%np
      diag(i, i) = diag(i, i) + vks(i)
    end do
    call pop_sub()
  end subroutine lead_diag


  ! ---------------------------------------------------------
  ! Calculate the offdiagonal block matrix, i. e. the Hamiltonian for
  ! entries not contained in the interface region, which is the matrix
  ! V^T_\alpha or H_{C\alpha}.
  subroutine lead_offdiag(lapl, intf, il, offdiag)
    type(nl_operator_t), intent(in)  :: lapl
    type(interface_t),   intent(in)  :: intf
    integer,             intent(in)  :: il
    CMPLX,               intent(out) :: offdiag(:, :)

    integer :: p_n(MAX_DIM), p_k(MAX_DIM), p_matr(MAX_DIM)
    integer :: j, k, k_stencil
    integer :: n, n_matr
    integer :: dir
    integer :: x_shift
    FLOAT   :: w_re, w_im

    call push_sub('ob_lead.lead_offdiag')

    ! Coupling direction.
    select case(il)
    case(LEFT)
      dir = 1
    case(RIGHT)
      dir = -1
    end select

    offdiag(:, :) = M_z0

    ! j iterates over rows of the block matrix.
    do j = 1, intf%np
      n = intf%index(j)
      ! k iterates over all stencil points.
      do k = 1, lapl%stencil%size
        ! Get point number of coupling point.
        k_stencil = lapl%m%idx%Lxyz_inv(            &
          lapl%m%idx%Lxyz(n, 1)+lapl%stencil%points(1, k), &
          lapl%m%idx%Lxyz(n, 2)+lapl%stencil%points(2, k), &
          lapl%m%idx%Lxyz(n, 3)+lapl%stencil%points(3, k))

        ! Get coordinates of current interface point n and current stencil point k_stencil.
        p_n = lapl%m%idx%Lxyz(n, :)
        p_k = lapl%m%idx%Lxyz(k_stencil, :)

        ! Now, we shift the stencil by the size of one unit cell (intf%extent)
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
        p_matr(TRANS_DIR) = p_matr(TRANS_DIR) + dir*(intf%extent-x_shift)
        n_matr            = lapl%m%idx%Lxyz_inv(p_matr(1), p_matr(2), p_matr(3))

        if(member_of_interface(n_matr, intf, n_matr)) then

          ! Multiply by the coefficient of the operator and sum up.
          if(lapl%cmplx_op) then
            if(lapl%const_w) then
              w_re = lapl%w_re(k, 1)
              w_im = lapl%w_im(k, 1)
            else
              w_re = lapl%w_re(k, n)
              w_im = lapl%w_im(k, n)
            end if
          else
            if(lapl%const_w) then
              w_re = lapl%w_re(k, 1)
            else
              w_re = lapl%w_re(k, n)
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
  end subroutine lead_offdiag


  ! ---------------------------------------------------------
  ! Calculate res <- offdiag^T matrix offdiag with all matrices np x np.
  ! If matrix is symmetric, so is the result.
  subroutine apply_coupling(matrix, offdiag, res, np, il)
    CMPLX,   intent(in)    :: matrix(np, np)
    CMPLX,   intent(in)    :: offdiag(np, np)
    integer, intent(in)    :: np, il
    CMPLX,   intent(out)   :: res(np, np)

    call push_sub('ob_lead.apply_coupling')

    res(1:np, 1:np) = matrix(1:np, 1:np)
    if(il.eq.LEFT) then
      call lalg_trmm(np, np, 'U', 'N', 'L', M_z1, offdiag, res)
      call lalg_trmm(np, np, 'U', 'T', 'R', M_z1, offdiag, res)
    else
      call lalg_trmm(np, np, 'L', 'N', 'L', M_z1, offdiag, res)
      call lalg_trmm(np, np, 'L', 'T', 'R', M_z1, offdiag, res)
    end if

    call pop_sub()
  end subroutine apply_coupling


  ! ---------------------------------------------------------
  ! Calculates the time-dependent lead potential from formula string.
  subroutine lead_td_pot(td_pot, formula, n_steps, tstep)
    integer,          intent(in)  :: n_steps
    FLOAT,            intent(out) :: td_pot(0:n_steps+1, 1:NLEADS)
    character(len=*), intent(in)  :: formula(:)
    FLOAT,            intent(in)  :: tstep

    integer             :: it, il
    FLOAT               :: t
    FLOAT, allocatable  :: pot_im(:)
    character(len=1024) :: tmp_c_string
    
    FLOAT :: zero(1)

    call push_sub('ob_lead.lead_td_pot')

    zero(1) = M_ZERO

    SAFE_ALLOCATE(pot_im(0:n_steps + 1))

    ! Calculate td potential.
    do il = 1, NLEADS
      tmp_c_string = formula(il)
      call conv_to_c_string(tmp_c_string)

      td_pot(0, il) = M_ZERO
      do it = 1, n_steps + 1
        t = it*tstep
        call loct_parse_expression(td_pot(it, il), pot_im(it), 1, zero, zero(1), t, tmp_c_string)
      end do
    end do

    SAFE_DEALLOCATE_A(pot_im)
    call pop_sub()
  end subroutine lead_td_pot

end module ob_lead_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
