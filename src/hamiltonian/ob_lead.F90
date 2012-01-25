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

module ob_lead_m
  use global_m
  use grid_m
  use lalg_adv_m
  use lalg_basic_m
  use parser_m
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
    lead_init_pot,  &
    lead_init_kin,  &
    lead_end,       &
    lead_diag,      &
    lead_offdiag,   &
    lead_td_pot,    &
    lead_resize,    &
    is_lead_transl_inv

contains

  ! ---------------------------------------------------------
  !> Calculate the diagonal block matrix, i.e. the Hamiltonian for
  !! entries contained in the interface region.
  subroutine lead_diag(lapl, vks, intf, diag)
    type(nl_operator_t), intent(in)  :: lapl
    FLOAT,               intent(in)  :: vks(:)
    type(interface_t),   intent(in)  :: intf
    CMPLX,               intent(out) :: diag(:, :)

    integer :: ip, kk, nn, k_stencil, intf_index
    FLOAT   :: w_re, w_im

    PUSH_SUB(lead_diag)

    diag = M_z0

    ! For all interface points...
    do ip = 1, intf%np_uc
      nn = intf%index(ip)
      ! ... find the points they are coupled to.
      do kk = 1, lapl%stencil%size
        k_stencil = lapl%mesh%idx%lxyz_inv(          &
          lapl%mesh%idx%lxyz(nn,1)+lapl%stencil%points(1,kk), &
          lapl%mesh%idx%lxyz(nn,2)+lapl%stencil%points(2,kk), &
          lapl%mesh%idx%lxyz(nn,3)+lapl%stencil%points(3,kk))
        ! If the coupling point is in the interface...
        if(k_stencil.le.lapl%np.and. &
          member_of_interface(k_stencil, intf, intf_index)) then
          ! ... get the operator coefficients.
          if(lapl%cmplx_op) then
            if(lapl%const_w) then
              w_re = lapl%w_re(kk, 1)
              w_im = lapl%w_im(kk, 1)
            else
              w_re = lapl%w_re(kk, nn)
              w_im = lapl%w_im(kk, nn)
            end if
          else
            if(lapl%const_w) then
              w_re = lapl%w_re(kk, 1)
            else
              w_re = lapl%w_re(kk, nn)
            end if
            w_im = M_ZERO
          end if
          ! Calculation if the kinetic term: -1/2 prefactor.
          w_im = -w_im/M_TWO
          w_re = -w_re/M_TWO
          ! ... write them into the right entry of the diagonal block.
          diag(ip, intf_index) = TOCMPLX(w_re, w_im)
        end if
      end do
    end do

    ! Add potential. FIXME: add vector potential (A^2)
    forall(ip = 1:intf%np_uc) diag(ip, ip) = diag(ip, ip) + vks(ip)

    POP_SUB(lead_diag)
  end subroutine lead_diag


  ! ---------------------------------------------------------
  !> Calculate the off-diagonal block matrix, i.e. the Hamiltonian for
  !! entries not contained in the interface region, which is the matrix
  !! \f$V^T_\alpha\f$ or \f$H_{C\alpha}\f$.
  subroutine lead_offdiag(lapl, intf, offdiag)
    type(nl_operator_t), intent(in)  :: lapl
    type(interface_t),   intent(in)  :: intf
    CMPLX,               intent(out) :: offdiag(:, :)

    integer :: p_n(MAX_DIM), p_k(MAX_DIM), p_matr(MAX_DIM)
    integer :: jp, kk, k_stencil, nn, n_matr, dir, tdir, shift, intf_idx
    FLOAT   :: w_re, w_im

    PUSH_SUB(lead_offdiag)

    ! Coupling direction.
    dir = (-1)**(intf%il+1)
    tdir = (intf%il+1)/2

    offdiag(:, :) = M_z0

    ! jp iterates over rows of the block matrix.
    do jp = 1, intf%np_uc
      nn = intf%index(jp)
      ! k iterates over all stencil points.
      do kk = 1, lapl%stencil%size
        ! Get point number of coupling point.
        k_stencil = lapl%mesh%idx%lxyz_inv(            &
          lapl%mesh%idx%lxyz(nn, 1)+lapl%stencil%points(1, kk), &
          lapl%mesh%idx%lxyz(nn, 2)+lapl%stencil%points(2, kk), &
          lapl%mesh%idx%lxyz(nn, 3)+lapl%stencil%points(3, kk))

        ! Get coordinates of current interface point nn and current stencil point k_stencil.
        p_n = lapl%mesh%idx%lxyz(nn, :)
        p_k = lapl%mesh%idx%lxyz(k_stencil, :)

        ! Now, we shift the stencil by the size of one unit cell (intf%extent_uc)
        ! and check if the coupling point with point number n_matr is in the interface.
        ! The sketch is for a 4x2 unit cell and a 5-point stencil in 2D:
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
        ! The point p_k is marked with an x, the point p_matr with a #.
        shift        = abs(p_k(tdir)-p_n(tdir))
        p_matr       = p_n
        p_matr(tdir) = p_matr(tdir) + dir*(intf%extent_uc-shift)
        n_matr       = lapl%mesh%idx%lxyz_inv(p_matr(1), p_matr(2), p_matr(3))

        if(member_of_interface(n_matr, intf, intf_idx)) then

          ! Multiply by the coefficient of the operator and sum up.
          if(lapl%cmplx_op) then
            if(lapl%const_w) then
              w_re = lapl%w_re(kk, 1)
              w_im = lapl%w_im(kk, 1)
            else
              w_re = lapl%w_re(kk, nn)
              w_im = lapl%w_im(kk, nn)
            end if
          else
            if(lapl%const_w) then
              w_re = lapl%w_re(kk, 1)
            else
              w_re = lapl%w_re(kk, nn)
            end if
            w_im = M_ZERO
          end if

          ! Calculation of the kinetic term: -1/2 prefactor.
          w_im = -w_im/M_TWO
          w_re = -w_re/M_TWO
          offdiag(jp, intf_idx) = TOCMPLX(w_re, w_im)
        end if
      end do
    end do

    POP_SUB(lead_offdiag)
  end subroutine lead_offdiag


  ! ---------------------------------------------------------
  !> Calculate res <- offdiag^T matrix offdiag with all matrices np x np.
  !! If matrix is symmetric, so is the result.
  subroutine apply_coupling(matrix, offdiag, res, np, il)
    CMPLX,   intent(in)    :: matrix(np, np)
    CMPLX,   intent(in)    :: offdiag(np, np)
    integer, intent(in)    :: np, il
    CMPLX,   intent(out)   :: res(np, np)

    PUSH_SUB(apply_coupling)

    res(1:np, 1:np) = matrix(1:np, 1:np)
    if(mod(il+1,2)+1.eq.1) then
      call lalg_trmm(np, np, 'U', 'N', 'L', M_z1, offdiag, res)
      call lalg_trmm(np, np, 'U', 'T', 'R', M_z1, offdiag, res)
    else
      call lalg_trmm(np, np, 'L', 'N', 'L', M_z1, offdiag, res)
      call lalg_trmm(np, np, 'L', 'T', 'R', M_z1, offdiag, res)
    end if

    POP_SUB(apply_coupling)
  end subroutine apply_coupling


  ! ---------------------------------------------------------
  !> Calculates the time-dependent lead potential from formula string.
  subroutine lead_td_pot(td_pot, formula, n_steps, tstep)
    integer,          intent(in)  :: n_steps
    FLOAT,            intent(out) :: td_pot(0:n_steps+1, 1:NLEADS)
    character(len=*), intent(in)  :: formula(:)
    FLOAT,            intent(in)  :: tstep

    integer             :: it, il
    FLOAT               :: tt
    FLOAT, allocatable  :: pot_im(:)
    character(len=1024) :: tmp_c_string

    FLOAT :: zero(1)

    PUSH_SUB(lead_td_pot)

    zero(1) = M_ZERO

    SAFE_ALLOCATE(pot_im(0:n_steps + 1))

    ! Calculate td potential.
    do il = 1, NLEADS
      tmp_c_string = formula(il)
      call conv_to_c_string(tmp_c_string)

      td_pot(0, il) = M_ZERO
      do it = 1, n_steps + 1
      ! FIXME: if a better algorithm for the gs is implemented this should change.
      ! Then the initial bias also matters
!      do it = 0, n_steps + 1
        tt = it*tstep
        call parse_expression(td_pot(it, il), pot_im(it), 1, zero, zero(1), tt, tmp_c_string)
      end do
    end do

    SAFE_DEALLOCATE_A(pot_im)
    POP_SUB(lead_td_pot)
  end subroutine lead_td_pot


  !> init the potentials of the lead
  subroutine lead_init_pot(lead, np, np_part, nspin)
    type(lead_t),        intent(out) :: lead
    integer,             intent(in)  :: np
    integer,             intent(in)  :: np_part ! including ghost and boundary points
    integer,             intent(in)  :: nspin

    PUSH_SUB(lead_init_pot)

    lead%np  = np
    lead%np_part = np_part
    SAFE_ALLOCATE(lead%vks(1:np, 1:nspin))
    SAFE_ALLOCATE(lead%vh(1:np))
    SAFE_ALLOCATE(lead%v0(1:np))

    POP_SUB(lead_init_pot)
  end subroutine lead_init_pot

  !> init the kinetic part (diag and offdiag) of the lead
  subroutine lead_init_kin(lead, np, np_part, dim)
    type(lead_t),        intent(out) :: lead
    integer,             intent(in)  :: np
    integer,             intent(in)  :: np_part ! including ghost and boundary points
    integer,             intent(in)  :: dim

    PUSH_SUB(lead_init_kin)

    lead%np  = np
    lead%np_part = np_part
    SAFE_ALLOCATE(lead%h_diag(1:np, 1:np, 1:dim))
    SAFE_ALLOCATE(lead%h_offdiag(1:np, 1:np))

    POP_SUB(lead_init_kin)
  end subroutine lead_init_kin


  ! ---------------------------------------------------------
  subroutine lead_end(lead)
    type(lead_t),        intent(inout)  :: lead

    PUSH_SUB(lead_end)

    SAFE_DEALLOCATE_P(lead%h_diag)
    SAFE_DEALLOCATE_P(lead%h_offdiag)
    SAFE_DEALLOCATE_P(lead%vks)
    SAFE_DEALLOCATE_P(lead%vh)
    SAFE_DEALLOCATE_P(lead%v0)

    POP_SUB(lead_end)
  end subroutine lead_end


  ! ---------------------------------------------------------
  subroutine lead_copy(src, dst, dim, nspin)
    type(lead_t),  intent(in)    :: src
    type(lead_t),  intent(inout) :: dst
    integer,       intent(in)    :: dim
    integer,       intent(in)    :: nspin

    integer :: np, ip

    PUSH_SUB(lead_copy)

    ASSERT(dst%np.le.src%np)
    np = dst%np

    forall(ip = 1:np) dst%h_diag(1:np, ip, 1:dim) = src%h_diag(1:np, ip, 1:dim)
    forall(ip = 1:np) dst%h_offdiag(1:np, ip)     = src%h_offdiag(1:np, ip)
    dst%vks(1:np, 1:nspin) = src%vks(1:np, 1:nspin)
    dst%vh(1:np)           = src%vh(1:np)
    dst%v0(1:np)           = src%v0(1:np)

    POP_SUB(lead_copy)
  end subroutine lead_copy


  ! ---------------------------------------------------------
  !> Resizes the lead unit cell according to the interface size.
  !!
  !! \todo if hamiltonian->lead is an allocatable array then pointers
  !! are the better choice instead of copying
  subroutine lead_resize(intf, lead, dim, nspin)
    type(interface_t),   intent(in)    :: intf
    type(lead_t),        intent(inout) :: lead
    integer,             intent(in)    :: dim
    integer,             intent(in)    :: nspin

    type(lead_t) :: old_lead
    integer :: np, np_part

    PUSH_SUB(lead_resize)

    np = intf%np_uc
    np_part = np + lead%np_part - lead%np

    if(np.ne.lead%np) then
      ! the new unit cell must be smaller
      ASSERT(np.le.lead%np)
      ! create temp lead
      call lead_init_pot(old_lead, np, np_part, nspin)
      call lead_init_kin(old_lead, np, np_part, dim)
      ! safe old lead
      call lead_copy(lead, old_lead, dim, nspin)
      ! delete lead
      call lead_end(lead)
      ! allocate new (smaller) lead
      call lead_init_pot(lead, np, np_part, nspin)
      call lead_init_kin(lead, np, np_part, dim)
      ! copy parts of the old leads
      call lead_copy(old_lead, lead, dim, nspin)
      ! delete old lead
      call lead_end(old_lead)
    end if

    POP_SUB(lead_resize)
  end subroutine lead_resize


  ! ---------------------------------------------------------
  !> Is the lead potential translationally invariant (in transport direction)?
  logical function is_lead_transl_inv(lapl, vks, intf)
    type(nl_operator_t), intent(in)  :: lapl
    FLOAT,               intent(in)  :: vks(:)
    type(interface_t),   intent(in)  :: intf

    PUSH_SUB(is_lead_transl_inv)

    ! FIXME: For now every potential is translationally invariant in transport direction.
    ! This should be tested and also be generalized to periodic potentials.
    is_lead_transl_inv = .true.

    POP_SUB(is_lead_transl_inv)
  end function is_lead_transl_inv


end module ob_lead_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
