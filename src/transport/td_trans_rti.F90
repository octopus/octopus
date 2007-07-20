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

! Modified Crank-Nicholson propagator with source and memory term.

#include "global.h"

module td_trans_rti_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lib_basic_alg_m
  use messages_m
  use sparskit_m
  use states_m
  use td_exp_m
  use td_rti_m
  use td_trans_intf_m
  use v_ks_m

  implicit none

#ifdef HAVE_SPARSKIT

  private
  public ::          &
    cn_src_mem_init, &
    cn_src_mem_end

  type(sparskit_solver_t), pointer :: tdsk
  type(td_exp_t)                   :: matr_exp
  CMPLX, allocatable               :: zpsi_tmp(:, :)

  type(hamiltonian_t), pointer :: h_p
  type(grid_t), pointer        :: gr_p
  CMPLX, pointer               :: mem_p(:, :, :, :)
  type(intface_t), pointer     :: intf_p
  FLOAT                        :: dt_op, t_op
  integer                      :: ik_op

contains

  ! ---------------------------------------------------------
  ! Initialize propagator.
  subroutine cn_src_mem_init(gr)
    type(grid_t), intent(in) :: gr

    call push_sub('td_trans_rti.cn_src_mem_init')

    call zsparskit_solver_init(NP, tdsk)

    ! Set up the matrix exponential to calculate 1-i\Delta t/2 H with
    ! first order Taylor expansion (like in Crank-Nicholson in td_rti.td_rti_dt).
    matr_exp%exp_method = TAYLOR
    matr_exp%exp_order  = 1

    call pop_sub()
  end subroutine cn_src_mem_init


  ! ---------------------------------------------------------
  ! Finish propagator.
  subroutine cn_src_mem_end()
    call push_sub('td_trans_rti.cn_src_mem_end')

    call zsparskit_solver_end()

    call pop_sub()
  end subroutine cn_src_mem_end


  ! ---------------------------------------------------------
  ! Crank-Nicholson timestep with source and memory.
  ! Only non-interacting electrons for the moment, so no
  ! predictor-corrector scheme.
  subroutine cn_src_mem_dt(intf, mem, st_intf, ks, h, st, gr, dt, t, max_iter)
    type(intface_t), target,     intent(in)    :: intf
    CMPLX, target,               intent(in)    :: mem(intf%np, intf%np, 0:max_iter, NLEADS)
    CMPLX,                       intent(inout) :: st_intf(:, :, :, :)
    type(v_ks_t),                intent(in)    :: ks
    type(hamiltonian_t), target, intent(inout) :: h
    type(states_t),              intent(inout) :: st
    type(grid_t), target,        intent(inout) :: gr
    FLOAT,                       intent(in)    :: dt
    FLOAT,                       intent(in)    :: t
    integer,                     intent(in)    :: max_iter

    integer :: ik, ist, il, it, n

    call push_sub('td_trans_rti.cn_src_mem_dt')

    ALLOCATE(zpsi_tmp(NP_PART, st%d%ispin), NP_PART*st%d%ispin)

    ! Set pointers to communicate with operators for sparskit solver.
    h_p    => h
    gr_p   => gr
    mem_p  => mem
    intf_p => intf
    dt_op  =  dt
    t_op   =  t

    n = t/dt

    ! Get right-hand side.
    ! 1. Apply effective Hamiltonian.
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        call apply_h_eff(h, gr, mem(:, :, 0, LEFT), mem(:, :, 0, RIGHT), intf, &
          -M_ONE, dt, t, ik, st%zpsi(:, :, ist, ik))
      end do
    end do

    ! 2. Add source term (to be implemented later).

    ! 3. Add memory term.
    ! Note: there is no potential for the leads at the moment. So, the
    ! prefactors of the potential are just 1.
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do il = 1, NLEADS
          do it = 0, n
            call apply_mem(mem(:, :, n-it, il), il, dt, intf, &
              st_intf(:, it+1, ist, ik), st%zpsi(:, :, ist, ik))
            call apply_mem(mem(:, :, n-it, il), il, dt, intf, &
              st_intf(:, it, ist, ik), st%zpsi(:, :, ist, ik))
            call apply_mem(mem(:, :, n-it-1, il), il, dt, intf, &
              st_intf(:, it+1, ist, ik), st%zpsi(:, :, ist, ik))
            call apply_mem(mem(:, :, n-it-1, il), il, dt, intf, &
              st_intf(:, it, ist, ik), st%zpsi(:, :, ist, ik))
          end do
        end do
      end do
    end do

    ! 4. Solve linear system (1 + i \Delta H_{eff}) zpsi_tmp = st%zpsi.
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        ik_op = ik
        call zsparskit_solver_run(tdsk, skop, skopt, zpsi_tmp(:, 1), st%zpsi(:, 1, ist, ik))
        call lalg_copy(NP, zpsi_tmp(:, 1), st%zpsi(:, 1, ist, ik))
      end do
    end do

    deallocate(zpsi_tmp)

    call pop_sub()

  end subroutine cn_src_mem_dt


  ! ---------------------------------------------------------
  ! Propagate forward/backward with effective Hamiltonian.
  subroutine apply_h_eff(h, gr, mem_l, mem_r, intf, sign, dt, t, ik, zpsi)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(inout) :: gr
    CMPLX, target,       intent(in)    :: mem_l(:, :), mem_r(:, :)
    type(intface_t),     intent(in)    :: intf
    FLOAT,               intent(in)    :: sign, dt, t
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)

    CMPLX, allocatable :: intf_wf(:)

    call push_sub('td_trans_rti.h_eff')

    ALLOCATE(intf_wf(intf%np), intf%np)

    ASSERT(sign.eq.M_ONE.or.sign.eq.-M_ONE)

    ! Apply Hamiltonian of central region.
    call td_exp_dt(matr_exp, gr, h, zpsi, ik, sign*dt/M_TWO, t-dt)

    ! Apply modification: - (\Delta t/2)^2 Q zpsi
    call get_intf_wf(intf, LEFT, zpsi, intf_wf)
    call apply_mem(mem_l, LEFT, dt, intf, intf_wf, zpsi)
    call get_intf_wf(intf, RIGHT, zpsi, intf_wf)
    call apply_mem(mem_r, RIGHT, dt, intf, intf_wf, zpsi)

    deallocate(intf_wf)

    call pop_sub()
  end subroutine apply_h_eff


  ! ---------------------------------------------------------
  ! Apply memory coefficient: zpsi <- zpsi - (\Delta t/2)^2 Q intf_wf
  subroutine apply_mem(mem, il, dt, intf, intf_wf, zpsi)
    CMPLX,           intent(in)    :: mem(:, :)
    integer,         intent(in)    :: il
    FLOAT,           intent(in)    :: dt
    type(intface_t), intent(in)    :: intf
    CMPLX,           intent(in)    :: intf_wf(:)
    CMPLX,           intent(inout) :: zpsi(:, :)

    CMPLX, allocatable :: mem_intf_wf(:)

    call push_sub('td_trans_rti.apply_mem')

    ALLOCATE(mem_intf_wf(intf%np), intf%np)

    call lalg_gemv(intf%np, intf%np, cmplx(-dt**M_TWO/M_FOUR, M_ZERO, REAL_PRECISION), &
      mem, intf_wf, cmplx(M_TWO, M_ZERO, REAL_PRECISION), mem_intf_wf)

    call lalg_axpy(intf%np, M_ONE, mem_intf_wf, &
      zpsi(intf%index_range(1, il):intf%index_range(2, il), 1))

    deallocate(mem_intf_wf)

    call pop_sub()
  end subroutine apply_mem
  
  ! ---------------------------------------------------------
  ! Operator for sparskit solver: 1 + i \Delta H_{eff}
  subroutine skop(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:), xim(:)
    FLOAT, intent(out) :: yre(:), yim(:)

    call push_sub('td_trans_rti.skop')

    zpsi_tmp(1:gr_p%m%np, 1) = xre(1:gr_p%m%np) + M_zI*xim(1:gr_p%m%np)

    ! Propagate backward.
    call apply_h_eff(h_p, gr_p, mem_p(:, :, 0, LEFT), mem_p(:, :, 0, RIGHT), intf_p, &
      M_ONE, dt_op, t_op, ik_op, zpsi_tmp(:, :))

    yre(1:gr_p%m%np) = real(zpsi_tmp(1:gr_p%m%np, 1))
    yim(1:gr_p%m%np) = aimag(zpsi_tmp(1:gr_p%m%np, 1))

    call pop_sub()
  end subroutine skop


  ! ---------------------------------------------------------
  ! Transposed operator for sparskit solver: (1 + i \Delta H_{eff})^T
  subroutine skopt(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:), xim(:)
    FLOAT, intent(out) :: yre(:), yim(:)

    call push_sub('td_trans_rti.skopt')

    ! To act with the transpose of H_{eff} on the wavefunction we
    ! apply H_{eff} to the conjugate of psi and conjugate the
    ! resulting H_{eff} psi (note that H_{eff} is not a purely real
    ! operator for scattering wavefunctions anymore).
    zpsi_tmp(1:gr_p%m%np, 1) = xre(1:gr_p%m%np) - M_zI*xim(1:gr_p%m%np)

    call apply_h_eff(h_p, gr_p, mem_p(:, :, 0, LEFT), mem_p(:, :, 0, RIGHT), intf_p, &
      M_ONE, dt_op, t_op, ik_op, zpsi_tmp(:, :))

    yre(1:gr_p%m%np) = real(zpsi_tmp(1:gr_p%m%np, 1))
    yim(1:gr_p%m%np) = -aimag(zpsi_tmp(1:gr_p%m%np, 1))

    call pop_sub()
  end subroutine skopt

#endif

end module td_trans_rti_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
