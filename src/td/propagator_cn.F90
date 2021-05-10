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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module propagator_cn_oct_m
  use density_oct_m
  use exponential_oct_m
  use gauge_field_oct_m
  use grid_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_base_oct_m
  use propagation_ops_elec_oct_m
  use solvers_oct_m
  use space_oct_m
  use sparskit_oct_m
  use states_elec_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                    &
    td_crank_nicolson

  type(namespace_t),        pointer, private :: namespace_p
  type(mesh_t),             pointer, private :: mesh_p
  type(hamiltonian_elec_t), pointer, private :: hm_p
  type(propagator_base_t),  pointer, private :: tr_p
  integer,                           private :: ik_op, ist_op, dim_op
  FLOAT,                             private :: t_op, dt_op

contains

  ! ---------------------------------------------------------
  !> Crank-Nicolson propagator
  subroutine td_crank_nicolson(hm, namespace, space, gr, st, tr, time, dt, ions_dyn, ions, use_sparskit)
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(namespace_t),        target, intent(in)    :: namespace
    type(space_t),                    intent(in)    :: space
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_base_t),  target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    type(ion_dynamics_t),             intent(inout) :: ions_dyn
    type(ions_t),                     intent(inout) :: ions
    logical,                          intent(in) :: use_sparskit

    CMPLX, allocatable :: zpsi_rhs(:,:), zpsi(:), rhs(:), inhpsi(:)
    integer :: ik, ist, idim, ip, np_part, np, iter
    FLOAT :: dres
    FLOAT :: cgtol = CNST(1.0e-12)
    logical :: converged

    PUSH_SUB(propagator_dt.td_crank_nicolson)

    !TODO: Add gauge field support
    ASSERT(.not.gauge_field_is_applied(hm%ep%gfield))

#ifndef HAVE_SPARSKIT
    if(use_sparskit) then
      message(1) = "Cannot use SPARSKIT in Crank-Nicolson propagator: not compiled with SPARSKIT support."
      call messages_fatal(1, namespace=namespace)
    end if
#endif

    np_part = gr%mesh%np_part
    np = gr%mesh%np

    ! define pointer and variables for usage in td_zop, td_zopt routines
    namespace_p => namespace
    mesh_p      => gr%mesh
    hm_p        => hm
    tr_p        => tr
    dt_op = dt
    t_op  = time - dt/M_TWO
    dim_op = st%d%dim

    ! we (ab)use exponential_apply to compute (1-i\delta t/2 H_n)\psi^n
    ! exponential order needs to be only 1
    tr%te%exp_method = EXP_TAYLOR
    tr%te%exp_order  = 1

    SAFE_ALLOCATE(zpsi_rhs(1:np_part, 1:st%d%dim))
    SAFE_ALLOCATE(zpsi(1:np*st%d%dim))
    SAFE_ALLOCATE(rhs(1:np*st%d%dim))

    !move the ions to time 'time - dt/2', and save the current status to return to it later.
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, space, ions_dyn, ions, &
                time - M_HALF*dt, M_HALF*dt, save_pos = .true.)

    if (family_is_mgga_with_exc(hm%xc)) then
      call potential_interpolation_interpolate(tr%vksold, 3, &
        time, dt, time -dt/M_TWO, hm%vhxc, vtau = hm%vtau)
    else 
      call potential_interpolation_interpolate(tr%vksold, 3, &
        time, dt, time -dt/M_TWO, hm%vhxc)
    end if

    call propagation_ops_elec_update_hamiltonian(namespace, space, st, gr%mesh, hm, time - dt*M_HALF)

    ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end

        call states_elec_get_state(st, gr%mesh, ist, ik, zpsi_rhs)
        call exponential_apply(tr%te, namespace, gr%mesh, hm, zpsi_rhs, ist, ik, dt/M_TWO)

        if(hamiltonian_elec_inh_term(hm)) then
          SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
          do idim = 1, st%d%dim
            call states_elec_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
            do ip = 1, gr%mesh%np
              zpsi_rhs(ip, idim) = zpsi_rhs(ip, idim) + dt*inhpsi(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(inhpsi)
        end if

        ! put the values in a continuous array
        do idim = 1, st%d%dim
          call states_elec_get_state(st, gr%mesh, idim, ist, ik, zpsi((idim - 1)*np+1:idim*np))
          rhs((idim - 1)*np + 1:idim*np) = zpsi_rhs(1:np, idim)
        end do

        ist_op = ist
        ik_op = ik

        if(use_sparskit) then
          call zsparskit_solver_run(tr%tdsk, td_zop, td_zopt, zpsi, rhs)
        else
          iter = 2000
          call zqmr_sym_gen_dotu(np*st%d%dim, zpsi, rhs, propagator_qmr_op, zmf_dotu_aux, zmf_nrm2_aux, &
            propagator_qmr_prec, iter, dres, cgtol, showprogress = .false., converged = converged)

          if(.not.converged) then
            write(message(1),'(a)')        'The linear solver used for the Crank-Nicolson'
            write(message(2),'(a,es14.4)') 'propagator did not converge: Residual = ', dres
            call messages_warning(2, namespace=namespace)
          end if

        end if

        do idim = 1, st%d%dim
          call states_elec_set_state(st, gr%mesh, idim, ist, ik, zpsi((idim-1)*np + 1:(idim - 1)*np + np))
        end do

      end do
    end do

    call density_calc(st, gr, st%rho)

    !restore to time 'time - dt'
    call propagation_ops_elec_restore_ions(tr%propagation_ops_elec, ions_dyn, ions)

    SAFE_DEALLOCATE_A(zpsi_rhs)
    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(rhs)
    POP_SUB(propagator_dt.td_crank_nicolson)

  end subroutine td_crank_nicolson
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> operators for Crank-Nicolson scheme
  subroutine td_zop(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(td_zop)

    SAFE_ALLOCATE(zpsi(1:mesh_p%np_part, 1:dim_op))
    zpsi = M_z0
    do idim = 1, dim_op
      zpsi(1:mesh_p%np, idim) = xre((idim-1)*mesh_p%np+1:idim*mesh_p%np) + &
        M_zI * xim((idim-1)*mesh_p%np+1:idim*mesh_p%np)
    end do

    call exponential_apply(tr_p%te, namespace_p, mesh_p, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO)

    do idim = 1, dim_op
      yre((idim-1)*mesh_p%np+1:idim*mesh_p%np) = TOFLOAT(zpsi(1:mesh_p%np, idim))
      yim((idim-1)*mesh_p%np+1:idim*mesh_p%np) = aimag(zpsi(1:mesh_p%np, idim))
    end do

    SAFE_DEALLOCATE_A(zpsi)

    POP_SUB(td_zop)
  end subroutine td_zop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Transpose of H (called e.g. by bi-conjugate gradient solver)
  subroutine td_zopt(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(td_zopt)

    ! To act with the transpose of H on the wfn we apply H to the conjugate of psi
    ! and conjugate the resulting hpsi (note that H is not a purely real operator
    ! for scattering wavefunctions anymore).
    SAFE_ALLOCATE(zpsi(1:mesh_p%np_part, 1:dim_op))
    zpsi = M_z0
    do idim = 1, dim_op
      zpsi(1:mesh_p%np, idim) = xre((idim-1)*mesh_p%np+1:idim*mesh_p%np) - &
        M_zI * xim((idim-1)*mesh_p%np+1:idim*mesh_p%np)
    end do

    call exponential_apply(tr_p%te, namespace_p, mesh_p, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO)

    do idim = 1, dim_op
      yre((idim-1)*mesh_p%np+1:idim*mesh_p%np) =   TOFLOAT(zpsi(1:mesh_p%np, idim))
      yim((idim-1)*mesh_p%np+1:idim*mesh_p%np) = - aimag(zpsi(1:mesh_p%np, idim))
    end do

    SAFE_DEALLOCATE_A(zpsi)
    POP_SUB(td_zopt)
  end subroutine td_zopt
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> operators for Crank-Nicolson scheme
  subroutine propagator_qmr_op(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(propagator_qmr_op)

    SAFE_ALLOCATE(zpsi(1:mesh_p%np_part, 1:dim_op))
    zpsi = M_z0
    do idim = 1, dim_op
      zpsi(1:mesh_p%np, idim) = x((idim-1)*mesh_p%np+1:idim*mesh_p%np)
    end do

    call exponential_apply(tr_p%te, namespace_p, mesh_p, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO)

    do idim = 1, dim_op
      y((idim-1)*mesh_p%np+1:idim*mesh_p%np) = zpsi(1:mesh_p%np, idim)
    end do

    SAFE_DEALLOCATE_A(zpsi)
    POP_SUB(propagator_qmr_op)
  end subroutine propagator_qmr_op
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_qmr_prec(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    PUSH_SUB(propagator_qmr_prec)
    y = x

    POP_SUB(propagator_qmr_prec)
  end subroutine propagator_qmr_prec
  ! ---------------------------------------------------------

end module propagator_cn_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
