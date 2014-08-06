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
!! $Id$

  ! ---------------------------------------------------------
  !> Crank-Nicolson scheme with source and memory term.
  subroutine td_crank_nicolson_src_mem(hm, gr, st, tr, max_iter, nt, time, dt)
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    integer,                         intent(in)    :: max_iter
    integer,                         intent(in)    :: nt
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt

    PUSH_SUB(td_crank_nicolson_src_mem)
    
    select case(tr%ob%mem_type)
    case(SAVE_CPU_TIME)
      call cn_src_mem_dt(tr%ob, st, hm, gr, dt, time, nt)
    case(SAVE_RAM_USAGE)
      call cn_src_mem_sp_dt(tr%ob, st, hm, gr, max_iter, dt, time, nt)
    end select

    if(.not. hm%cmplxscl%space) then
      call density_calc(st, gr, st%rho)
    else
      call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
    end if

    POP_SUB(td_crank_nicolson_src_mem)
  end subroutine td_crank_nicolson_src_mem
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> Crank-Nicolson propagator
  subroutine td_crank_nicolson(ks, hm, gr, st, tr, time, dt, ions, geo, use_sparskit)
    type(v_ks_t), target,            intent(inout) :: ks
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    logical, intent(in) :: use_sparskit

    CMPLX, allocatable :: zpsi_rhs(:,:), zpsi(:), rhs(:), inhpsi(:)
    integer :: ik, ist, idim, ip, np_part, np, iter
    FLOAT :: dres
    FLOAT :: cgtol = CNST(1.0e-12)
    logical :: converged
    type(ion_state_t) :: ions_state

    PUSH_SUB(propagator_dt.td_crank_nicolson)

#ifndef HAVE_SPARSKIT
    if(use_sparskit) then
      message(1) = "Cannot use SPARSKIT in Crank-Nicolson propagator: not compiled with SPARSKIT support."
      call messages_fatal(1)
    endif
#endif

    np_part = gr%mesh%np_part
    np = gr%mesh%np

    ! define pointer and variables for usage in td_zop, td_zopt routines
    grid_p    => gr
    hm_p      => hm
    tr_p      => tr
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
    if(ion_dynamics_ions_move(ions)) then
      call ion_dynamics_save_state(ions, geo, ions_state)
      call ion_dynamics_propagate(ions, gr%sb, geo, time - dt/M_TWO, M_HALF*dt)
      call hamiltonian_epot_generate(hm, gr, geo, st, time = time - dt/M_TWO)
    end if

    if(hm%cmplxscl%space) then
    call vksinterp_interpolate(tr%vksold, 3, gr%mesh%np, st%d%nspin, &
      time, dt, time -dt/M_TWO, hm%vhxc, hm%imvhxc)
    else
    call vksinterp_interpolate(tr%vksold, 3, gr%mesh%np, st%d%nspin, &
      time, dt, time -dt/M_TWO, hm%vhxc)
    end if
    
    call hamiltonian_update(hm, gr%mesh, time = time - dt/M_TWO)
     
    ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end

        call states_get_state(st, gr%mesh, ist, ik, zpsi_rhs)
        call exponential_apply(tr%te, gr%der, hm, zpsi_rhs, ist, ik, dt/M_TWO, time - dt/M_TWO)

        if(hamiltonian_inh_term(hm)) then
          SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
          do idim = 1, st%d%dim
            call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
            forall(ip = 1:gr%mesh%np) zpsi_rhs(ip, idim) = zpsi_rhs(ip, idim) + dt*inhpsi(ip)
          end do
          SAFE_DEALLOCATE_A(inhpsi)
        end if

        ! put the values in a continuous array
        do idim = 1, st%d%dim
          call states_get_state(st, gr%mesh, idim, ist, ik, zpsi((idim - 1)*np+1:idim*np))
          rhs((idim - 1)*np + 1:idim*np) = zpsi_rhs(1:np, idim)
        end do

        ist_op = ist
        ik_op = ik

        if(use_sparskit) then
#ifdef HAVE_SPARSKIT
          call zsparskit_solver_run(tr%tdsk, td_zop, td_zopt, zpsi, rhs)
#endif
        else
          iter = 2000
          call zqmr_sym_gen_dotu(np*st%d%dim, zpsi, rhs, propagator_qmr_op, zmf_dotu_aux, zmf_nrm2_aux, &
            propagator_qmr_prec, iter, dres, cgtol, showprogress = .false., converged = converged)

          if(.not.converged) then
            write(message(1),'(a)')        'The linear solver used for the Crank-Nicolson'
            write(message(2),'(a,es14.4)') 'propagator did not converge: Residual = ', dres
            call messages_warning(2)
          end if

        endif

        do idim = 1, st%d%dim
          call states_set_state(st, gr%mesh, idim, ist, ik, zpsi((idim-1)*np + 1:(idim - 1)*np + np))
        end do          

      end do
    end do

    if(hm%cmplxscl%space) then !Left states
        
      dt_op = - dt !propagate backwards
      t_op  = time + dt/M_TWO
        
      call vksinterp_interpolate(tr%vksold, 3,  &
        gr%mesh%np, st%d%nspin, time, dt, time + dt/M_TWO, hm%vhxc, hm%imvhxc)

      call hamiltonian_update(hm, gr%mesh, time = time + dt/M_TWO)
      
      
      ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end

          call states_get_state(st, gr%mesh, ist, ik, zpsi_rhs,left = .true. )
          call exponential_apply(tr%te, gr%der, hm, zpsi_rhs, ist, ik, -dt/M_TWO, time + dt/M_TWO)

          if(hamiltonian_inh_term(hm)) then
            SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
            do idim = 1, st%d%dim
              call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi, left = .true.)
              forall(ip = 1:gr%mesh%np) zpsi_rhs(ip, idim) = zpsi_rhs(ip, idim) - dt*inhpsi(ip)
            end do
            SAFE_DEALLOCATE_A(inhpsi)
          end if

          ! put the values in a continuous array
          do idim = 1, st%d%dim
            call states_get_state(st, gr%mesh, idim, ist, ik, zpsi((idim - 1)*np+1:idim*np), left = .true.)
            rhs((idim - 1)*np + 1:idim*np) = zpsi_rhs(1:np, idim)
          end do

          ist_op = ist
          ik_op = ik
          if(use_sparskit) then
#ifdef HAVE_SPARSKIT
            call zsparskit_solver_run(tr%tdsk, td_zop, td_zopt, zpsi, rhs)
#endif
          else
            iter = 2000
            call zqmr_sym_gen_dotu(np*st%d%dim, zpsi, rhs, propagator_qmr_op, zmf_dotu_aux, zmf_nrm2_aux, &
              propagator_qmr_prec, iter, dres, cgtol, showprogress = .false., converged = converged)

            if(.not.converged) then
              write(message(1),'(a)')        'The linear solver used for the Crank-Nicolson'
              write(message(2),'(a,es14.4)') 'propagator did not converge for left states: Residual = ', dres
              call messages_warning(2)
            end if

          endif

          do idim = 1, st%d%dim
            call states_set_state(st, gr%mesh, idim, ist, ik, zpsi((idim-1)*np + 1:(idim - 1)*np + np), left = .true.)
          end do

        end do
      end do
        
    end if

    if(.not. hm%cmplxscl%space) then
      call density_calc(st, gr, st%rho)
    else
      call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
    end if

    !restore to time 'time - dt'
    if(ion_dynamics_ions_move(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

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

    SAFE_ALLOCATE(zpsi(1:grid_p%mesh%np_part, 1:dim_op))
    zpsi = M_z0
    forall(idim = 1:dim_op)
      zpsi(1:grid_p%mesh%np, idim) = xre((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) + &
                                     M_zI * xim((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np)
    end forall

    call exponential_apply(tr_p%te, grid_p%der, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO, t_op)

    forall(idim = 1:dim_op)
      yre((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = real(zpsi(1:grid_p%mesh%np, idim))
      yim((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = aimag(zpsi(1:grid_p%mesh%np, idim))
    end forall

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
    SAFE_ALLOCATE(zpsi(1:grid_p%mesh%np_part, 1:dim_op))
    zpsi = M_z0
    forall(idim = 1:dim_op)
      zpsi(1:grid_p%mesh%np, idim) = xre((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) - &
                                     M_zI * xim((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np)
    end forall

    call exponential_apply(tr_p%te, grid_p%der, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO, t_op)

    forall(idim = 1:dim_op)
      yre((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) =    real(zpsi(1:grid_p%mesh%np, idim))
      yim((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = - aimag(zpsi(1:grid_p%mesh%np, idim))
    end forall

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

    SAFE_ALLOCATE(zpsi(1:grid_p%mesh%np_part, 1:dim_op))
    zpsi = M_z0
    forall(idim = 1:dim_op)
      zpsi(1:grid_p%mesh%np, idim) = x((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np)
    end forall

    call exponential_apply(tr_p%te, grid_p%der, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO, t_op)

    forall(idim = 1:dim_op)
      y((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = zpsi(1:grid_p%mesh%np, idim)
    end forall

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


