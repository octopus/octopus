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
!! $Id: propagator.F90 5406 2009-05-16 18:17:47Z xavier $


  ! ---------------------------------------------------------
  !> Propagator specifically designed for the QOCT+TDDFT problem
  subroutine td_qoct_tddft_propagator(hm, gr, st, tr, t, dt)!, gauge_force, ions, geo)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    type(propagator_t),  intent(inout) :: tr
    FLOAT,               intent(in)    :: t, dt

    PUSH_SUB(td_qoct_tddft_propagator)
    
    if(hm%theory_level .ne. INDEPENDENT_PARTICLES) then
      call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), t-dt/M_TWO, hm%vhxc(:, :))
    end if

    call hamiltonian_update(hm, gr%mesh, time = t-dt/M_TWO)
    call exponential_apply_all(tr%te, gr%der, hm, st, dt, t - dt/M_TWO)

    POP_SUB(td_qoct_tddft_propagator)
  end subroutine td_qoct_tddft_propagator
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Propagator specifically designed for the QOCT+TDDFT problem (2nd version)
  !!
  !! \warning  This is *not* going to work if parallelization in states is used.
  subroutine td_qoct_tddft_propagator_2(hm, gr, st, tr, t, dt)
    type(hamiltonian_t), target, intent(inout) :: hm
    type(grid_t), target,        intent(inout) :: gr
    type(states_t),              intent(inout) :: st
    type(propagator_t), target,      intent(inout) :: tr
    FLOAT,                       intent(in)    :: t, dt

    CMPLX, allocatable :: zpsi(:), rhs(:)
    integer :: ik, ist, idim, isize, np_part, np, iter
    FLOAT :: dres
    FLOAT :: cgtol = CNST(1.0e-8)
    logical :: converged

    PUSH_SUB(td_qoct_tddft_propagator_2)

    np_part = gr%mesh%np_part
    np = gr%mesh%np
    isize = np_part*st%lnst*st%d%kpt%nlocal*st%d%dim

    ! define pointer and variables for usage in td_zop, td_zopt routines
    grid_p    => gr
    hm_p      => hm
    tr_p      => tr
    dt_op = dt
    t_op  = t-dt
    dim_op = st%d%dim
    nst_op = st%nst
    call states_copy(st_op, st)

    ! we (ab)use exponential_apply to compute (1-i\delta t/2 H_n)\psi^n
    ! exponential order needs to be only 1
    tr%te%exp_method = EXP_TAYLOR
    tr%te%exp_order  = 1

    SAFE_ALLOCATE(zpsi(1:np*st%d%dim*st%nst))
    SAFE_ALLOCATE(rhs(1:np*st%d%dim*st%nst))
        
    call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), t - dt/M_TWO, hm%vhxc(:, :))
    call hamiltonian_update(hm, gr%mesh, time = t - dt/M_TWO)

    call exponential_apply_all(tr%te, gr%der, hm, st_op, dt/M_TWO, t-dt)

    ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          call states_get_state(st, gr%mesh, idim, ist, ik, &
            zpsi((ist - 1)*np*st%d%dim + (idim - 1)*np + 1:(ist - 1)*np*st%d%dim + idim*np))
          rhs((ist-1)*np*st%d%dim + (idim-1)*np +1:(ist-1)*np*st%d%dim + idim*np) = &
            st_op%zpsi(1:np, idim, ist, ik)
        end do
      end do
    end do

    ik_op = 1
    iter = 2000
    call zqmr_sym(st%nst*np*st%d%dim, zpsi, rhs, propagator_qmr2_op, propagator_qmr_prec, iter, dres, cgtol, &
                  showprogress = .false., converged = converged)

    if(.not.converged) then
      write(message(1),'(a)')        'The linear solver used for the Crank-Nicholson'
      write(message(2),'(a,es14.4)') 'propagator did not converge: Residual = ', dres
      call messages_warning(2)
    end if

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          call states_set_state(st, gr%mesh, idim, ist, ik, &
            zpsi((ist - 1)*np*st%d%dim + (idim - 1)*np + 1:(ist - 1)*np*st%d%dim + idim*np))
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(rhs)
    call states_end(st_op)
    POP_SUB(td_qoct_tddft_propagator_2)
  end subroutine td_qoct_tddft_propagator_2
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> operator for Crank-Nicholson scheme
  subroutine propagator_qmr2_op(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)
    integer :: idim, ist, np

    np = grid_p%mesh%np
    forall(ist = 1:nst_op, idim = 1:dim_op)
      st_op%zpsi(1:np, idim, ist, ik_op)= x((ist-1)*np*dim_op + (idim-1)*np+1: &
                                            (ist-1)*np*dim_op + (idim-1)*np+np)
    end forall

    call exponential_apply_all(tr_p%te, grid_p%der, hm_p, st_op, -dt_op/M_TWO, t_op)

    forall(ist = 1:nst_op, idim = 1:dim_op)
      y((ist-1)*np*dim_op + (idim-1)*np+1:(ist-1)*np*dim_op + (idim-1)*np+np) = &
        st_op%zpsi(1:np, idim, ist, ik_op)
    end forall

  end subroutine propagator_qmr2_op
  ! ---------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
