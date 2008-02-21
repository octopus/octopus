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
  use datasets_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_parser_m
  use math_m
  use mesh_function_m
  use messages_m
  use sparskit_m
  use states_m
  use td_exp_m
  use td_trans_mem_m
  use td_trans_lead_m
  use td_trans_src_m
  use td_trans_intf_m
  use v_ks_m
  use nl_operator_m

  implicit none

  private

  integer, parameter :: &
    mem_term_flag  = 1, &
    src_term_flag  = 2, &
    save_cpu_time  = 1, &
    save_ram_usage = 2

  ! Pointers for the h_eff_backward(t) operator for the iterative linear solver.
  type(hamiltonian_t), pointer :: h_p
  type(grid_t), pointer        :: gr_p
  CMPLX, pointer               :: mem_p(:, :, :), sp_mem_p(:, :), mem_s_p(:, :, :, :)
  type(intface_t), pointer     :: intf_p
  FLOAT, pointer               :: dt_p, t_p
  integer, pointer             :: ist_p, ik_p, mem_type_p, mapping_p(:)

  public ::          &
    cn_src_mem_init, &
    cn_src_mem_end,  &
    cn_src_mem_dt,   &
    transport_t

  type transport_t
    integer          :: mem_type              ! 1: fast/lots of memory   2: slow/little memory
    integer          :: sp_length             ! length of the sparse array (as 1d array)
    integer          :: additional_terms      ! Shall we add source and memory term?
    CMPLX, pointer   :: mem_coeff(:, :, :, :) ! i, j, t, il (if mem_type=1)
    CMPLX, pointer   :: mem_sp_coeff(:, :, :) ! sp_length, t, il (if mem_type=1)
    integer, pointer :: sp2full_map(:)        ! the mapping indices from sparse to full matrices
    CMPLX, pointer   :: src_mem_u(:, :)       ! t, il     u(m)
    CMPLX, pointer   :: src_factor(:)         ! iter
    CMPLX, pointer   :: offdiag(:, :, :)      ! hopping operator V^T (np,np,nleads)
    CMPLX, pointer   :: mem_s(:, :, :, :)     ! the matrices to diagonalize coeff0

    type(intface_t)  :: intface

    CMPLX, pointer   :: st_intface(:, :, :, :, :) ! np_intf, nst, nik, nleads, max_iter
    CMPLX, pointer   :: st_sincos(:, :, :) ! np, (sin=1;cos=2,modsin=3,modcos=4), nleads
    FLOAT, pointer   :: st_phase(:, :, :)   ! nst, nik, nleads; phaseshift of lead eigenfunction at the interface

    type(trans_lead_t) :: lead
  end type transport_t

  ! We use 1st order Taylor expansion of exp(-i dt/2 H)|psi>
  ! to calculate (1 - i dt H)|psi>.
  type(td_exp_t) :: taylor_1st

  ! Parameters to the BiCG in Crank-Nicholson.
  integer          :: cg_max_iter
  FLOAT            :: cg_tol

contains

  ! ---------------------------------------------------------
  ! Initialize propagator.
  subroutine cn_src_mem_init(st, gr, trans, dt, max_iter)
    type(states_t),    intent(in)  :: st
    type(grid_t),      intent(in)  :: gr
    type(transport_t), intent(out) :: trans
    FLOAT,             intent(in)  :: dt
    integer,           intent(in)  :: max_iter

    integer :: allocsize, id
    FLOAT   :: energy, q
    type(nl_operator_t)  :: op

    call push_sub('td_trans_rti.cn_src_mem_init')

    taylor_1st%exp_method = TAYLOR
    taylor_1st%exp_order = 1

    call intface_init(gr, trans%intface)

    !%Variable TDTransBiCGMaxIter
    !%Type integer
    !%Default 100
    !%Section Transport
    !%Description
    !% Sets the maximum iteration number for the BiCG linear solver in
    !% the Crank-Nicholson procedure.
    !%End
    call loct_parse_int(check_inp('TDTransBiCGMaxIter'), 100, cg_max_iter)
    if(cg_max_iter.le.0) then
      write(message(1), '(a,i6,a)') "Input : '", cg_max_iter, "' is not a valid TDTransBiCGMaxIter."
      message(2) = '(0 < TDTransBiCGMaxIter)'
      call write_fatal(2)
    end if

    !%Variable TDTransBiCGTol
    !%Type integer
    !%Default 1e-12
    !%Section Transport
    !%Description
    !% Sets the convergence tolerance for the residue in the BiCG linear solver
    !% in the Crank-Nicholson procedure.
    !%End
    call loct_parse_float(check_inp('TDTransBiCGMaxIter'), CNST(1e-12), cg_tol)
    if(cg_tol.le.M_ZERO) then
      write(message(1), '(a,f14.6,a)') "Input : '", cg_tol, "' is not a valid TDTransBiCGTol."
      message(2) = '(0 < TDTransBiCGTol)'
      call write_fatal(2)
    end if

    !%Variable TDTransMemType
    !%Type integer
    !%Default save_cpu_time
    !%Section Transport
    !%Description
    !% Decides whether the memory coefficients use lots of RAM (default)
    !% or uses a more compact scheme but with the need of more CPU-cycles.
    !%
    !%Option save_cpu_time 1
    !% Use the memory intensive procedure
    !%Option save_ram_usage 2
    !% Use the RAM saving, CPU intensive procedure
    !%End
    call loct_parse_int(check_inp('TDTransMemType'), save_cpu_time, trans%mem_type)
    if(trans%mem_type.lt.save_cpu_time.or.trans%mem_type.gt.save_ram_usage) then
      message(1) = 'TDTransMemType may either be "save_cpu_time" or "save_ram_usage".'
      call write_fatal(1)
    end if

    !%Variable TDTransAdditionalTerms
    !%Type integer
    !%Default mem_term + src_term
    !%Section Transport
    !%Description
    !% The TD transport propagator inserts to additional terms in
    !% the Crank-Nicholson scheme: source and memory. With this variable,
    !% one or both of them can be switched off.
    !%
    !%Option mem_term 1
    !% If present, include memory term in propagator
    !%Option src_term 2
    !% If present, include source term in propagator
    !%End
    call loct_parse_int(check_inp('TDTransAdditionalTerms'), mem_term_flag+src_term_flag, trans%additional_terms)

    call lead_init(trans%lead, max_iter, dt)

    call td_trans_rti_write_info(trans, st, max_iter, gr%f_der%der_discr%order, trans%mem_type)

    call memory_init(trans%intface, dt/M_TWO, max_iter, gr%f_der%der_discr%lapl,&
                      trans%mem_coeff, trans%mem_sp_coeff, trans%mem_s, trans%offdiag, gr%sb%dim, &
                      gr%sb%h(1), trans%mem_type, gr%f_der%der_discr%order, trans%sp2full_map)

    ! FIXME  to be set in the inp file (fermie-energy)
    q = M_PI/(gr%sb%lsize(1)-(gr%f_der%der_discr%order-1)*gr%sb%h(1))
    energy = M_HALF*q**2
    op = gr%f_der%der_discr%lapl
    ! calculate energy with tight binding model
    !energy = M_HALF*op%w_re(1, 1)
    !do id = 1, gr%f_der%der_discr%order
    !  energy = energy + op%w_re(2*id, 1)*cos(id*q)
    !end do
    !energy = -energy*gr%sb%h(1)**2

!write(*,*) energy, M_HALF*q**2
    !energy = M_ONE-cos(M_PI/(gr%sb%lsize(1)-(gr%f_der%der_discr%order-1)*gr%sb%h(1)))
    do id=2, gr%sb%dim ! FIXME: when the last gridpoint is not the border, recalculate transversal energy
      energy = energy + M_HALF*(M_PI/(M_TWO*(gr%sb%lsize(id)+gr%sb%h(id))))**2
!      energy = energy + M_ONE-cos(M_PI/(M_TWO*(gr%sb%lsize(id)+gr%sb%h(id))))
    end do

    call source_init(trans%src_factor, trans%st_sincos, trans%st_phase, dt, energy, q, &
                     max_iter, trans%intface%np, st%nst, st%d%nik, gr, gr%f_der%der_discr%order)

    !call calculate_initial_states() in init_wfs()
    !call match_initial_states_to_lead() in cn_src_mem_dt
    !write(*,*) 'psi(0)'
    !write(*,*) st%zpsi(:, :, 1, 1)


    ! Allocate memory for the interface wave functions.
    allocsize = trans%intface%np*st%lnst*st%d%nik*(max_iter+1)*NLEADS
    ALLOCATE(trans%st_intface(trans%intface%np, st%st_start:st%st_end, st%d%nik, NLEADS, 0:max_iter), allocsize)
    trans%st_intface = M_z0
    
    ! Allocate memory for the src_mem_u (needed for source and memory term)
    ALLOCATE(trans%src_mem_u(0:max_iter, NLEADS), (max_iter+1)*NLEADS)
    !            /      dt        \   / /       dt        \
    ! u(m, il) = |1 - i -- U(m,il)|  /  | 1 + i -- U(m,il) |
    !            \      4         / /   \       4         /
    trans%src_mem_u(:, :) = (M_z1 - M_zI*dt/M_FOUR*trans%lead%td_pot(:, :)) / &
      (M_z1 + M_zI*dt/M_FOUR*trans%lead%td_pot(:, :))

    call pop_sub()
  end subroutine cn_src_mem_init


  ! ---------------------------------------------------------
  ! Finish propagator.
  subroutine cn_src_mem_end(trans)
    type(transport_t), intent(inout) :: trans

    call push_sub('td_trans_rti.cn_src_mem_end')

    if(associated(trans%st_intface)) then
      deallocate(trans%st_intface)
      nullify(trans%st_intface)
    end if

    if(associated(trans%src_mem_u)) then
      deallocate(trans%src_mem_u)
      nullify(trans%src_mem_u)
    end if

    call intface_end(trans%intface)
    call lead_end(trans%lead)
    call source_end(trans%src_factor, trans%st_sincos, trans%st_phase)
    call memory_end(trans%mem_coeff, trans%mem_sp_coeff, trans%offdiag, trans%mem_s, trans%sp2full_map)

    call pop_sub()
  end subroutine cn_src_mem_end


  ! ---------------------------------------------------------
  ! Write some status information to stdout.
  subroutine td_trans_rti_write_info(trans, st, max_iter, order, mem_type)
    type(transport_t), intent(in) :: trans
    type(states_t),    intent(in) :: st
    integer,           intent(in) :: max_iter
    integer,           intent(in) :: order
    integer,           intent(in) :: mem_type

    character(len=64) :: terms, mem_type_name

    call push_sub('td_transport.td_transport_write_info')

    call messages_print_stress(stdout, 'Transport')

    if(trans%mem_type.eq.save_cpu_time) then
      mem_type_name = 'full'
    end if
    if(trans%mem_type.eq.save_ram_usage) then
      mem_type_name = 'sparse'
    end if

    terms = ''
    if(iand(trans%additional_terms, mem_term_flag).ne.0) then
      terms = trim(terms)//'memory'
    end if
    if(iand(trans%additional_terms, src_term_flag).ne.0) then
      terms = trim(terms)//' source'
    end if
    
    write(message(1), '(a,i10,i10)') 'Points in interface regions:     ', &
      trans%intface%np, trans%intface%np
    write(message(2), '(a,a10)') 'Type of memory coefficients:     ', trim(mem_type_name)
    write(message(3), '(a,f10.3)') 'MBytes required for memory term: ', &
      mbytes_memory_term(max_iter, trans%intface%np, NLEADS, st, trans%mem_type, order)
    write(message(4), '(a,i10)') 'Maximum BiCG iterations:         ', cg_max_iter
    write(message(5), '(a,es10.1)') 'BiCG residual tolerance:         ', cg_tol
    write(message(6), '(a,a20)') 'Included additional terms:       ', trim(terms)
    write(message(7), '(a,2a10)') 'Lead types (L/R):                ', &
      trim(lead_type_names(trans%lead%lead_type(LEFT))),               &
      trim(lead_type_names(trans%lead%lead_type(RIGHT)))
    write(message(8), '(a,2a10)') 'TD lead potential (L/R):         ', &
      trim(trans%lead%td_pot_formula(LEFT)), trim(trans%lead%td_pot_formula(RIGHT))
    call write_info(8, stdout)

    call messages_print_stress(stdout)

    call pop_sub()
  end subroutine td_trans_rti_write_info

  ! ---------------------------------------------------------
  ! calculate phase shift at interface
  subroutine calculate_phase_shift(order, np, dx, lsize, st, u, phase_shift)
    integer,         intent(in)  :: order
    integer,         intent(in)  :: np
    FLOAT,           intent(in)  :: dx
    FLOAT,           intent(in)  :: lsize
    type(states_t),  intent(in)  :: st
    FLOAT,           intent(in)  :: u(:) ! the lead potential at time 0
    FLOAT,           intent(inout) :: phase_shift(:, :, :) ! nst, nik, NLEADS

    integer         :: lshift, length, ik, ist
    FLOAT           :: q, denom, energy
    FLOAT, allocatable :: coeffs(:,:) ! the weights

    call push_sub('td_transport.calculate_phase_shift')

  ! Matching algorithm
  ! 2. calculate phaseshift of the lead eigenfunction (eq. (41) in paper)
  ! 2.1 calculate the second derivatives at the boundaries
  !write(*,*) 'vorher'

    ! the number of derivative coefficients
    length = 2*order+1
    ! with increasing order the interface boundary shifts towards the center
    lshift = order-1
    ALLOCATE(coeffs(length, NLEADS),length*NLEADS)

    q = M_PI/(lsize-lshift*dx)
    energy = M_HALF*q**2

    call deriv_coeffs(order, dx, coeffs(:,LEFT), +1)
    call deriv_coeffs(order, dx, coeffs(:,RIGHT), -1)
    !write(*,*) 'left stencil weights', coeffs(:,LEFT)
    !write(*,*) 'right stencil weights', coeffs(:,RIGHT)

    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        ! left lead
        if (abs(st%zpsi(1+lshift, 1, ist, ik)).lt.CNST(1e-15)) then
          phase_shift(ist, ik, LEFT) = M_ZERO
        else
          ! FIXME: use a two-sided but asymmetric derivative (exacter)
          ! calculate the forward derivative
          denom = sum(log(abs(st%zpsi(1+lshift:length+lshift, 1, ist, ik))**2)*coeffs(:, LEFT))
          if (abs(denom).lt.CNST(1e-15)) then
            phase_shift(ist, ik, LEFT) = M_HALF*M_PI
          else
            ! no tight binding yet
            phase_shift(ist, ik, LEFT) = atan(M_TWO*sqrt(M_TWO*(energy-u(LEFT))) / denom)
          end if
        end if

        ! right lead
        if (abs(st%zpsi(np-lshift, 1, ist, ik)).lt.CNST(1e-15)) then
          phase_shift(ist, ik, RIGHT) = M_ZERO
        else
          ! FIXME: use a two-sided but asymmetric derivative (exacter)
          ! calculate the backward derivative
            denom = sum(log(abs(st%zpsi(np-length+1-lshift:np-lshift, 1, ist, ik))**2)*coeffs(:, RIGHT))
          if (abs(denom).lt.CNST(1e-15)) then
            phase_shift(ist, ik, RIGHT) = sign(M_HALF*M_PI, denom)
          else
            ! no tight binding yet
            phase_shift(ist, ik, RIGHT) = atan(M_TWO*sqrt(M_TWO*(energy-u(RIGHT))) / denom)
          end if
        end if
      end do
    end do

!  st_phase = M_ZERO
    deallocate(coeffs)

    call pop_sub()
  end subroutine calculate_phase_shift

  
  ! ---------------------------------------------------------
  ! Crank-Nicholson timestep with source and memory.
  ! Only non-interacting electrons for the moment, so no
  ! predictor-corrector scheme.
  subroutine cn_src_mem_dt(intf, mem_type, mem, sp_mem, st_intf, ks, h, st, gr, dt,  &
    t, max_iter, sm_u, u, timestep, src_factor, offdiag, st_sincos, st_phase, mem_s, &
    additional_terms, mapping)
    type(intface_t), target,     intent(in)    :: intf
    type(states_t),              intent(inout) :: st
    integer,                     intent(in)    :: max_iter
    CMPLX,                       intent(inout) :: st_intf(intf%np, st%st_start:st%st_end, st%d%nik, NLEADS, 0:max_iter)
    type(v_ks_t),                intent(in)    :: ks
    type(hamiltonian_t), target, intent(inout) :: h
    type(grid_t), target,        intent(inout) :: gr
    FLOAT, target,               intent(in)    :: dt
    FLOAT, target,               intent(in)    :: t
    integer, target,             intent(in)    :: mem_type
    CMPLX, target,               intent(in)    :: mem(intf%np, intf%np, 0:max_iter, NLEADS)
    CMPLX, target,               intent(in)    :: sp_mem(intf%np*gr%f_der%der_discr%order, 0:max_iter, NLEADS)
    CMPLX, target,               intent(in)    :: sm_u(0:max_iter, NLEADS)
    FLOAT, target,               intent(in)    :: u(0:max_iter, NLEADS)
    integer,                     intent(in)    :: timestep
    CMPLX, target,               intent(in)    :: src_factor(0:max_iter) ! max_iter
    CMPLX, target,               intent(in)    :: offdiag(intf%np, intf%np, NLEADS)   ! hopping operator V^T
    CMPLX, target,               intent(in)    :: st_sincos(:, :, :) ! sin & cos vectors
    FLOAT, target,               intent(inout) :: st_phase(:, :, :)   ! phaseshift at interface
    CMPLX, target,               intent(in)    :: mem_s(intf%np, intf%np, 2, NLEADS)
    integer,                     intent(in)    :: additional_terms
    integer, target,             intent(in)    :: mapping(:)   ! the mapping

    integer            :: il, it, m, cg_iter, j, order
    integer, target    :: ist, ik
    FLOAT              :: tol, q, energy, c
    CMPLX              :: factor
    CMPLX, allocatable :: tmp(:, :), tmp_wf(:), tmp_mem(:, :)

    call push_sub('td_trans_rti.cn_src_mem_dt')

    order = gr%f_der%der_discr%order
    ALLOCATE(tmp(NP, st%d%ispin), NP*st%d%ispin)
    ALLOCATE(tmp_wf(intf%np), intf%np)
    if (mem_type.eq.1) then
      ALLOCATE(tmp_mem(intf%np, intf%np), intf%np**2)
    else
      ALLOCATE(tmp_mem(intf%np*order, 1), intf%np*order)
    end if

    ! Set pointers to communicate with with backward propagator passed
    ! to iterative linear solver.
    h_p        => h
    gr_p       => gr
    mem_p      => mem(:, :, 0, :)
    sp_mem_p   => sp_mem(:, 0, :)
    intf_p     => intf
    dt_p       => dt
    t_p        => t
    ist_p      => ist
    ik_p       => ik
    mem_s_p    => mem_s(:, :, :, :)
    mem_type_p => mem_type
    mapping_p  => mapping
    ! For the dot product passed to BiCG routine.
    call mesh_init_mesh_aux(gr%m)

    m = timestep-1
    cg_iter = cg_max_iter

! FIXME: this does NOT belong here, since the extended groundstate has to be calculated and not set as it is now
if (m.eq.0) then
  ! 1. calculate the extended eigenstate (unbounded)

  ! 2. calculate phase shift of the plane waves in the leads
  call calculate_phase_shift(order, NP, gr%sb%h(1), gr%sb%lsize(1), st, u(0,:), st_phase)

  ! Matching algorithm
  ! 2. calculate phaseshift of the lead eigenfunction (eq. (41) in paper)
  ! 2.1 calculate the second derivatives at the boundaries
!  write(*,*) 'vorher'
!      write(*,*) st_phase(LEFT), sin(st_phase(LEFT))**2, abs(st%zpsi(1, 1, ist, ik))**2
!      write(*,*) st_phase(RIGHT), sin(st_phase(RIGHT))**2, log(abs(st%zpsi(NP, 1, ist, ik))**2)
!      c = M_ONE/sqrt(gr%sb%h(1)+abs(st%zpsi(1, 1, ist, ik))**2/(M_TWO*q)+abs(st%zpsi(NP, 1, ist, ik))**2/(M_TWO*q))
!write(*,*) c
!      st%zpsi(:, :, ist, ik) = c*st%zpsi(:, :, ist, ik)

!  st_phase(:) = M_ZERO ! FIXME
end if
!write(*,*) 'left  center', st%zpsi(1, 1, 1, 1)
!write(*,*) 'right center', st%zpsi(NP, 1, 1, 1)

    ! Save interface part of wavefunctions for subsequent iterations
    ! before we overwrite them with the values for the new timestep.
    ! (Copying before the propagation gets the saving right for the
    ! initial state also.)
    call save_intf_wf(intf, timestep-1, st, max_iter, st_intf)

    ! Get right-hand side.

    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        ! 1. Apply effective Hamiltonian.
        if (mem_type.eq.1) then
          call apply_h_eff(h, gr, mem(:, :, 0, :), intf, -M_ONE, dt, t, ist, ik, st%zpsi(:, :, ist, ik))
        else
          call apply_h_eff_sp(h, gr, sp_mem(:, 0, :), intf, -M_ONE, dt, t, ist, ik, &
                              st%zpsi(:, :, ist, ik), mem_s(:, :, :, :), mapping)
        end if

        do il = 1, NLEADS
          ! 2. Add source term
          if(iand(additional_terms, src_term_flag).ne.0) then
            if (mem_type.eq.1) then
              call calc_source_wf(max_iter, m, il, st_phase(ist, ik, il), offdiag(:, :, il), src_factor(:), &
                mem(:, :, :, il), dt, intf%np, st_sincos(:, :, il), tmp_wf(:))
            else
              ! TODO
            end if
            factor = -M_zI*dt*lambda(m, 0, il, max_iter, sm_u) / sm_u(m, il)
            call apply_src(intf, il, factor, tmp_wf, st%zpsi(:, :, ist, ik))
          end if

          ! 3. Add memory term.
          if(iand(additional_terms, mem_term_flag).ne.0) then
            do it = 0, m-1
              factor = -dt**2/M_FOUR*lambda(m, it, il, max_iter, sm_u) / &
                (sm_u(m, il)*sm_u(it, il))
              tmp_wf(:) = st_intf(:, ist, ik, il, it+1) + st_intf(:, ist, ik, il, it)
              if (mem_type.eq.1) then
                tmp_mem = mem(:, :, m-it, il) + mem(:, :, m-it-1, il)
                call apply_mem(tmp_mem, il, intf, tmp_wf, st%zpsi(:, :, ist, ik), factor)
              else
                tmp_mem(:,1) = sp_mem(:, m-it, il) + sp_mem(:, m-it-1, il)
                call apply_sp_mem(tmp_mem(:,1), il, intf, tmp_wf, st%zpsi(:, :, ist, ik),&
                  factor, mem_s(:,:,:,il),order,gr%sb%dim, mapping)
              end if
            end do
          end if
        end do

        ! TODO: for many states parallel
        ! 4. Solve linear system (1 + i \delta H_{eff}) st%zpsi = tmp.
        tmp(1:NP, 1) = st%zpsi(1:NP, 1, ist, ik)
        if (mem_type.eq.1) then
          call zconjugate_gradients(NP, st%zpsi(1:NP, 1, ist, ik), tmp(1:NP, 1), &
            h_eff_backward, h_eff_backwardt, zmf_dotp_aux, cg_iter, threshold=cg_tol)
        else
          call zconjugate_gradients(NP, st%zpsi(1:NP, 1, ist, ik), tmp(1:NP, 1), &
            h_eff_backward_sp, h_eff_backwardt_sp, zmf_dotp_aux, cg_iter, threshold=cg_tol)
        end if
        ! Write warning if BiCG did not converge.
        if(cg_iter.gt.cg_max_iter) then
          write(message(1), '(a,i5,a)') 'BiCG did not converge in ', cg_max_iter, ' iterations'
          message(2) = 'when propagating with H_eff.'
          call write_warning(2)
        end if

      end do
    end do

    ! Save interface part of wavefunction for subsequent iterations.
    call save_intf_wf(intf, timestep, st, max_iter, st_intf)

    deallocate(tmp, tmp_wf, tmp_mem)
    call pop_sub()

  contains
    CMPLX function lambda(m, k, il, max_iter, sm_u) result(res)
      integer, intent(in) :: m, k, il, max_iter
      CMPLX,   intent(in) :: sm_u(0:max_iter, NLEADS)

      integer :: ij
      res = M_z1

      do ij = k, m
        res = res*sm_u(ij,il)**2
      end do
    end function lambda

  end subroutine cn_src_mem_dt


  ! ---------------------------------------------------------
  ! output of a matrix (for debugging)
  subroutine write_matrix(matr,np)
    CMPLX,     intent(in)  :: matr(:, :)
    integer,   intent(in)  :: np
    integer            :: i, j

    do j = 1, np
      do i = 1, np
        write(*,'(a,i3,a,i3,a,2e20.9)') "m(",j,",",i,")=",matr(j,i)
      end do
    end do
   end subroutine write_matrix


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff}
  ! Used by the iterative linear solver.
  subroutine h_eff_backward(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('td_trans_rti.h_eff_backward')
    
    ALLOCATE(tmp(gr_p%m%np, 1),gr_p%m%np)
    ! Propagate backward.
    call lalg_copy(gr_p%m%np, x, tmp(:, 1))
    call apply_h_eff(h_p, gr_p, mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp)
    call lalg_copy(gr_p%m%np, tmp(:, 1), y)

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backward

  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff}
  ! Used by the iterative linear solver.
  subroutine h_eff_backward_sp(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('td_trans_rti.h_eff_backward')
    
    ALLOCATE(tmp(gr_p%m%np, 1),gr_p%m%np)
    ! Propagate backward.
    call lalg_copy(gr_p%m%np, x, tmp(:, 1))
    call apply_h_eff_sp(h_p, gr_p, sp_mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp, mem_s_p, mapping_p)
    call lalg_copy(gr_p%m%np, tmp(:, 1), y)

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backward_sp


  ! ---------------------------------------------------------
  ! Progagate backwards with: (1 + i \delta H_{eff})^\dagger
  ! Used by the iterative linear solver.
  subroutine h_eff_backwardt(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('td_trans_rti.h_eff_backwardt')
    
    ALLOCATE(tmp(gr_p%m%np, 1),gr_p%m%np)
    ! Propagate backward.
    ! To act with the Hermitian conjugate of H_{eff} on the wavefunction
    ! we apply H_{eff} to the conjugate of psi and conjugate the
    ! resulting H_{eff} psi (note that H_{eff} is not a purely real
    ! operator for scattering wavefunctions anymore).
    tmp(1:gr_p%m%np, 1) = conjg(x(1:gr_p%m%np))
    call apply_h_eff(h_p, gr_p, mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp)
    y(1:gr_p%m%np) = conjg(tmp(1:gr_p%m%np, 1))

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backwardt

  ! ---------------------------------------------------------
  ! Progagate backwards with: (1 + i \delta H_{eff})^\dagger
  ! Used by the iterative linear solver.
  subroutine h_eff_backwardt_sp(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('td_trans_rti.h_eff_backwardt')
    
    ALLOCATE(tmp(gr_p%m%np, 1),gr_p%m%np)
    ! Propagate backward.
    ! To act with the Hermitian conjugate of H_{eff} on the wavefunction
    ! we apply H_{eff} to the conjugate of psi and conjugate the
    ! resulting H_{eff} psi (note that H_{eff} is not a purely real
    ! operator for scattering wavefunctions anymore).
    tmp(1:gr_p%m%np, 1) = conjg(x(1:gr_p%m%np))
    call apply_h_eff_sp(h_p, gr_p, sp_mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp, mem_s_p, mapping_p)
    y(1:gr_p%m%np) = conjg(tmp(1:gr_p%m%np, 1))

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backwardt_sp


  ! ---------------------------------------------------------
  ! Save the interface part of all states st for timestep into
  ! the st_intf array.
  subroutine save_intf_wf(intf, timestep, st, max_iter, st_intf)
    type(intface_t), intent(in)    :: intf
    integer,         intent(in)    :: timestep
    type(states_t),  intent(in)    :: st
    integer,         intent(in)    :: max_iter
    CMPLX,           intent(inout) :: st_intf(1:intf%np, st%st_start:st%st_end, 1:st%d%nik, 1:NLEADS, 0:max_iter)

    integer :: ik, ist, il
    
    call push_sub('td_trans_rti.save_intf_wf')

    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do il = 1, NLEADS
          call get_intf_wf(intf, il, st%zpsi(:, :, ist, ik), st_intf(:, ist, ik, il, timestep))
        end do
      end do
    end do

    call pop_sub()
  end subroutine save_intf_wf


  ! ---------------------------------------------------------
  ! Propagate forward/backward with effective Hamiltonian.
  subroutine apply_h_eff(h, gr, mem, intf, sign, dt, t, ist, ik, zpsi)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(inout) :: gr
    CMPLX, target,       intent(in)    :: mem(:, :, :)
    type(intface_t),     intent(in)    :: intf
    FLOAT,               intent(in)    :: sign, dt, t
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)

    CMPLX, allocatable :: intf_wf(:, :)

    call push_sub('td_trans_rti.apply_h_eff')

    ALLOCATE(intf_wf(intf%np, NLEADS), intf%np*NLEADS)

    ASSERT(sign.eq.M_ONE.or.sign.eq.-M_ONE)

    call get_intf_wf(intf, LEFT, zpsi, intf_wf(:, LEFT))
    call get_intf_wf(intf, RIGHT, zpsi, intf_wf(:, RIGHT))

    ! Apply (1 - i\delta H_{CC}^{(m)}) to zpsi.
    ! td_exp_dt with Taylor expansion calculates exp(-i dt H), i. e. the
    ! minus is already built in.
    call td_exp_dt(taylor_1st, gr, h, zpsi, ist, ik, -sign*dt/M_TWO, t-dt)

    ! Apply modification: sign \delta^2 Q zpsi
    call apply_mem(mem(:,:,LEFT), LEFT, intf, intf_wf(:, LEFT), zpsi, TOCMPLX(sign*dt**2/M_FOUR, M_ZERO))
    call apply_mem(mem(:,:,RIGHT), RIGHT, intf, intf_wf(:, RIGHT), zpsi, TOCMPLX(sign*dt**2/M_FOUR, M_ZERO))

    deallocate(intf_wf)

    call pop_sub()
  end subroutine apply_h_eff

  ! ---------------------------------------------------------
  ! Propagate forward/backward with effective Hamiltonian.
  subroutine apply_h_eff_sp(h, gr, sp_mem, intf, sign, dt, t, ist, ik, zpsi, mem_s, mapping)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(inout) :: gr
    CMPLX, target,       intent(in)    :: sp_mem(:, :)
    type(intface_t),     intent(in)    :: intf
    FLOAT,               intent(in)    :: sign, dt, t
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)
    CMPLX,               intent(in)    :: mem_s(:, :, :, :)
    integer,             intent(in)    :: mapping(:)   ! the mapping

    CMPLX, allocatable :: intf_wf(:, :)

    call push_sub('td_trans_rti.apply_h_eff_sp')

    ALLOCATE(intf_wf(intf%np, NLEADS), intf%np*NLEADS)

    ASSERT(sign.eq.M_ONE.or.sign.eq.-M_ONE)

    call get_intf_wf(intf, LEFT, zpsi, intf_wf(:, LEFT))
    call get_intf_wf(intf, RIGHT, zpsi, intf_wf(:, RIGHT))

    ! Apply (1 - i\delta H_{CC}^{(m)}) to zpsi.
    ! td_exp_dt with Taylor expansion calculates exp(-i dt H), i. e. the
    ! minus is already built in.
    call td_exp_dt(taylor_1st, gr, h, zpsi, ist, ik, -sign*dt/M_TWO, t-dt)

    ! Apply modification: sign \delta^2 Q zpsi
    call apply_sp_mem(sp_mem(:,LEFT), LEFT, intf, intf_wf(:, LEFT), zpsi, &
              TOCMPLX(sign*dt**2/M_FOUR, M_ZERO), mem_s(:,:,:,LEFT),gr%f_der%der_discr%order,gr%sb%dim, mapping)
    call apply_sp_mem(sp_mem(:,RIGHT), RIGHT, intf, intf_wf(:, RIGHT), zpsi, &
              TOCMPLX(sign*dt**2/M_FOUR, M_ZERO), mem_s(:,:,:,RIGHT),gr%f_der%der_discr%order,gr%sb%dim, mapping)

    deallocate(intf_wf)

    call pop_sub()
  end subroutine apply_h_eff_sp


  ! ---------------------------------------------------------
  ! Apply source coefficient: zpsi <- zpsi + factor src_wf
  subroutine apply_src(intf, il, factor, src_wf, zpsi)
    type(intface_t), intent(in)    :: intf
    integer,         intent(in)    :: il
    CMPLX,           intent(in)    :: factor
    CMPLX,           intent(in)    :: src_wf(:)
    CMPLX,           intent(inout) :: zpsi(:, :)

    call push_sub('td_trans_rti.apply_src')

    ! Do not use use BLAS here.
    zpsi(intf%index_range(1, il):intf%index_range(2, il), 1) = &
      zpsi(intf%index_range(1, il):intf%index_range(2, il), 1) + factor*src_wf(:)

    call pop_sub()
  end subroutine apply_src

  ! ---------------------------------------------------------
  ! Apply memory coefficient: zpsi <- zpsi + factor Q intf_wf
  ! Use symmetric matrix-vector multiply, because Q is symmetric
  subroutine apply_mem(mem, il, intf, intf_wf, zpsi, factor)
    CMPLX,           intent(in)    :: mem(:, :)
    integer,         intent(in)    :: il
    type(intface_t), intent(in)    :: intf
    CMPLX,           intent(in)    :: intf_wf(:)
    CMPLX,           intent(inout) :: zpsi(:, :)
    CMPLX,           intent(in)    :: factor

    CMPLX, allocatable :: mem_intf_wf(:)
    CMPLX, allocatable :: tmem(:, :, :)

    call push_sub('td_trans_rti.apply_mem')

    ALLOCATE(mem_intf_wf(intf%np), intf%np)
    ALLOCATE(tmem(intf%np, intf%np, 2), 2*intf%np**2)


    mem_intf_wf = M_z0

    call zsymv('U', intf%np, factor, mem, intf%np, intf_wf, 1, M_z0, mem_intf_wf, 1)

    ! Do not use use BLAS here.
    zpsi(intf%index_range(1, il):intf%index_range(2, il), 1) = &
      zpsi(intf%index_range(1, il):intf%index_range(2, il), 1) + mem_intf_wf(:)

    deallocate(mem_intf_wf, tmem)
    call pop_sub()
  end subroutine apply_mem

  ! ---------------------------------------------------------
  ! Apply memory coefficient: zpsi <- zpsi + factor Q intf_wf
  ! Use symmetric matrix-vector multiply, because Q is symmetric
  subroutine apply_sp_mem(sp_mem, il, intf, intf_wf, zpsi, factor, mem_s, order, dim, mapping)
    CMPLX,           intent(in)    :: sp_mem(:)
    integer,         intent(in)    :: il
    type(intface_t), intent(in)    :: intf
    CMPLX,           intent(in)    :: intf_wf(:)
    CMPLX,           intent(inout) :: zpsi(:, :)
    CMPLX,           intent(in)    :: factor
    CMPLX,           intent(in)    :: mem_s(:, :, :)
    integer,         intent(in)    :: order
    integer,         intent(in)    :: dim
    integer,         intent(in)    :: mapping(:)   ! the mapping

    CMPLX, allocatable :: tmem(:, :)

    call push_sub('td_trans_rti.apply_sp_mem')

    ALLOCATE(tmem(intf%np, intf%np), intf%np**2)
    call make_full_matrix(intf%np, order, dim, sp_mem, mem_s, tmem, mapping)
    call apply_mem(tmem, il, intf, intf_wf, zpsi, factor)
    deallocate(tmem)

    call pop_sub()
  end subroutine apply_sp_mem

end module td_trans_rti_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
