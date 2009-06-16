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

! Implementation of the propagator for open boundaries, i. e. the
! modified Crank-Nicholson with source and memory terms.

#include "global.h"

module ob_rti_m
  use datasets_m
  use exponential_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_basic_m
  use lasers_m
  use loct_parser_m
  use math_m
  use mesh_function_m
  use messages_m
  use ob_interface_m
  use ob_lead_m
  use ob_mem_m
  use ob_src_m
  use ob_terms_m
  use profiling_m
  use simul_box_m
  use solvers_m
  use states_m
  use varinfo_m
  use v_ks_m

  implicit none

  private
  public ::         &
    ob_rti_init,    &
    ob_rti_end,     &
    cn_src_mem_dt,  &
    cn_src_mem_sp_dt

  ! We use 1st order Taylor expansion of exp(-i dt/2 H)|psi>
  ! to calculate (1 - i dt H)|psi>.
  type(exponential_t) :: taylor_1st

  ! Parameters to the BiCG in Crank-Nicholson.
  integer :: cg_max_iter
  FLOAT   :: cg_tol

  ! is the effective Hamiltonian complex symmetric?
  ! true if no magnetic fields or vector potentials are present
  ! if true some simplifications and speedups are possible
  logical :: heff_sym

  ! Pointers for the h_eff_backward(t) operator for the iterative linear solver.
  type(hamiltonian_t), pointer :: hm_p
  type(grid_t), pointer        :: gr_p
  CMPLX, pointer               :: mem_p(:, :, :), sp_mem_p(:, :), mem_s_p(:, :, :, :)
  type(interface_t), pointer   :: intf_p(:)
  FLOAT, pointer               :: dt_p, t_p
  integer, pointer             :: ist_p, ik_p, mem_type_p, mapping_p(:)
  
contains

  ! ---------------------------------------------------------
  ! Initialize propagator.
  subroutine ob_rti_init(st, gr, hm, ob, dt, max_iter)
    type(states_t),      intent(in)  :: st
    type(grid_t),        intent(in)  :: gr
    type(hamiltonian_t), intent(in)  :: hm
    type(ob_terms_t),    intent(out) :: ob
    FLOAT,               intent(in)  :: dt
    integer,             intent(in)  :: max_iter
    
    integer            :: order, it, np, il
    CMPLX, allocatable :: um(:)
    FLOAT, allocatable :: td_pot(:, :)

    call push_sub('ob_rti.ob_rti_init')

    taylor_1st%exp_method = TAYLOR
    taylor_1st%exp_order  = 1
    order                 = gr%der%order

    !%Variable OpenBoundariesBiCGMaxIter
    !%Type integer
    !%Default 100
    !%Section Open Boundaries
    !%Description
    !% Sets the maximum iteration number for the BiCG linear solver in
    !% the Crank-Nicholson procedure.
    !%End
    call loct_parse_int(datasets_check('OpenBoundariesBiCGMaxIter'), 100, cg_max_iter)
    if(cg_max_iter.le.0) then
      call input_error('OpenBoundariesBiCGMaxIter')
    end if

    !%Variable OpenBoundariesBiCGTol
    !%Type integer
    !%Default 1e-12
    !%Section Open Boundaries
    !%Description
    !% Sets the convergence tolerance for the residue in the BiCG linear solver
    !% in the Crank-Nicholson procedure.
    !%End
    call loct_parse_float(datasets_check('OpenBoundariesBiCGTol'), CNST(1e-12), cg_tol)
    if(cg_tol.le.M_ZERO) then
      call input_error('OpenBoundariesBiCGTol')
    end if

    !%Variable OpenBoundariesMemType
    !%Type integer
    !%Default save_cpu_time
    !%Section Open Boundaries
    !%Description
    !% Decides whether the memory coefficients use lots of RAM (default)
    !% or uses a more compact scheme but with the need of more CPU-cycles.
    !%
    !%Option save_cpu_time 1
    !% Use the memory-intensive procedure
    !%Option save_ram_usage 2
    !% Use the RAM-saving, CPU-intensive procedure
    !%End
    call loct_parse_int(datasets_check('OpenBoundariesMemType'), SAVE_CPU_TIME, ob%mem_type)
    if(.not.varinfo_valid_option('OpenBoundariesMemType', ob%mem_type)) then
      call input_error('OpenBoundariesMemType')
    end if

    !%Variable OpenBoundariesAdditionalTerms
    !%Type flag
    !%Default mem_term + src_term
    !%Section Open Boundaries
    !%Description
    !% The open-boundaries propagator inserts two additional terms in
    !% the Crank-Nicholson scheme: source and memory. With this variable,
    !% one or both of them can be switched off.
    !%
    !%Option mem_term 1
    !% If present, include memory term in propagator
    !%Option src_term 2
    !% If present, include source term in propagator
    !%End
    call loct_parse_int(datasets_check('OpenBoundariesAdditionalTerms'), MEM_TERM_FLAG+SRC_TERM_FLAG, ob%additional_terms)
    if(.not.varinfo_valid_option('OpenBoundariesAdditionalTerms', ob%additional_terms, is_flag=.true.)) then
      call input_error('OpenBoundariesAdditionalTerms')
    end if

    !%Variable OpenBoundariesMaxMemCoeffs
    !%Type integer
    !%Default TDMaximumIter
    !%Section Calculation Modes::Transport
    !%Description
    !% Sets the maximum number of used memory coefficients.
    !% Can only be less than TDMaximumIter if source term is switched off.
    !%End
    call loct_parse_int(datasets_check('OpenBoundariesMaxMemCoeffs'), max_iter, ob%max_mem_coeffs)
    if(ob%max_mem_coeffs.le.0) then
      write(message(1), '(a,i6,a)') "Input : '", ob%max_mem_coeffs, "' is not a valid OpenBoundariesMaxMemCoeffs."
      message(2) = '(0 < OpenBoundariesMaxMemCoeffs <= TDMaximumIter)'
      call write_fatal(2)
    end if
    ob%max_mem_coeffs = min(ob%max_mem_coeffs, max_iter)
    if((iand(ob%additional_terms, SRC_TERM_FLAG).ne.0).and.(ob%max_mem_coeffs.ne.max_iter)) then
      write(message(1), '(a,i6,a)') "Input : '", ob%max_mem_coeffs, "' is not a valid OpenBoundariesMaxMemCoeffs."
      message(2) = 'If in OpenBoundariesAdditionalTerms src_term is present this has to be equal to TDMaximumIter.'
      message(3) = 'If an open system should be simulated the source term should be switched of.'
      call write_fatal(3)
    end if

    ! Calculate td-potential.
    SAFE_ALLOCATE(td_pot(0:max_iter+1, 1:NLEADS))
    call lead_td_pot(td_pot, gr%sb%lead_td_pot_formula, max_iter, dt)
    ! Allocate memory for the src_mem_u (needed for source and memory term.
    SAFE_ALLOCATE(ob%src_mem_u(0:max_iter, 1:NLEADS))
    SAFE_ALLOCATE(um(1:NLEADS))
    !            /      dt        \   / /       dt        \
    ! u(m, il) = |1 - i -- U(m,il)|  /  | 1 + i -- U(m,il) |
    !            \      4         / /   \       4         /
    do it = 0, max_iter
      um(1:NLEADS)               = M_HALF*(td_pot(it+1, 1:NLEADS) + td_pot(it, 1:NLEADS))
      ob%src_mem_u(it, 1:NLEADS) = (M_z1 - M_zI*dt/M_FOUR*um(1:NLEADS)) / (M_z1 + M_zI*dt/M_FOUR*um(1:NLEADS))
    end do
    SAFE_DEALLOCATE_A(um)
    SAFE_DEALLOCATE_A(td_pot)

    call ob_rti_write_info(ob, st, gr, max_iter, order)

    ! Initialize source and memory terms.
    call ob_mem_init(gr%intf, hm, ob, dt/M_TWO, ob%max_mem_coeffs, gr%der%lapl, &
      gr%sb%h(TRANS_DIR), order, st%d%kpt%mpi_grp)
    call ob_src_init(ob, st, gr%intf(LEFT)%np)

    ! Allocate memory for the interface wave functions of previous
    ! timesteps.
    np = gr%intf(LEFT)%np
    SAFE_ALLOCATE(ob%st_intface(1:np, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end, 1:NLEADS, 0:ob%max_mem_coeffs))
    ob%st_intface = M_z0

    ! check the symmetry of the effective hamiltonian
    ! if no magnetic field or vector potential is present then
    ! the h-matrix is complex symmetric, otherwise general complex
    heff_sym = .true.
    if(associated(hm%ep%A_static)) heff_sym = .false.
    do il=1, hm%ep%no_lasers
      if(laser_kind(hm%ep%lasers(il)).eq.E_FIELD_MAGNETIC) heff_sym = .false.
      if(laser_kind(hm%ep%lasers(il)).eq.E_FIELD_VECTOR_POTENTIAL) heff_sym = .false.
    end do

    nullify(mem_p, sp_mem_p, mem_s_p, intf_p, dt_p, t_p, ist_p, ik_p, mem_type_p, mapping_p)

    call pop_sub()
  end subroutine ob_rti_init


  CMPLX function lambda(m, k, max_iter, sm_u) result(res)
    integer, intent(in) :: m, k, max_iter
    CMPLX,   intent(in) :: sm_u(0:max_iter)

    integer :: j
    res = M_z1

    do j = k, m
      res = res*sm_u(j)**2
    end do
  end function lambda

  ! ---------------------------------------------------------
  ! Crank-Nicholson timestep with source and memory.
  ! Only non-interacting electrons for the moment, so no
  ! predictor-corrector scheme.
  subroutine cn_src_mem_dt(ob, st, ks, hm, gr, max_iter, dt, t, timestep)
    type(ob_terms_t), target,    intent(inout) :: ob
    type(states_t),              intent(inout) :: st
    type(v_ks_t),                intent(in)    :: ks
    type(hamiltonian_t), target, intent(inout) :: hm
    type(grid_t), target,        intent(inout) :: gr
    integer,                     intent(in)    :: max_iter
    FLOAT, target,               intent(in)    :: dt
    FLOAT, target,               intent(in)    :: t
    integer,                     intent(in)    :: timestep

    integer            :: il, it, m, cg_iter, order, inp
    integer, target    :: ist, ik
    CMPLX              :: factor, fac, f0
    CMPLX, allocatable :: tmp(:, :), tmp_wf(:), tmp_mem(:, :)
    FLOAT              :: dres
    logical            :: conv
    
    call push_sub('ob_rti.cn_src_mem_dt')

    inp = gr%intf(LEFT)%np ! Assuming symmetric leads.

    order = gr%der%order
    SAFE_ALLOCATE(tmp(1:gr%mesh%np, 1:st%d%ispin))
    SAFE_ALLOCATE(tmp_wf(1:inp))
    SAFE_ALLOCATE(tmp_mem(1:inp, 1:inp))

    ! Set pointers to communicate with with backward propagator passed
    ! to iterative linear solver.
    hm_p       => hm
    gr_p       => gr
    intf_p     => gr%intf
    dt_p       => dt
    t_p        => t
    ist_p      => ist
    ik_p       => ik
    mem_type_p => ob%mem_type
    mem_p      => ob%mem_coeff(:, :, 0, :)

    ! For the dot product passed to BiCG routine.
    call mesh_init_mesh_aux(gr%mesh)

    m = timestep-1

    ! save the initial state
    if (m.eq.0) then
      do il = 1, NLEADS
        call save_intf_wf(gr%intf(il), st, ob%st_intface(1:inp, :, :, il, 0))
      end do
    end if

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        ! Get right-hand side.
        !   1. Apply effective Hamiltonian.
        call apply_h_eff(hm, gr, ob%mem_coeff(:, :, 0, :), gr%intf, -M_ONE, dt, t, &
          ist, ik, st%zpsi(:, :, ist, ik))

        do il = 1, NLEADS
          ! 2. Add source term
          if(iand(ob%additional_terms, SRC_TERM_FLAG).ne.0) then
            f0  = M_z1/(M_z1+M_zI*M_HALF*dt*st%ob_eigenval(ist, ik))
            fac = (M_z1-M_zI*M_HALF*dt*st%ob_eigenval(ist, ik))*f0
            tmp_mem(:, :) = ob%mem_coeff(:, :, m, il)
            if(m.gt.0) tmp_mem(:, :) = tmp_mem(:, :) + ob%mem_coeff(:, :, m-1, il)
            call calc_source_wf(max_iter, m, inp, il, hm%lead_h_offdiag(:, :, il), tmp_mem, dt, &
              st%ob_intf_psi(:, :, 1, ist, ik, il), ob%src_mem_u(:, il), f0, fac,              &
              lambda(m, 0, max_iter, ob%src_mem_u(:, il)), ob%src_prev(:, 1, ist, ik, il))
            call apply_src(gr%intf(il), ob%src_prev(1:inp, 1, ist, ik, il), st%zpsi(:, :, ist, ik))
          end if
          ! 3. Add memory term.
          if(iand(ob%additional_terms, MEM_TERM_FLAG).ne.0) then
            do it = 0, m-1
!            do it = max(m-1), m-1
              factor = -dt**2/M_FOUR*lambda(m, it, max_iter, ob%src_mem_u(:, il)) / &
                (ob%src_mem_u(m, il)*ob%src_mem_u(it, il))
              tmp_wf(:) = ob%st_intface(:, ist, ik, il, it+1) + ob%st_intface(:, ist, ik, il, it)
              tmp_mem(:, :) = ob%mem_coeff(:, :, m-it, il)
              if((m-it).gt.0) tmp_mem(:, :) = tmp_mem(:, :) + ob%mem_coeff(:, :, m-it-1, il)
              call apply_mem(tmp_mem, gr%intf(il), tmp_wf, st%zpsi(:, :, ist, ik), factor)
            end do
          end if
        end do
        if (.not.heff_sym) then
          ! 4. if H_eff is not complex symmetric do the following for solving the linear equation:
          ! A.x = b --> (A^T).A.x = (A^T).b
          ! In this case (A^T).A is complex symmetric and we can use the QMR_sym solver
          call apply_h_eff(hm, gr, ob%mem_coeff(:, :, 0, :), gr%intf, M_ONE, dt, t, &
            ist, ik, st%zpsi(:, :, ist, ik), .true.)
        end if

        ! Solve linear system (1 + i \delta H_{eff}) st%zpsi = tmp.
        cg_iter      = cg_max_iter
        tmp(1:gr%mesh%np, 1) = st%zpsi(1:gr%mesh%np, 1, ist, ik)
        ! Use the stable symmetric QMR solver
        ! h_eff_backward must be a complex symmetric operator !
        call zqmr_sym(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(:, 1), h_eff_backward, precond_prop, &
          cg_iter, residue=dres, threshold=cg_tol, showprogress=.false., converged=conv)
        !if (.not.conv) then
        !  write(*,*) 'ik, residue', ik, dres
        !end if
      end do
    end do

    ! Save interface part of wavefunction for subsequent iterations.
    do il = 1, NLEADS
      call save_intf_wf(gr%intf(il), st, ob%st_intface(:, :, :, il, timestep))
    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_wf)
    SAFE_DEALLOCATE_A(tmp_mem)
    call pop_sub()
  end subroutine cn_src_mem_dt


  ! ---------------------------------------------------------
  ! Crank-Nicholson timestep with source and memory - sparse version.
  ! Only non-interacting electrons for the moment, so no
  ! predictor-corrector scheme.
  subroutine cn_src_mem_sp_dt(ob, st, ks, hm, gr, max_iter, dt, t, timestep)
    type(ob_terms_t), target,    intent(inout) :: ob
    type(states_t),              intent(inout) :: st
    type(v_ks_t),                intent(in)    :: ks
    type(hamiltonian_t), target, intent(inout) :: hm
    type(grid_t), target,        intent(inout) :: gr
    integer,                     intent(in)    :: max_iter
    FLOAT, target,               intent(in)    :: dt
    FLOAT, target,               intent(in)    :: t
    integer,                     intent(in)    :: timestep

    integer            :: il, it, m, cg_iter, order, inp
    integer, target    :: ist, ik
    CMPLX              :: factor, fac, f0
    CMPLX, allocatable :: tmp(:, :), tmp_wf(:), tmp_mem(:)
    FLOAT              :: dres
    
    call push_sub('ob_rti.cn_src_mem_dt')

    inp = gr%intf(LEFT)%np ! Assuming symmetric leads.

    order = gr%der%order
    SAFE_ALLOCATE(tmp(1:gr%mesh%np, 1:st%d%ispin))
    SAFE_ALLOCATE(tmp_wf(1:inp))
    SAFE_ALLOCATE(tmp_mem(1:inp*order))

    ! Set pointers to communicate with with backward propagator passed
    ! to iterative linear solver.
    hm_p       => hm
    gr_p       => gr
    sp_mem_p   => ob%mem_sp_coeff(:, 0, :)
    intf_p     => gr%intf
    dt_p       => dt
    t_p        => t
    ist_p      => ist
    ik_p       => ik
    mem_s_p    => ob%mem_s(:, :, :, :)
    mem_type_p => ob%mem_type
    mapping_p  => ob%sp2full_map

    ! For the dot product passed to BiCG routine.
    call mesh_init_mesh_aux(gr%mesh)

    m = timestep-1


    ! save the initial state
    if (m.eq.0) then
      do il = 1, NLEADS
        call save_intf_wf(gr%intf(il), st, ob%st_intface(1:inp, :, :, il, 0))
      end do
    end if

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        ! Get right-hand side.
        !   1. Apply effective Hamiltonian.
        call apply_h_eff_sp(hm, gr, ob%mem_sp_coeff(:, 0, :), gr%intf, -M_ONE, dt, t, &
          ist, ik, st%zpsi(:, :, ist, ik), ob%mem_s(:, :, :, :), ob%sp2full_map)

        do il = 1, NLEADS
          ! 2. Add source term
          if(iand(ob%additional_terms, SRC_TERM_FLAG).ne.0) then
            f0  = M_z1/(M_z1+M_zI*M_HALF*dt*st%ob_eigenval(ist, ik))
            fac = (M_z1-M_zI*M_HALF*dt*st%ob_eigenval(ist, ik))*f0
            tmp_mem(:) = ob%mem_sp_coeff(:, m, il)
            if(m.gt.0) tmp_mem(:) = tmp_mem(:) + ob%mem_sp_coeff(:, m-1, il)
            call calc_source_wf_sp(max_iter, m, inp, il, hm%lead_h_offdiag(:, :, il),     &
              tmp_mem(:), dt, order, gr%sb%dim, st%ob_intf_psi(:, :, 1, ist, ik, il), &
              ob%mem_s(:, :, :, il), ob%sp2full_map, ob%src_mem_u(:, il), f0, fac,   &
              lambda(m, 0, max_iter, ob%src_mem_u(:, il)), ob%src_prev(:, 1, ist, ik, il))
            call apply_src(gr%intf(il), ob%src_prev(1:inp, 1, ist, ik, il), st%zpsi(:, :, ist, ik))
          end if
          ! 3. Add memory term.
          if(iand(ob%additional_terms, MEM_TERM_FLAG).ne.0) then
            do it = 0, m-1
              factor = -dt**2/M_FOUR*lambda(m, it, max_iter, ob%src_mem_u(:, il)) / &
                (ob%src_mem_u(m, il)*ob%src_mem_u(it, il))
              tmp_wf(:) = ob%st_intface(:, ist, ik, il, it+1) + ob%st_intface(:, ist, ik, il, it)
              tmp_mem(:) = ob%mem_sp_coeff(:, m-it, il)
              if((m-it).gt.0) tmp_mem(:) = tmp_mem(:) + ob%mem_sp_coeff(:, m-it-1, il)
              call apply_sp_mem(tmp_mem(:), il, gr%intf(il), tmp_wf, st%zpsi(:, :, ist, ik), &
                factor, ob%mem_s(:, :, :, il), order, gr%sb%dim, ob%sp2full_map)
            end do
          end if
        end do
        if (.not.heff_sym) then
          ! 4. if H_eff is not complex symmetric do the following for solving the linear equation:
          ! A.x = b --> (A^T).A.x = (A^T).b
          ! In this case (A^T).A is complex symmetric and we can use the QMR_sym solver
          call apply_h_eff_sp(hm, gr, ob%mem_sp_coeff(:, 0, :), gr%intf, M_ONE, dt, t, &
            ist, ik, st%zpsi(:, :, ist, ik), ob%mem_s(:, :, :, :), ob%sp2full_map, .true.)
        end if

        ! Solve linear system (1 + i \delta H_{eff}) st%zpsi = tmp.
        cg_iter      = cg_max_iter
        tmp(1:gr%mesh%np, 1) = st%zpsi(1:gr%mesh%np, 1, ist, ik)
        call zqmr_sym(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(:, 1), h_eff_backward_sp, precond_prop, &
          cg_iter, residue=dres, threshold=cg_tol, showprogress=.false.)
      end do
    end do

    ! Save interface part of wavefunction for subsequent iterations.
    do il = 1, NLEADS
      call save_intf_wf(gr%intf(il), st, ob%st_intface(1:inp, :, :, il, timestep))
    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_wf)
    SAFE_DEALLOCATE_A(tmp_mem)
    call pop_sub()
  end subroutine cn_src_mem_sp_dt


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff}
  ! Used by the iterative linear solver.
  subroutine h_eff_backward(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('ob_rti.h_eff_backward')
    
    SAFE_ALLOCATE(tmp(1:gr_p%mesh%np_part, 1:1))
    ! Propagate backward.
    tmp(1:gr_p%mesh%np, 1) = x(1:gr_p%mesh%np)
    call apply_h_eff(hm_p, gr_p, mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp)
    if (.not.heff_sym) then ! make operator complex symmetric
      call apply_h_eff(hm_p, gr_p, mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp, .true.)
    end if
    y(1:gr_p%mesh%np) = tmp(1:gr_p%mesh%np, 1)

    SAFE_DEALLOCATE_A(tmp)
    call pop_sub()
  end subroutine h_eff_backward


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff}
  ! Used by the iterative linear solver.
  subroutine h_eff_backward_sp(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('ob_rti.h_eff_backward_sp')
    
    SAFE_ALLOCATE(tmp(1:gr_p%mesh%np_part, 1:1))
    ! Propagate backward.
    tmp(1:gr_p%mesh%np, 1) = x(1:gr_p%mesh%np)
    call apply_h_eff_sp(hm_p, gr_p, sp_mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp, mem_s_p, mapping_p)
    if (.not.heff_sym) then ! make operator complex symmetric
      call apply_h_eff_sp(hm_p, gr_p, sp_mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp, mem_s_p, mapping_p, .true.)
    end if
    y(1:gr_p%mesh%np) = tmp(1:gr_p%mesh%np, 1)

    SAFE_DEALLOCATE_A(tmp)
    call pop_sub()
  end subroutine h_eff_backward_sp


  ! ---------------------------------------------------------
  ! Save the interface part of all states st for timestep into
  ! the st_intf array.
  subroutine save_intf_wf(intf, st, st_intf)
    type(interface_t), intent(in)    :: intf
    type(states_t),    intent(in)    :: st
    CMPLX,             intent(inout) :: st_intf(1:intf%np, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end)

    integer :: ik, ist
    
    call push_sub('ob_rti.save_intf_wf')

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        call get_intf_wf(intf, st%zpsi(:, 1, ist, ik), st_intf(:, ist, ik))
      end do
    end do

    call pop_sub()
  end subroutine save_intf_wf


  ! ---------------------------------------------------------
  ! Propagate forward/backward with effective Hamiltonian.
  subroutine apply_h_eff(hm, gr, mem, intf, sign, dt, t, ist, ik, zpsi, transposed)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    CMPLX, target,       intent(in)    :: mem(:, :, :)
    type(interface_t),   intent(in)    :: intf(1:NLEADS)
    FLOAT,               intent(in)    :: sign, dt, t
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)
    logical, optional,   intent(in)    :: transposed

    integer            :: il
    CMPLX, allocatable :: intf_wf(:, :)
    logical            :: transposed_

    call push_sub('ob_rti.apply_h_eff')

    SAFE_ALLOCATE(intf_wf(1:intf(LEFT)%np, 1:NLEADS))

    ASSERT(sign.eq.M_ONE.or.sign.eq.-M_ONE)

    transposed_ = .false.
    if(present(transposed)) then
      transposed_ = transposed
    end if


    call get_intf_wf(intf(LEFT), zpsi(:, 1), intf_wf(:, LEFT))
    call get_intf_wf(intf(RIGHT), zpsi(:, 1), intf_wf(:, RIGHT))

    ! To act with the transposed of H on the wavefunction
    ! we apply H to the conjugate of psi and conjugate the
    ! resulting H psi 
    if(transposed_) zpsi = conjg(zpsi)

    ! Apply (1 - i\delta H_{CC}^{(m)}) to zpsi.
    ! td_exp_dt with Taylor expansion calculates exp(-i dt H), i. e. the
    ! minus is already built in.
    call exponential_apply(taylor_1st, gr, hm, zpsi, ist, ik, -sign*dt/M_TWO, t-dt)
    if(transposed_) zpsi = conjg(zpsi)

    ! Apply modification: sign \delta^2 Q zpsi
    do il = 1, NLEADS
      call apply_mem(mem(:, :, il), intf(il), intf_wf(:, il), zpsi, &
                     TOCMPLX(sign*dt**2/M_FOUR, M_ZERO), transposed_)
    end do

    SAFE_DEALLOCATE_A(intf_wf)

    call pop_sub()
  end subroutine apply_h_eff


  ! ---------------------------------------------------------
  ! Propagate forward/backward with effective Hamiltonian.
  subroutine apply_h_eff_sp(hm, gr, sp_mem, intf, sign, dt, t, ist, ik, zpsi, mem_s, mapping, transposed)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    CMPLX, target,       intent(in)    :: sp_mem(:, :)
    type(interface_t),   intent(in)    :: intf(1:NLEADS)
    FLOAT,               intent(in)    :: sign, dt, t
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)
    CMPLX,               intent(in)    :: mem_s(:, :, :, :)
    integer,             intent(in)    :: mapping(:)   ! the mapping
    logical, optional,   intent(in)    :: transposed

    integer            :: il
    CMPLX, allocatable :: intf_wf(:, :)
    logical            :: transposed_

    call push_sub('ob_rti.apply_h_eff_sp')

    SAFE_ALLOCATE(intf_wf(1:intf(LEFT)%np, 1:NLEADS))

    ASSERT(sign.eq.M_ONE.or.sign.eq.-M_ONE)

    transposed_ = .false.
    if(present(transposed)) then
      transposed_ = transposed
    end if

    call get_intf_wf(intf(LEFT), zpsi(:, 1), intf_wf(:, LEFT))
    call get_intf_wf(intf(RIGHT), zpsi(:, 1), intf_wf(:, RIGHT))

    if(transposed_) zpsi = conjg(zpsi)

    ! Apply (1 - i\delta H_{CC}^{(m)}) to zpsi.
    ! td_exp_dt with Taylor expansion calculates exp(-i dt H), i. e. the
    ! minus is already built in.
    call exponential_apply(taylor_1st, gr, hm, zpsi, ist, ik, -sign*dt/M_TWO, t-dt)

    if(transposed_) zpsi = conjg(zpsi)

    ! Apply modification: sign \delta^2 Q zpsi
    do il = 1, NLEADS
      call apply_sp_mem(sp_mem(:, il), il, intf(il), intf_wf(:, il), zpsi, &
        TOCMPLX(sign*dt**2/M_FOUR, M_ZERO), mem_s(:, :, :, il), gr%der%order, gr%sb%dim, mapping)
    end do

    SAFE_DEALLOCATE_A(intf_wf)

    call pop_sub()
  end subroutine apply_h_eff_sp


  ! ---------------------------------------------------------
  ! Apply source coefficient: zpsi <- zpsi + src_wf
  subroutine apply_src(intf, src_wf, zpsi)
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: src_wf(:)
    CMPLX,             intent(inout) :: zpsi(:, :)

    call push_sub('ob_rti.apply_src')

    ! Do not use use BLAS here.
    zpsi(intf%index_range(1):intf%index_range(2), 1) = &
      zpsi(intf%index_range(1):intf%index_range(2), 1) + src_wf(1:intf%np)

    call pop_sub()
  end subroutine apply_src


  ! ---------------------------------------------------------
  ! Apply memory coefficient: zpsi <- zpsi + factor Q intf_wf
  ! Use symmetric matrix-vector multiply, because Q is symmetric
  subroutine apply_mem(mem, intf, intf_wf, zpsi, factor, transposed)
    CMPLX,             intent(in)    :: mem(:, :)
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: intf_wf(:)
    CMPLX,             intent(inout) :: zpsi(:, :)
    CMPLX,             intent(in)    :: factor
    logical, optional, intent(in)    :: transposed

    CMPLX, allocatable :: mem_intf_wf(:)

    call push_sub('ob_rti.apply_mem')

    SAFE_ALLOCATE(mem_intf_wf(1:intf%np))

    mem_intf_wf = M_z0

    if(present(transposed)) then
      ! FIXME: transpose if mem is not symmetric
    end if

    call zsymv('U', intf%np, factor, mem, intf%np, intf_wf, 1, M_z0, mem_intf_wf, 1)

    ! Do not use use BLAS here.
    zpsi(intf%index_range(1):intf%index_range(2), 1) = &
      zpsi(intf%index_range(1):intf%index_range(2), 1) + mem_intf_wf(:)

    SAFE_DEALLOCATE_A(mem_intf_wf)
    call pop_sub()
  end subroutine apply_mem

  ! ---------------------------------------------------------
  ! Apply memory coefficient: zpsi <- zpsi + factor Q intf_wf
  ! Use symmetric matrix-vector multiply, because Q is symmetric
  subroutine apply_sp_mem(sp_mem, il, intf, intf_wf, zpsi, factor, mem_s, order, dim, mapping, transposed)
    CMPLX,             intent(in)    :: sp_mem(:)
    integer,           intent(in)    :: il
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: intf_wf(:)
    CMPLX,             intent(inout) :: zpsi(:, :)
    CMPLX,             intent(in)    :: factor
    CMPLX,             intent(in)    :: mem_s(:, :, :)
    integer,           intent(in)    :: order
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: mapping(:)   ! the mapping
    logical, optional, intent(in)    :: transposed

    CMPLX, allocatable :: tmem(:, :)

    call push_sub('ob_rti.apply_sp_mem')

    SAFE_ALLOCATE(tmem(1:intf%np, 1:intf%np))
    call make_full_matrix(intf%np, order, dim, sp_mem, mem_s, tmem, mapping)
    call apply_mem(tmem, intf, intf_wf, zpsi, factor, transposed)
    SAFE_DEALLOCATE_A(tmem)

    call pop_sub()
  end subroutine apply_sp_mem


  ! ---------------------------------------------------------
  subroutine precond_prop(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

!    integer            :: np
!    CMPLX, allocatable :: diag(:, :)

    call push_sub('ob_rti.preconditioner')
!    np = gr_p%mesh%np

    y(:) = x(:) ! no preconditioner
!    SA FE_ALLOCATE(diag(1:np, 1:1))

!    call zhamiltonian_diagonal(hm_p, gr_p, diag, 1)

!    diag(:, 1) = M_z1 + M_HALF*dt_p*M_zI*diag(:,1) 
!    if (.not.heff_sym) diag(:, 1) = diag(:, 1)**2
!    y(1:np)    = x(1:np)/diag(1:np, 1)

!    SA FE_DEALLOCATE_A(diag)

    call pop_sub()
  end subroutine precond_prop


  ! ---------------------------------------------------------
  ! Write some status information to stdout.
  subroutine ob_rti_write_info(ob, st, gr, max_iter, order)
    type(ob_terms_t), intent(in) :: ob
    type(states_t),   intent(in) :: st
    type(grid_t),     intent(in) :: gr
    integer,          intent(in) :: max_iter
    integer,          intent(in) :: order

    character(len=64) :: terms, mem_type_name

    call push_sub('ob_rti.ob_rti_write_info')

    call messages_print_stress(stdout, 'Open Boundaries')

    if(ob%mem_type.eq.save_cpu_time) then
      mem_type_name = 'full'
    end if
    if(ob%mem_type.eq.save_ram_usage) then
      mem_type_name = 'sparse'
    end if

    terms = ''
    if(iand(ob%additional_terms, mem_term_flag).ne.0) then
      terms = trim(terms)//'memory'
    end if
    if(iand(ob%additional_terms, src_term_flag).ne.0) then
      terms = trim(terms)//' source'
    end if
    
    write(message(1), '(a,a10)')    'Type of memory coefficients:     ', trim(mem_type_name)
    write(message(2), '(a,f10.3)')  'MBytes required for memory term: ', &
      mbytes_memory_term(ob%max_mem_coeffs, gr%intf(:)%np, NLEADS, st, ob%mem_type, order)
    write(message(3), '(a,i10)')    'Maximum BiCG iterations:         ', cg_max_iter
    write(message(4), '(a,es10.1)') 'BiCG residual tolerance:         ', cg_tol
    write(message(5), '(a,a20)')    'Included additional terms:       ', trim(terms)
    write(message(6), '(a,a10)')    'TD left lead potential:          ', &
      trim(gr%sb%lead_td_pot_formula(LEFT))
    write(message(7), '(a,a10)')    'TD right lead potential:         ', &
      trim(gr%sb%lead_td_pot_formula(RIGHT))
    call write_info(7, stdout)

    call messages_print_stress(stdout)

    call pop_sub()
  end subroutine ob_rti_write_info


  ! ---------------------------------------------------------
  subroutine ob_rti_end(ob)
    type(ob_terms_t), intent(inout) :: ob

    call push_sub('ob_rti.ob_rti_end')

    call ob_mem_end(ob)
    call ob_src_end(ob)

    SAFE_DEALLOCATE_P(ob%src_mem_u)
    SAFE_DEALLOCATE_P(ob%st_intface)

    call pop_sub()
  end subroutine ob_rti_end


end module ob_rti_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
