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
  public ::      &
    ob_rti_init, &
    ob_rti_end,  &
    cn_src_mem_dt

  ! We use 1st order Taylor expansion of exp(-i dt/2 H)|psi>
  ! to calculate (1 - i dt H)|psi>.
  type(exponential_t) :: taylor_1st

  ! Parameters to the BiCG in Crank-Nicholson.
  integer :: cg_max_iter
  FLOAT   :: cg_tol

  ! Pointers for the h_eff_backward(t) operator for the iterative linear solver.
  type(hamiltonian_t), pointer :: hm_p
  type(grid_t), pointer        :: gr_p
  CMPLX, pointer               :: mem_p(:, :, :), sp_mem_p(:, :), mem_s_p(:, :, :, :), green_l_p(:, :, :)
  type(interface_t), pointer   :: intf_p(:)
  FLOAT, pointer               :: dt_p, t_p, energy_p
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
    
    integer            :: order, it, allocsize
    CMPLX              :: um(NLEADS)
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

    ! Calculate td-potential.
    ALLOCATE(td_pot(0:max_iter+1, NLEADS), (max_iter+2)*NLEADS)
    call lead_td_pot(td_pot, gr%sb%lead_td_pot_formula, max_iter, dt)
    ! Allocate memory for the src_mem_u (needed for source and memory term.
    ALLOCATE(ob%src_mem_u(0:max_iter, NLEADS), (max_iter+1)*NLEADS)
    !            /      dt        \   / /       dt        \
    ! u(m, il) = |1 - i -- U(m,il)|  /  | 1 + i -- U(m,il) |
    !            \      4         / /   \       4         /
    do it = 0, max_iter
      um(1:NLEADS)               = M_HALF*(td_pot(it+1, 1:NLEADS) + td_pot(it, 1:NLEADS))
      ob%src_mem_u(it, 1:NLEADS) = (M_z1 - M_zI*dt/M_FOUR*um(1:NLEADS)) / (M_z1 + M_zI*dt/M_FOUR*um(1:NLEADS))
    end do
    deallocate(td_pot)

    call ob_rti_write_info(ob, st, gr, max_iter, order)

    ! Initialize source and memory terms.
    call ob_mem_init(gr%intf, hm, ob, dt/M_TWO, max_iter, gr%der%lapl, &
      gr%sb%h(TRANS_DIR), order, st%mpi_grp)
    call ob_src_init(ob, st, gr%intf(LEFT)%np)

    ! Allocate memory for the interface wave functions of previous
    ! timesteps.
    allocsize = gr%intf(LEFT)%np*st%lnst*st%d%kpt%nlocal*(max_iter+1)*NLEADS
    ALLOCATE(ob%st_intface(gr%intf(LEFT)%np, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end, NLEADS, 0:max_iter), allocsize)
    ob%st_intface = M_z0

    call pop_sub()
  end subroutine ob_rti_init


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

    integer            :: il, it, m, cg_iter, j, order, ierr, inp
    integer, target    :: ist, ik, intf_np
    CMPLX              :: factor, alpha, fac, f0
    CMPLX, allocatable :: tmp(:, :), tmp_wf(:), tmp_mem(:, :)
    CMPLX, allocatable :: ext_wf(:, :, :, :) ! (gr%mesh%np+2*np, ndim, nst, nik)
    character(len=100) :: filename
    FLOAT              :: dres
    
    call push_sub('ob_rti.cn_src_mem_dt')

    intf_np = gr%intf(LEFT)%np ! Assuming symmetric leads.

    order = gr%der%order
    ALLOCATE(tmp(gr%mesh%np, st%d%ispin), gr%mesh%np*st%d%ispin)
    ALLOCATE(tmp_wf(intf_np), intf_np)
    if(ob%mem_type.eq.save_cpu_time) then
      ALLOCATE(tmp_mem(intf_np, intf_np), intf_np**2)
    else
      ALLOCATE(tmp_mem(intf_np*order, 1), intf_np*order)
    end if

    ! Set pointers to communicate with with backward propagator passed
    ! to iterative linear solver.
    hm_p       => hm
    gr_p       => gr
    mem_p      => ob%mem_coeff(:, :, 0, :)
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

    m   = timestep-1
    inp = intf_np

    ! Save interface part of wavefunctions for subsequent iterations
    ! before we overwrite them with the values for the new timestep.
    ! (Copying before the propagation gets the saving right for the
    ! initial state also.)
    do il = 1, NLEADS
      call save_intf_wf(gr%intf(il), st, ob%st_intface(:, :, :, il, timestep-1))
    end do

    ! Get right-hand side.
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        ! 1. Apply effective Hamiltonian.
        if(ob%mem_type.eq.SAVE_CPU_TIME) then
          call apply_h_eff(hm, gr, ob%mem_coeff(:, :, 0, :), gr%intf, -M_ONE, dt, t, &
            ist, ik, st%zpsi(:, :, ist, ik))
        else
          call apply_h_eff_sp(hm, gr, ob%mem_sp_coeff(:, 0, :), gr%intf, -M_ONE, dt, t, &
            ist, ik, st%zpsi(:, :, ist, ik), ob%mem_s(:, :, :, :), ob%sp2full_map)
        end if

        do il = 1, NLEADS
          ! 2. Add source term
          if(iand(ob%additional_terms, SRC_TERM_FLAG).ne.0) then
            f0  = M_z1/(M_z1+M_zI*M_HALF*dt*st%ob_eigenval(ist, ik))
            fac = (M_z1-M_zI*M_HALF*dt*st%ob_eigenval(ist, ik))*f0
            if(ob%mem_type.eq.SAVE_CPU_TIME) then
              tmp_mem(:, :) = ob%mem_coeff(:, :, m, il)
              if(m.gt.0) tmp_mem(:, :) = tmp_mem(:, :) + ob%mem_coeff(:, :, m-1, il)
              call calc_source_wf(max_iter, m, inp, il, hm%lead_h_offdiag(:, :, il), tmp_mem, dt, &
                st%ob_intf_psi(:, :, 1, ist, ik, il), ob%src_mem_u(:, il), f0, fac,              &
                lambda(m, 0, max_iter, ob%src_mem_u(:, il)), ob%src_prev(:, 1, ist, ik, il))
            else
              tmp_mem(:, 1) = ob%mem_sp_coeff(:, m, il)
              if(m.gt.0) tmp_mem(:, 1) = tmp_mem(:, 1) + ob%mem_sp_coeff(:, m-1, il)
              call calc_source_wf_sp(max_iter, m, inp, il, hm%lead_h_offdiag(:, :, il),     &
                tmp_mem(:, 1), dt, order, gr%sb%dim, st%ob_intf_psi(:, :, 1, ist, ik, il), &
                ob%mem_s(:, :, :, il), ob%sp2full_map, ob%src_mem_u(:, il), f0, fac,   &
                lambda(m, 0, max_iter, ob%src_mem_u(:, il)), ob%src_prev(:, 1, ist, ik, il))
            end if
            call apply_src(gr%intf(il), ob%src_prev(:, 1, ist, ik, il), st%zpsi(:, :, ist, ik))
          end if
          ! 3. Add memory term.
          if(iand(ob%additional_terms, MEM_TERM_FLAG).ne.0) then
            do it = 0, m-1
              factor = -dt**2/M_FOUR*lambda(m, it, max_iter, ob%src_mem_u(:, il)) / &
                (ob%src_mem_u(m, il)*ob%src_mem_u(it, il))
              tmp_wf(:) = ob%st_intface(:, ist, ik, il, it+1) + ob%st_intface(:, ist, ik, il, it)
              select case(ob%mem_type)
              case(SAVE_CPU_TIME)
                tmp_mem(:, :) = ob%mem_coeff(:, :, m-it, il)
                if((m-it).gt.0) tmp_mem(:, :) = tmp_mem(:, :) + ob%mem_coeff(:, :, m-it-1, il)
                call apply_mem(tmp_mem, gr%intf(il), tmp_wf, st%zpsi(:, :, ist, ik), factor)
              case(SAVE_RAM_USAGE)
                tmp_mem(:, 1) = ob%mem_sp_coeff(:, m-it, il)
                if((m-it).gt.0) tmp_mem(:, 1) = tmp_mem(:, 1) + ob%mem_sp_coeff(:, m-it-1, il)
                call apply_sp_mem(tmp_mem(:, 1), il, gr%intf(il), tmp_wf, st%zpsi(:, :, ist, ik), &
                  factor, ob%mem_s(:, :, :, il), order, gr%sb%dim, ob%sp2full_map)
              end select
            end do
          end if
        end do

        ! 4. Solve linear system (1 + i \delta H_{eff}) st%zpsi = tmp.
        cg_iter      = cg_max_iter
        tmp(1:gr%mesh%np, 1) = st%zpsi(1:gr%mesh%np, 1, ist, ik)
        ! tmp(1:gr%mesh%np, 1) = M_z0
        ! Use QMR-solver since it is faster and more stable than BiCG.
        select case(ob%mem_type)
        case(SAVE_CPU_TIME)
          ! Use for non-magnetic case.
          ! Either the stable symmetric QMR solver
          call zqmr_sym(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(:, 1), h_eff_backward, precond_prop, &
            cg_iter, residue=dres, threshold=cg_tol, showprogress=.false.)
          ! or the slightly faster but less stable CG solver.
          ! call zconjugate_gradients(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(:, 1), h_eff_backward, &
          !   zdot_productu, cg_iter, threshold=cg_tol)

          ! Use for magnetic fields a general solver like BiCG (working)
          ! call zconjugate_gradients(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(:, 1), &
          !   h_eff_backward, h_eff_backward_dagger, zmf_dotp_aux, cg_iter, threshold=cg_tol)
          ! or the general QMR solver (not working yet)
          ! call zqmr(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(:, 1), h_eff_backward, h_eff_backward_t, &
          !   precond_prop, precond_prop, cg_iter, residue=dres, threshold=cg_tol, showprogress=.true.)
        case(SAVE_RAM_USAGE)
          call zqmr_sym(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(:, 1), h_eff_backward_sp, precond_prop, &
            cg_iter, residue=dres, threshold=cg_tol, showprogress=.false.)
          ! call zqmr(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(:, 1), h_eff_backward_sp, h_eff_backward_sp_t, &
          !   precond_prop, precond_prop, cg_iter, residue=dres, threshold=cg_tol, showprogress=.false.)
          ! call zconjugate_gradients(gr%mesh%np, st%zpsi(:, 1, ist, ik), tmp(1:gr%mesh%np, 1), &
          !   h_eff_backward_sp, h_eff_backwardt_sp, zmf_dotp_aux, cg_iter, threshold=cg_tol)
        end select
        ! Write warning if BiCG did not converge.
        ! if(cg_iter.gt.cg_max_iter) then
        !   write(message(1), '(a,i5,a)') 'BiCG did not converge in ', cg_max_iter, ' iterations'
        !   message(2) = 'when propagating with H_eff.'
        !   call write_warning(2)
        ! end if
      end do
    end do

    ! Save interface part of wavefunction for subsequent iterations.
    do il = 1, NLEADS
      call save_intf_wf(gr%intf(il), st, ob%st_intface(:, :, :, il, timestep))
    end do

    deallocate(tmp, tmp_wf, tmp_mem)
    call pop_sub()

  contains

    CMPLX function lambda(m, k, max_iter, sm_u) result(res)
      integer, intent(in) :: m, k, max_iter
      CMPLX,   intent(in) :: sm_u(0:max_iter)

      integer :: j
      res = M_z1

      do j = k, m
        res = res*sm_u(j)**2
      end do
    end function lambda
  end subroutine cn_src_mem_dt


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff}
  ! Used by the iterative linear solver.
  subroutine h_eff_backward(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('ob_rti.h_eff_backward')
    
    ALLOCATE(tmp(gr_p%mesh%np_part, 1),gr_p%mesh%np_part)
    ! Propagate backward.
    call lalg_copy(gr_p%mesh%np, x, tmp(:, 1))
    call apply_h_eff(hm_p, gr_p, mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp)
    call lalg_copy(gr_p%mesh%np, tmp(:, 1), y)

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backward


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff} (transposed)
  ! Used by the iterative linear solver.
  subroutine h_eff_backward_t(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('ob_rti.h_eff_backward_t')
    
    ALLOCATE(tmp(gr_p%mesh%np_part, 1),gr_p%mesh%np_part)
    ! Propagate backward.
    call lalg_copy(gr_p%mesh%np, x, tmp(:, 1))
    call apply_h_eff(hm_p, gr_p, mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp, .true.)
    call lalg_copy(gr_p%mesh%np, tmp(:, 1), y)

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backward_t


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff} (daggered)
  ! Used by the iterative linear solver.
  subroutine h_eff_backward_dagger(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('ob_rti.h_eff_backward_dagger')
    
    ALLOCATE(tmp(gr_p%mesh%np_part, 1),gr_p%mesh%np_part)
    ! Propagate backward.
    call lalg_copy(gr_p%mesh%np, x, tmp(:, 1))
    call apply_h_eff_dagger(hm_p, gr_p, mem_p, intf_p, dt_p, t_p, ist_p, ik_p, tmp)
    call lalg_copy(gr_p%mesh%np, tmp(:, 1), y)

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backward_dagger


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff}
  ! Used by the iterative linear solver.
  subroutine h_eff_backward_sp(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('ob_rti.h_eff_backward_sp')
    
    ALLOCATE(tmp(gr_p%mesh%np_part, 1),gr_p%mesh%np_part)
    ! Propagate backward.
    call lalg_copy(gr_p%mesh%np, x, tmp(:, 1))
    call apply_h_eff_sp(hm_p, gr_p, sp_mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp, mem_s_p, mapping_p)
    call lalg_copy(gr_p%mesh%np, tmp(:, 1), y)

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backward_sp


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff} (transposed)
  ! Used by the iterative linear solver.
  subroutine h_eff_backward_sp_t(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('ob_rti.h_eff_backward_sp_t')
    
    ALLOCATE(tmp(gr_p%mesh%np_part, 1),gr_p%mesh%np_part)
    ! Propagate backward.
    call lalg_copy(gr_p%mesh%np, x, tmp(:, 1))
    call apply_h_eff_sp(hm_p, gr_p, sp_mem_p, intf_p, M_ONE, dt_p, t_p, ist_p, ik_p, tmp, mem_s_p, mapping_p, .true.)
    call lalg_copy(gr_p%mesh%np, tmp(:, 1), y)

    deallocate(tmp)
    call pop_sub()
  end subroutine h_eff_backward_sp_t


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
    type(interface_t),   intent(in)    :: intf(NLEADS)
    FLOAT,               intent(in)    :: sign, dt, t
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)
    logical, optional,   intent(in)    :: transposed

    integer            :: il
    CMPLX, allocatable :: intf_wf(:, :)

    call push_sub('ob_rti.apply_h_eff')

    ALLOCATE(intf_wf(intf(LEFT)%np, NLEADS), intf(LEFT)%np*NLEADS)

    ASSERT(sign.eq.M_ONE.or.sign.eq.-M_ONE)

    call get_intf_wf(intf(LEFT), zpsi(:, 1), intf_wf(:, LEFT))
    call get_intf_wf(intf(RIGHT), zpsi(:, 1), intf_wf(:, RIGHT))

    ! To act with the transposed of H on the wavefunction
    ! we apply H to the conjugate of psi and conjugate the
    ! resulting H psi 
    if(present(transposed)) then
      if(transposed) zpsi = conjg(zpsi)
    end if

    ! Apply (1 - i\delta H_{CC}^{(m)}) to zpsi.
    ! td_exp_dt with Taylor expansion calculates exp(-i dt H), i. e. the
    ! minus is already built in.
    call exponential_apply(taylor_1st, gr, hm, zpsi, ist, ik, -sign*dt/M_TWO, t-dt)
    if(present(transposed)) then
      if(transposed) zpsi = conjg(zpsi)
    end if

    ! Apply modification: sign \delta^2 Q zpsi
    do il = 1, NLEADS
      call apply_mem(mem(:, :, il), intf(il), intf_wf(:, il), zpsi, TOCMPLX(sign*dt**2/M_FOUR, M_ZERO))
    end do

    deallocate(intf_wf)

    call pop_sub()
  end subroutine apply_h_eff


  ! ---------------------------------------------------------
  ! Propagate forward/backward with effective Hamiltonian.
  subroutine apply_h_eff_dagger(hm, gr, mem, intf, dt, t, ist, ik, zpsi)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    CMPLX, target,       intent(in)    :: mem(:, :, :)
    type(interface_t),   intent(in)    :: intf(NLEADS)
    FLOAT,               intent(in)    :: dt, t
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)

    integer            :: il
    CMPLX, allocatable :: intf_wf(:, :)

    call push_sub('ob_rti.apply_h_eff_dagger')

    ALLOCATE(intf_wf(intf(LEFT)%np, NLEADS), intf(LEFT)%np*NLEADS)

    call get_intf_wf(intf(LEFT), zpsi(:, 1), intf_wf(:, LEFT))
    call get_intf_wf(intf(RIGHT), zpsi(:, 1), intf_wf(:, RIGHT))

    ! Apply (1 - i\delta H_{CC}^{(m)}) to zpsi.
    ! td_exp_dt with Taylor expansion calculates exp(-i dt H), i. e. the
    ! minus is already built in.
    call exponential_apply(taylor_1st, gr, hm, zpsi, ist, ik, dt/M_TWO, t-dt)

    ! Apply modification: sign \delta^2 Q zpsi
    do il = 1, NLEADS
      call apply_mem(conjg(mem(:, :, il)), intf(il), intf_wf(:, il), zpsi, TOCMPLX(dt**2/M_FOUR, M_ZERO))
    end do

    deallocate(intf_wf)

    call pop_sub()
  end subroutine apply_h_eff_dagger


  ! ---------------------------------------------------------
  ! Propagate forward/backward with effective Hamiltonian.
  subroutine apply_h_eff_sp(hm, gr, sp_mem, intf, sign, dt, t, ist, ik, zpsi, mem_s, mapping, transposed)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    CMPLX, target,       intent(in)    :: sp_mem(:, :)
    type(interface_t),   intent(in)    :: intf(NLEADS)
    FLOAT,               intent(in)    :: sign, dt, t
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)
    CMPLX,               intent(in)    :: mem_s(:, :, :, :)
    integer,             intent(in)    :: mapping(:)   ! the mapping
    logical, optional,   intent(in)    :: transposed

    integer            :: il
    CMPLX, allocatable :: intf_wf(:, :)

    call push_sub('ob_rti.apply_h_eff_sp')

    ALLOCATE(intf_wf(intf(LEFT)%np, NLEADS), intf(LEFT)%np*NLEADS)

    ASSERT(sign.eq.M_ONE.or.sign.eq.-M_ONE)

    call get_intf_wf(intf(LEFT), zpsi(:, 1), intf_wf(:, LEFT))
    call get_intf_wf(intf(RIGHT), zpsi(:, 1), intf_wf(:, RIGHT))

    if(present(transposed)) then 
      if(transposed) zpsi = conjg(zpsi)
    end if

    ! Apply (1 - i\delta H_{CC}^{(m)}) to zpsi.
    ! td_exp_dt with Taylor expansion calculates exp(-i dt H), i. e. the
    ! minus is already built in.
    call exponential_apply(taylor_1st, gr, hm, zpsi, ist, ik, -sign*dt/M_TWO, t-dt)

    if(present(transposed)) then
      if(transposed) zpsi = conjg(zpsi)
    end if

    ! Apply modification: sign \delta^2 Q zpsi
    do il = 1, NLEADS
      call apply_sp_mem(sp_mem(:, il), il, intf(il), intf_wf(:, il), zpsi, &
        TOCMPLX(sign*dt**2/M_FOUR, M_ZERO), mem_s(:, :, :, il), gr%der%order, gr%sb%dim, mapping)
    end do

    deallocate(intf_wf)

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
      zpsi(intf%index_range(1):intf%index_range(2), 1) + src_wf(:)

    call pop_sub()
  end subroutine apply_src


  ! ---------------------------------------------------------
  ! Apply memory coefficient: zpsi <- zpsi + factor Q intf_wf
  ! Use symmetric matrix-vector multiply, because Q is symmetric
  subroutine apply_mem(mem, intf, intf_wf, zpsi, factor)
    CMPLX,             intent(in)    :: mem(:, :)
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: intf_wf(:)
    CMPLX,             intent(inout) :: zpsi(:, :)
    CMPLX,             intent(in)    :: factor

    CMPLX, allocatable :: mem_intf_wf(:)
    CMPLX, allocatable :: tmem(:, :, :)

    call push_sub('ob_rti.apply_mem')

    ALLOCATE(mem_intf_wf(intf%np), intf%np)
    ALLOCATE(tmem(intf%np, intf%np, 2), 2*intf%np**2)

    mem_intf_wf = M_z0

    call zsymv('U', intf%np, factor, mem, intf%np, intf_wf, 1, M_z0, mem_intf_wf, 1)

    ! Do not use use BLAS here.
    zpsi(intf%index_range(1):intf%index_range(2), 1) = &
      zpsi(intf%index_range(1):intf%index_range(2), 1) + mem_intf_wf(:)

    deallocate(mem_intf_wf, tmem)
    call pop_sub()
  end subroutine apply_mem

  ! ---------------------------------------------------------
  ! Apply memory coefficient: zpsi <- zpsi + factor Q intf_wf
  ! Use symmetric matrix-vector multiply, because Q is symmetric
  subroutine apply_sp_mem(sp_mem, il, intf, intf_wf, zpsi, factor, mem_s, order, dim, mapping)
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

    CMPLX, allocatable :: tmem(:, :)

    call push_sub('ob_rti.apply_sp_mem')

    ALLOCATE(tmem(intf%np, intf%np), intf%np**2)
    call make_full_matrix(intf%np, order, dim, sp_mem, mem_s, tmem, mapping)
    call apply_mem(tmem, intf, intf_wf, zpsi, factor)
    deallocate(tmem)

    call pop_sub()
  end subroutine apply_sp_mem


  ! ---------------------------------------------------------
  subroutine precond_prop(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: diag(:, :)
    integer            :: np

    call push_sub('ob_rti.preconditioner')
    np = gr_p%mesh%np

    y(:) = x(:) ! no preconditioner
!     AL LOCATE(diag(np, 1), np)

!     call zhamiltonian_diagonal(hm_p, gr_p, diag, 1)

!     diag(:, 1) = M_z1 + M_HALF*dt_p*M_zI*diag(:,1) 
!     y(1:np)    = x(1:np)/diag(1:np, 1)

!     deallocate(diag)

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
    
    write(message(1), '(a,a10)') 'Type of memory coefficients:     ', trim(mem_type_name)
    write(message(2), '(a,f10.3)') 'MBytes required for memory term: ', &
      mbytes_memory_term(max_iter, gr%intf(:)%np, NLEADS, st, ob%mem_type, order)
    write(message(3), '(a,i10)') 'Maximum BiCG iterations:         ', cg_max_iter
    write(message(4), '(a,es10.1)') 'BiCG residual tolerance:         ', cg_tol
    write(message(5), '(a,a20)') 'Included additional terms:       ', trim(terms)
    write(message(6), '(a,2a10)') 'TD lead potential (L/R):         ', &
      trim(gr%sb%lead_td_pot_formula(LEFT)), trim(gr%sb%lead_td_pot_formula(RIGHT))
    call write_info(6, stdout)

    call messages_print_stress(stdout)

    call pop_sub()
  end subroutine ob_rti_write_info


  ! ---------------------------------------------------------
  subroutine ob_rti_end(ob)
    type(ob_terms_t), intent(inout) :: ob

    call push_sub('ob_rti.ob_rti_end')

    call ob_mem_end(ob)
    call ob_src_end(ob)

    DEALLOC(ob%src_mem_u)
    DEALLOC(ob%st_intface)

    call pop_sub()
  end subroutine ob_rti_end


end module ob_rti_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
