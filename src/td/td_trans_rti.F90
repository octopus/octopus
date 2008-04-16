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
  use states_m
  use td_exp_m
  use td_trans_mem_m
  use td_trans_lead_m
  use td_trans_src_m
  use td_trans_intf_m
  use v_ks_m
  use nl_operator_m
  use loct_m

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
  CMPLX, pointer               :: mem_p(:, :, :), sp_mem_p(:, :), mem_s_p(:, :, :, :), green_l_p(:, :, :)
  type(intface_t), pointer     :: intf_p
  FLOAT, pointer               :: dt_p, t_p, energy_p
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
    FLOAT            :: energy                !
    CMPLX, pointer   :: src_prev(:, :, :, :, :)! intf%np, ndim, nst, nik, NLEADS
    CMPLX, pointer   :: diag(:, :, :)         ! h
    CMPLX, pointer   :: offdiag(:, :, :)      ! hopping operator V^T (np,np,nleads)
    CMPLX, pointer   :: mem_s(:, :, :, :)     ! the matrices to diagonalize coeff0

    type(intface_t)  :: intface

    CMPLX, pointer   :: st_intface(:, :, :, :, :) ! np_intf, nst, nik, nleads, max_iter
    CMPLX, pointer   :: st_psi0(:, :, :, :, :, :) ! np, (lead=1;center=2), ndim, nst, nik, nleads

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
    type(nl_operator_t)  :: op
    FLOAT, allocatable   :: um(:)

    ALLOCATE(um(NLEADS), NLEADS)

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
    call loct_parse_float(check_inp('TDTransBiCGTol'), CNST(1e-12), cg_tol)
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

    !%Variable TDTransEnergy
    !%Type float
    !%Default 0.1
    !%Section Transport
    !%Description
    !% The energy for the wave function.
    !%End
    call loct_parse_float(check_inp('TDTransEnergy'), CNST(0.1), trans%energy)

    call lead_init(trans%lead, max_iter, dt)

    call td_trans_rti_write_info(trans, st, max_iter, gr%f_der%der_discr%order, trans%mem_type)

    call memory_init(trans%intface, dt/M_TWO, max_iter, gr%f_der%der_discr%lapl,&
                      trans%mem_coeff, trans%mem_sp_coeff, trans%mem_s, trans%diag, trans%offdiag, gr%sb%dim, &
                      gr%sb%h(1), trans%mem_type, gr%f_der%der_discr%order, trans%sp2full_map)

    call source_init(st, trans%src_prev, trans%st_psi0, dt, trans%intface%np, gr)

    ! Allocate memory for the interface wave functions.
    allocsize = trans%intface%np*st%lnst*st%d%nik*(max_iter+1)*NLEADS
    ALLOCATE(trans%st_intface(trans%intface%np, st%st_start:st%st_end, st%d%nik, NLEADS, 0:max_iter), allocsize)
    trans%st_intface = M_z0
    
    ! Allocate memory for the src_mem_u (needed for source and memory term)
    ALLOCATE(trans%src_mem_u(0:max_iter, NLEADS), (max_iter+1)*NLEADS)
    !            /      dt        \   / /       dt        \
    ! u(m, il) = |1 - i -- U(m,il)|  /  | 1 + i -- U(m,il) |
    !            \      4         / /   \       4         /
    do id=0, max_iter
      um(:) = M_HALF*(trans%lead%td_pot(id+1, :) + trans%lead%td_pot(id, :))
      !um(:) = trans%lead%td_pot(id, :)
      trans%src_mem_u(id, :) = (M_z1 - M_zI*dt/M_FOUR*um(:)) / (M_z1 + M_zI*dt/M_FOUR*um(:))
    end do
    deallocate(um)
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
    call source_end(trans%src_prev, trans%st_psi0)
    call memory_end(trans%mem_coeff, trans%mem_sp_coeff, trans%diag, trans%offdiag, trans%mem_s, trans%sp2full_map)

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

    call push_sub('td_trans_rti.td_transport_write_info')

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
  ! compute the semi-infinite surface green function
  ! algorithm taken from the paper
  ! Highly convergent schemes for the calculation of bulk and surface Green function
  ! M P Lopez Sanco, J M Sancho and J Rubio (1984)
  ! J. Phys. F: Met. Phys 15 (1985) 851-858
  subroutine green_lead(energy, diag, offdiag, np, il, green, dx)
    FLOAT,               intent(in)  :: energy
    CMPLX,               intent(in)  :: diag(:, :)
    CMPLX,               intent(in)  :: offdiag(:, :)
    integer,             intent(in)  :: np ! number of interface points
    integer,             intent(in)  :: il
    CMPLX,               intent(out) :: green(:, :)
    FLOAT,               intent(out) :: dx

    CMPLX, allocatable   :: e(:, :), es(:, :), a(:, :), b(:, :), inv(:, :)
    CMPLX, allocatable   :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
    integer              :: i, j, ilog_thr, ilog_res
    FLOAT                :: det, old_norm, norm, dx2, threshold, log_thr, log_res, res

    call push_sub('td_trans_rti.green_lead2')

    ALLOCATE(e(np,np),np**2)
    ALLOCATE(es(np,np),np**2)
    ALLOCATE(a(np,np),np**2)
    ALLOCATE(b(np,np),np**2)
    ALLOCATE(inv(np,np),np**2)
    ALLOCATE(tmp1(np,np),np**2)
    ALLOCATE(tmp2(np,np),np**2)
    ALLOCATE(tmp3(np,np),np**2)
    dx2 = dx**2
    threshold = CNST(1e-12)

    write(message(1), '(a,es10.3)') '  Calculating surface Green function for '//trim(lead_name(il))// &
    ' lead: Energy = ', energy
    call write_info(1, stdout)

    ! initialize progress bar
    log_thr = -log(threshold)
    ilog_thr = M_TEN**2*log_thr
    call loct_progress_bar(-1, ilog_thr)

    ! fill with start values
    e(:, :) = diag(:, :)
    a(:, :) = offdiag(:, :)
    do i=1, np
      b(i,1:np) = offdiag(1:np,i)
    end do
    es(:, :) = diag(:, :)
    old_norm = M_ZERO

    ! start the iteration
    do i=1, 1000! 2^1000 efective layers
      inv(:, :) = -e(:, :)
      do j=1, np
        inv(j, j) = inv(j, j) + energy + CNST(1e-10)*M_zI
      end do
      det = lalg_inverter(np, inv, invert = .true.)
      call zgemm('N','N',np,np,np,M_z1,a,np,inv,np,M_z0,tmp2,np)
      call zgemm('N','N',np,np,np,M_z1,tmp2,np,b,np,M_z0,tmp1,np)
      es(:, :) = es(:, :) + tmp1(:, :)
      norm = infinity_norm(es)*dx2
      res = abs(norm-old_norm)
      log_res = -log(res)
      if (log_res > log_thr) log_res = log_thr
      ilog_res = M_TEN**2*log_res
      call loct_progress_bar(ilog_res, ilog_thr)
      if(res.lt.threshold) then
        !write(*,*) 'green-iter:', i
        exit
      end if
      old_norm = norm
      e(:, :) = e(:, :) + tmp1(:, :)
      call zgemm('N','N',np,np,np,M_z1,b,np,inv,np,M_z0,tmp1,np)
      call zgemm('N','N',np,np,np,M_z1,tmp1,np,a,np,M_z1,e,np)
      tmp3(:, :) = a(:, :)
      call zgemm('N','N',np,np,np,M_z1,tmp2,np,tmp3,np,M_z0,a,np)
      tmp3(:, :) = b(:, :)
      call zgemm('N','N',np,np,np,M_z1,tmp1,np,tmp3,np,M_z0,b,np)
    end do
    green(:, :) = -es(:, :)
    do j=1, np
      green(j, j) = green(j, j) + energy + threshold*M_zI
    end do
    det = lalg_inverter(np, green, invert = .true.)
    call make_symmetric_average(green, np)
    det = aimag(green(1,1))
    do i=2,np
      det = det + aimag(green(i,i))
    end do
    if (det.gt.M_ZERO) then
      green = conjg(green)
      det = -det
    end if
    !write(*,*) 'dos', -det

    deallocate(e, es, a, b, inv, tmp1, tmp2, tmp3)
    call pop_sub()
  end subroutine green_lead

  ! ---------------------------------------------------------
  ! the left side of the Lippmann-Schwinger equation
  ! e-H_cc+V(lead)-sum(a)[H_ca*g_a*H_ac]
  ! Used by the iterative linear solver.
  subroutine h_eff_lip_sch(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmpx(:, :)
    CMPLX, allocatable :: tmpy(:, :)
    integer            :: inp, np

    call push_sub('td_trans_rti.h_eff_lip_sch')

    ALLOCATE(tmpx(gr_p%m%np_part, 1), gr_p%m%np_part)
    ALLOCATE(tmpy(gr_p%m%np_part, 1), gr_p%m%np_part)
    inp = intf_p%np
    np = gr_p%m%np

    call lalg_copy(np, x, tmpx(:, 1))
    ! calculate right hand side (e-T-V(lead)-sum(a)[H_ca*g_a*H_ac]
    call zhpsi(h_p, gr_p, tmpx, tmpy, 1, 1)

    y(1:np) = energy_p*x(1:np) - tmpy(1:np,1)
    ! TODO: the static potential of the lead
    call zsymv('U', inp, -M_z1, green_l_p(:, :, LEFT), inp, x(1:inp), 1, M_z1, y(1:inp), 1)
    call zsymv('U', inp, -M_z1, green_l_p(:, :, RIGHT), inp, x(np-inp+1:np), &
               1, M_z1, y(np-inp+1:np), 1)


    deallocate(tmpx, tmpy)
    call pop_sub()
  end subroutine h_eff_lip_sch

  ! ---------------------------------------------------------
  subroutine preconditioner(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: diag(:, :)
    integer            :: np, i

    call push_sub('td_rti.preconditioner')
    np = gr_p%m%np

    ALLOCATE(diag(np, 1), np)

    call zhpsi_diag(h_p, gr_p, diag, 1)

    diag(:, 1) = energy_p - diag(:,1) 
    y(1:np) = x(1:np)/diag(1:np, 1)
  !  y(:) = x(:) ! no preconditioner

    deallocate(diag)

    call pop_sub()
  end subroutine preconditioner

  ! ---------------------------------------------------------
  ! compute the extended eigenstate(s) 
  subroutine ext_eigenstate_lip_sch(h, gr, np, diag, offdiag, order, energy, st)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),   intent(inout)      :: gr
    integer,        intent(in)      :: np                      ! number of interface points
    CMPLX,          intent(in)      :: diag(np, np, NLEADS)    ! lead hamiltonian (diagonal part)
    CMPLX,          intent(in)      :: offdiag(np, np, NLEADS) ! hopping operator V^T
    integer,        intent(in)      :: order                   ! discretization order
    FLOAT, target,  intent(in)      :: energy                  ! total energy
    type(states_t), intent(inout)   :: st                      ! states

    CMPLX, target, allocatable :: green_l(:,:,:)
    CMPLX, allocatable :: tmp(:, :)
    integer            :: i, id, iter, n(3), ist, nst2
    integer, pointer   :: lxyz(:,:)
    FLOAT, target      :: en
    FLOAT              :: lsize(3), q(3), dres, emax, qmin, qmax, dx(3)

    call push_sub('td_trans_rti.ext_eigenstate_lip_sch')

    ALLOCATE(green_l(np, np, NLEADS),np**2*NLEADS)
    ALLOCATE(tmp(NP_PART, 1), NP_PART)

    ! now generate the extended eigenstate with separated x-wavefunction
    ! yet only for a box psi=psi(x)*psi(y)*psi(z)
    ! e_tot = e_x + e_y + e_z, we need e_x for q (for box-sized leads)
    ! FIXME: (probably) run trough all possible linear combinations
    ! yet only the highest state of the separated function is multiplied

    ! 1. compute the unpertubed eigenstates (groundstate calculation)
    ! for now set manually psi(x)*phi(y,z)
    energy_p => en
    lxyz => gr%m%lxyz
    lsize = gr%sb%lsize
    dx(:) = gr%sb%h(:)
    emax = energy
    do id=2, gr%sb%dim ! FIXME: when the last gridpoint is not the border, recalculate transversal energy
      q(id) = M_PI/(M_TWO*(lsize(id)+dx(id)))
      if (order.eq.1) then ! tight binding
        emax = emax - (M_ONE-cos(dx(id)*q(id)))/dx(id)**2
      else ! continuous relation
        emax = emax - M_HALF*q(id)**2
      end if
    end do
    if(emax < M_ZERO) then
      write(message(1), '(a,f14.6,a)') "The input energy : '", energy, &
                  "' is smaller than the groundstate energy for the given system."
      call write_fatal(1)
    end if
    qmin = CNST(1e-5)
    if (order.eq.1) then ! tight binding
      qmax = acos(M_ONE-emax*dx(1)**2)/dx(1)
    else ! continuous relation
      qmax = sqrt(M_TWO*emax)
    end if

    nst2 = (st%lnst-1)/2
    if (nst2.eq.0) nst2 = 1
    ! sample the k integral (k=q), start with highest k
    do ist=st%st_start, st%st_end
      ! always two states (ist=(1,2); (3,4), ...) have the same energy
      ! start with en > emin
      i = (ist-1)/2
      q(1) = qmax - i*(qmax-qmin)/nst2
      ! now calculate the q's for each state
      ! TODO: run over all states, now only the lowest transversal mode is used
      if (order.eq.1) then ! tight binding
        en = (M_ONE-cos(dx(1)*q(1)))/dx(1)**2
      else ! continuous relation
        en = q(1)**2/M_TWO
      end if
      ! now set the unscattered states (flat leads with box geometry)
      do i=1, NP
        if (mod(ist-1,2).eq.0) then
          st%zpsi(i, 1, ist, 1) = exp(M_zI*q(1)*lxyz(i,1)*dx(1))
        else
          st%zpsi(i, 1, ist, 1) = exp(-M_zI*q(1)*lxyz(i,1)*dx(1))
        end if
        do id=2, gr%sb%dim
          if (mod(n(id),2).eq.0) then
            st%zpsi(i, 1, ist, 1) = st%zpsi(i, 1, ist, 1)*sin(q(id)*lxyz(i,id)*dx(id))
          else
            st%zpsi(i, 1, ist, 1) = st%zpsi(i, 1, ist, 1)*cos(q(id)*lxyz(i,id)*dx(id))
          end if
        end do
      end do
      ! calculate right hand side (e-T-V(lead)-sum(a)[H_ca*g_a*H_ac]
      ! only calculate once for a given energy
      if (mod(ist-1,2).eq.0) then
        call green_lead(en, diag(:, :, LEFT), offdiag(:, :, LEFT), np, LEFT, green_l(:, :, LEFT), dx(1))
        call apply_coupling(green_l(:, :, LEFT), offdiag(:, :, LEFT), green_l(:, :, LEFT), np, LEFT)
        call green_lead(en, diag(:, :, RIGHT), offdiag(:, :, RIGHT), np, RIGHT, green_l(:, :, RIGHT), dx(1))
        call apply_coupling(green_l(:, :, RIGHT), offdiag(:, :, RIGHT), green_l(:, :, RIGHT), np, RIGHT)
      end if

      tmp(:,:) = M_z0
      call zkinetic(h, gr, st%zpsi(:, :, ist, 1), tmp(:,:))
      tmp(1:NP, :) = en*st%zpsi(1:NP, :, ist, 1) - tmp(1:NP, :)
      ! TODO: the static potential of the lead
      call zsymv('U', np, -M_z1, green_l(:, :, LEFT), np, st%zpsi(1:np, 1, ist, 1), 1, M_z1, tmp(1:np, 1), 1)
      call zsymv('U', np, -M_z1, green_l(:, :, RIGHT), np, &
        st%zpsi(NP - np + 1:NP, 1, ist, 1), 1, M_z1, tmp(NP - np + 1:NP, 1), 1)

      ! now solve the equation
      green_l_p => green_l
      ! now solve the linear system to get the extended eigenstate
      write(message(1), '(a,es10.3)') '  Solving Lippmann Schwinger equation for energy ', energy
      call write_info(1, stdout)
      iter = 10000
      ! zconjugate_gradients fails in some cases (wrong solution) and takes longer
      ! therefore take the symmetric quasi-minimal residual solver (QMR)
      call zqmr_sym(NP, st%zpsi(:, 1, ist, 1), tmp(:, 1), h_eff_lip_sch, preconditioner, &
                      iter, residue=dres, threshold=cg_tol, showprogress = .true.)
      !write(*,*) 'iter =',iter, 'residue =', dres
    end do !ist

    deallocate(green_l, tmp)
    call pop_sub()
  end subroutine ext_eigenstate_lip_sch

  ! ---------------------------------------------------------
  ! compute the extended eigenstate(s)
  subroutine calculate_ext_eigenstate(h, gr, np, diag, offdiag, order, energy, lead_pot, st)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),   intent(inout)      :: gr
    integer,        intent(in)      :: np                      ! number of interface points
    CMPLX,          intent(in)      :: diag(np, np, NLEADS)    ! lead hamiltonian (diagonal part)
    CMPLX,          intent(in)      :: offdiag(np, np, NLEADS) ! hopping operator V^T
    integer,        intent(in)      :: order                   ! discretization order
    FLOAT,          intent(in)      :: energy                  ! total energy
    FLOAT,          intent(in)      :: lead_pot(NLEADS)        ! lead potential at t=0
    type(states_t), intent(inout)   :: st                      ! states

    integer  :: eigenstate_type

    call push_sub('td_trans_rti.calculate_extended_eigenstate')

    eigenstate_type = 2

    select case(eigenstate_type)
    case(1) ! reference implementation, only feasible in 1D (maybe 2D)
      !call ext_eigenst_imag_green(h, gr, np, diag, offdiag, order, energy, lead_pot, dx, st, phase, ik)
    case(2) ! method by Florian Lorenzen and Heiko Appel
      call ext_eigenstate_lip_sch(h, gr, np, diag, offdiag, order, energy, st)
    case(3) ! method by Gianluca Stefannuci
     !call ext_eigenstate_stefanucci(h, gr, np, diag, offdiag, order, energy, lead_pot, dx, st, phase, ik)
    end select

    call pop_sub()
  end subroutine calculate_ext_eigenstate

  ! ---------------------------------------------------------
  ! Crank-Nicholson timestep with source and memory.
  ! Only non-interacting electrons for the moment, so no
  ! predictor-corrector scheme.
  subroutine cn_src_mem_dt(intf, mem_type, mem, sp_mem, st_intf, ks, h, st, gr, energy, dt,  &
    t, max_iter, sm_u, td_pot, timestep, src_prev, diag, offdiag, st_psi0, mem_s, &
    additional_terms, mapping)
    type(intface_t), target,     intent(in)    :: intf
    type(states_t),              intent(inout) :: st
    integer,                     intent(in)    :: max_iter
    CMPLX,                       intent(inout) :: st_intf(intf%np, st%st_start:st%st_end, st%d%nik, NLEADS, 0:max_iter)
    type(v_ks_t),                intent(in)    :: ks
    type(hamiltonian_t), target, intent(inout) :: h
    type(grid_t), target,        intent(inout) :: gr
    FLOAT, target,               intent(in)    :: energy
    FLOAT, target,               intent(in)    :: dt
    FLOAT, target,               intent(in)    :: t
    integer, target,             intent(in)    :: mem_type
    CMPLX, target,               intent(in)    :: mem(intf%np, intf%np, 0:max_iter, NLEADS)
    CMPLX, target,               intent(in)    :: sp_mem(intf%np*gr%f_der%der_discr%order, 0:max_iter, NLEADS)
    CMPLX, target,               intent(in)    :: sm_u(0:max_iter, NLEADS)
    FLOAT, target,               intent(in)    :: td_pot(0:max_iter, NLEADS)
    integer,                     intent(in)    :: timestep
    CMPLX, target,               intent(inout) :: src_prev(intf%np, 1, st%st_start:st%st_end, st%d%nik, NLEADS)
    CMPLX, target,               intent(in)    :: diag(intf%np, intf%np, NLEADS)
    CMPLX, target,               intent(in)    :: offdiag(intf%np, intf%np, NLEADS)   ! hopping operator V^T
    CMPLX, target,               intent(in)    :: st_psi0(:, :, :, :, :, :)
    CMPLX, target,               intent(in)    :: mem_s(intf%np, intf%np, 2, NLEADS)
    integer,                     intent(in)    :: additional_terms
    integer, target,             intent(in)    :: mapping(:)   ! the mapping

    integer            :: il, it, m, cg_iter, j, order, ierr, inp, groundstate
    integer, target    :: ist, ik
    CMPLX              :: factor, alpha, fac, f0
    CMPLX, allocatable :: tmp(:, :), tmp_wf(:), tmp_mem(:, :)
    CMPLX, allocatable :: ext_wf(:,:,:,:) ! NP+2*np, ndim, nst, nik

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
    inp = intf%np

    f0 = M_z1/(M_z1+M_zI*M_HALF*dt*energy)
    fac= (M_z1-M_zI*M_HALF*dt*energy)*f0

    ! hack to do the calculation of the groundstate
    if (m.eq.0) then
      !%Variable TDTransGroundState
      !%Type integer
      !%Default -1
      !%Section Transport
      !%Description
      !% Decides whether we have the groundstate calculation or the td run.
      !%End
      call loct_parse_int(check_inp('TDTransGroundState'), -1, groundstate)

      select case(groundstate)
      case(-1) ! do nothing
      case(0) ! td run with reading extended eigenstate
        j = (NP+2*inp)*st%d%dim*st%lnst*st%d%nik
        ALLOCATE(ext_wf(NP+2*inp, st%d%dim, st%st_start:st%st_end, st%d%nik), j)
        call read_binary(j, ext_wf(1:NP+2*inp, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik), 3, ierr, 'ext_eigenstate.obf')
        st%zpsi(1:NP, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik) = &
                      ext_wf(inp+1:NP+inp, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik)
        deallocate(ext_wf)
      case(1) ! calculate the extended eigenstate
        ! DON'T FORGET TO MAKE THE SIMULATION BOX BIGGER IN THE INP FILE
        j = NP*st%d%dim*st%lnst*st%d%nik
        call calculate_ext_eigenstate(h, gr, inp, diag, offdiag, order, energy, td_pot(0,:), st)
        call write_binary(j, st%zpsi(1:NP, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik), 3, ierr, 'ext_eigenstate.obf')
      end select
    end if

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
              tmp_mem(:, :) = mem(:, :, m, il)
              if (m.gt.0) tmp_mem(:, :) = tmp_mem(:, :) + mem(:, :, m-1, il)
              call calc_source_wf(max_iter, m, inp, il, offdiag(:,:,il), tmp_mem, dt, &
                                    st_psi0(:,:,1,ist,ik,il), sm_u, f0, fac, &
                                    lambda(m, 0, il, max_iter, sm_u), src_prev(:, 1, ist, ik, il))
            else
              tmp_mem(:, 1) = sp_mem(:, m, il)
              if (m.gt.0) tmp_mem(:, 1) = tmp_mem(:, 1) + sp_mem(:, m-1, il)
              call calc_source_wf_sp(max_iter, m, inp, il, offdiag(:, :, il), tmp_mem(:, 1), dt, order, &
                       gr%sb%dim, st_psi0(:,:,1,ist,ik,il), mem_s(:,:,:,il), mapping, sm_u, f0, fac, &
                       lambda(m, 0, il, max_iter, sm_u), src_prev(:, 1, ist, ik, il))
            end if
            call apply_src(intf, il, src_prev(:, 1, ist, ik, il), st%zpsi(:, :, ist, ik))
          end if
          ! 3. Add memory term.
          if(iand(additional_terms, mem_term_flag).ne.0) then
            do it = 0, m-1
              factor = -dt**2/M_FOUR*lambda(m, it, il, max_iter, sm_u) / &
                (sm_u(m, il)*sm_u(it, il))
              tmp_wf(:) = st_intf(:, ist, ik, il, it+1) + st_intf(:, ist, ik, il, it)
              if (mem_type.eq.1) then
                tmp_mem(:, :) = mem(:, :, m-it, il)
                if ((m-it).gt.0) tmp_mem(:, :) = tmp_mem(:, :) + mem(:, :, m-it-1, il)
                call apply_mem(tmp_mem, il, intf, tmp_wf, st%zpsi(:, :, ist, ik), factor)
              else
                tmp_mem(:, 1) = sp_mem(:, m-it, il)
                if ((m-it).gt.0) tmp_mem(:, 1) = tmp_mem(:, 1) + sp_mem(:, m-it-1, il)
                call apply_sp_mem(tmp_mem(:,1), il, intf, tmp_wf, st%zpsi(:, :, ist, ik),&
                  factor, mem_s(:,:,:,il),order,gr%sb%dim, mapping)
              end if
            end do
          end if
        end do

        ! TODO: for many states parallel
        ! 4. Solve linear system (1 + i \delta H_{eff}) st%zpsi = tmp.
        cg_iter = cg_max_iter
        tmp(1:NP, 1) = st%zpsi(1:NP, 1, ist, ik)
        if (mem_type.eq.1) then
          call zconjugate_gradients(NP, st%zpsi(:, 1, ist, ik), tmp(:, 1), &
            h_eff_backward, h_eff_backwardt, zmf_dotp_aux, cg_iter, threshold=cg_tol)
        else
          call zconjugate_gradients(NP, st%zpsi(:, 1, ist, ik), tmp(1:NP, 1), &
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

!write(*,*) calc_current(gr, st)

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

    character(len=10) :: c_np
    write(c_np,'(i10)') np
    do j = 1, np
      write(*,'(a,i3,a,'//trim(c_np)//'e20.9)') "m(",j,",i)=",real(matr(j,1:np))
    end do
!    do j = 1, np
!      do i = 1, np
!        write(*,'(a,i3,a,i3,a,2e20.9)') "m(",j,",",i,")=",matr(j,i)
!      end do
!    end do
   end subroutine write_matrix


  ! ---------------------------------------------------------
  ! Progagate backwards with: 1 + i \delta H_{eff}
  ! Used by the iterative linear solver.
  subroutine h_eff_backward(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp(:, :)
    call push_sub('td_trans_rti.h_eff_backward')
    
    ALLOCATE(tmp(gr_p%m%np_part, 1),gr_p%m%np_part)
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
    
    ALLOCATE(tmp(gr_p%m%np_part, 1),gr_p%m%np_part)
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
    
    ALLOCATE(tmp(gr_p%m%np_part, 1),gr_p%m%np_part)
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
    
    ALLOCATE(tmp(gr_p%m%np_part, 1), gr_p%m%np_part)
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
  ! Apply source coefficient: zpsi <- zpsi + src_wf
  subroutine apply_src(intf, il, src_wf, zpsi)
    type(intface_t), intent(in)    :: intf
    integer,         intent(in)    :: il
    CMPLX,           intent(in)    :: src_wf(:)
    CMPLX,           intent(inout) :: zpsi(:, :)

    call push_sub('td_trans_rti.apply_src')

    ! Do not use use BLAS here.
    zpsi(intf%index_range(1, il):intf%index_range(2, il), 1) = &
      zpsi(intf%index_range(1, il):intf%index_range(2, il), 1) + src_wf(:)

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
