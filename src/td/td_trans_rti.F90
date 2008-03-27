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
    CMPLX, pointer   :: src_factor(:)         ! iter
    CMPLX, pointer   :: diag(:, :, :)         ! h
    CMPLX, pointer   :: offdiag(:, :, :)      ! hopping operator V^T (np,np,nleads)
    CMPLX, pointer   :: mem_s(:, :, :, :)     ! the matrices to diagonalize coeff0

    type(intface_t)  :: intface

    CMPLX, pointer   :: st_intface(:, :, :, :, :) ! np_intf, nst, nik, nleads, max_iter
    CMPLX, pointer   :: st_sincos(:, :, :) ! np, (sin=1;cos=2,modsin=3,modcos=4), nleads
    CMPLX, pointer   :: st_psi0(:, :, :, :, :, :) ! np, (lead=1;center=2), ndim, nst, nik, nleads
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

    call lead_init(trans%lead, max_iter, dt)

    call td_trans_rti_write_info(trans, st, max_iter, gr%f_der%der_discr%order, trans%mem_type)

    call memory_init(trans%intface, dt/M_TWO, max_iter, gr%f_der%der_discr%lapl,&
                      trans%mem_coeff, trans%mem_sp_coeff, trans%mem_s, trans%diag, trans%offdiag, gr%sb%dim, &
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
energy = 0.2
    call source_init(st, trans%src_factor, trans%st_sincos, trans%st_psi0, trans%st_phase, dt, energy, q, &
                     max_iter, trans%intface%np, st%nst, st%d%nik, gr, gr%f_der%der_discr%order, NP)

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
    call source_end(trans%src_factor, trans%st_sincos, trans%st_phase, trans%st_psi0)
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

!  subroutine calculate_green_function()
!  end subroutine calculate_green_function

  ! ---------------------------------------------------------
  ! calculate phase shift at interface
  subroutine calculate_phase_shift(order, np, q, st, u, der_coeff, phase_shift)
    integer,         intent(in)  :: order
    integer,         intent(in)  :: np
    FLOAT,           intent(in)  :: q
    type(states_t),  intent(in)  :: st
    FLOAT,           intent(in)  :: u(:) ! the lead potential at time 0
    FLOAT,           intent(in)  :: der_coeff(:, :) ! weights of the derivative
    FLOAT,           intent(inout) :: phase_shift(:, :, :) ! nst, nik, NLEADS

    integer         :: lshift, length, ik, ist
    FLOAT           :: denom, energy

    call push_sub('td_trans_rti.calculate_phase_shift')

    ! the number of derivative coefficients
    length = 2*order+1
    ! with increasing order the interface boundary shifts towards the center
    lshift = order-1

    ! no tight binding yet
    energy = M_HALF*q**2

    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        ! left lead
        if (abs(st%zpsi(1+lshift, 1, ist, ik)).lt.M_EPSILON) then
          phase_shift(ist, ik, LEFT) = M_ZERO
        else
          ! FIXME: use a two-sided but asymmetric derivative
          ! for now, calculate the forward derivative
          denom = sum(log(abs(st%zpsi(1+lshift:length+lshift, 1, ist, ik))**2)*der_coeff(:, LEFT))
          if (abs(denom).lt.M_EPSILON) then
            phase_shift(ist, ik, LEFT) = M_HALF*M_PI
          else
            ! no tight binding yet
            phase_shift(ist, ik, LEFT) = atan(M_TWO*sqrt(M_TWO*(energy-u(LEFT))) / denom)
          end if
        end if

        ! right lead
        if (abs(st%zpsi(np-lshift, 1, ist, ik)).lt.M_EPSILON) then
          phase_shift(ist, ik, RIGHT) = M_ZERO
        else
          ! FIXME: use a two-sided but asymmetric derivative smaller error
          ! calculate the backward derivative
            denom = sum(log(abs(st%zpsi(np-length+1-lshift:np-lshift, 1, ist, ik))**2)*der_coeff(:, RIGHT))
          if (abs(denom).lt.M_EPSILON) then
            phase_shift(ist, ik, RIGHT) = sign(M_HALF*M_PI, denom)
          else
            ! no tight binding yet
            phase_shift(ist, ik, RIGHT) = atan(M_TWO*sqrt(M_TWO*(energy-u(RIGHT))) / denom)
          end if
        end if
      end do
    end do

    call pop_sub()
  end subroutine calculate_phase_shift

  ! ---------------------------------------------------------
  ! calculate scaling factor for matching the central wave function
  ! to the lead wave function of the form sqrt(2)*sin(q*x+delta)
  FLOAT function scaling_factor(phase, zpsi, np, lshift, der_coeff, length, q)
    FLOAT,      intent(in) :: phase(:)   ! NLEADS; phase shift at each interface
    CMPLX,      intent(in) :: zpsi(:)! np
    integer,    intent(in) :: np     ! number of gridpoints
    integer,    intent(in) :: lshift   ! shift (discr.-order dependent)
    FLOAT,      intent(in) :: der_coeff(:, :) ! derivative coefficients
    integer,    intent(in) :: length ! size of the coeffs array
    FLOAT,      intent(in) :: q   ! q = sqrt(2*(en-u))

    FLOAT   :: l_fac, r_fac
    call push_sub('td_trans_rti.scaling_factor')

    l_fac = M_ONE
    r_fac = M_ONE

    ! left lead
    if (abs(zpsi(lshift)).gt.M_EPSILON) then ! scale to match 1st derivatives
      l_fac = sqrt(M_TWO)*abs(sin(phase(LEFT))/zpsi(lshift))
    end if
! TODO: check if all is correct when zpsi==0
    ! right lead
    if (abs(zpsi(np-lshift)).gt.M_EPSILON) then
       r_fac = sqrt(M_TWO)*abs(sin(phase(RIGHT))/zpsi(np-lshift))
    end if

    scaling_factor = sqrt(l_fac**2 + r_fac**2)

    call pop_sub()
  end function scaling_factor
  
  ! ---------------------------------------------------------
  ! compute the semi-infinite surface green function (half) analyticly
  ! only the diagonalization process is numerical done
  subroutine green_lead(energy, diag, offdiag, np, il, green)
    FLOAT,               intent(in)  :: energy
    CMPLX,               intent(in)  :: diag(:, :)
    CMPLX,               intent(in)  :: offdiag(:, :)
    integer,             intent(in)  :: np ! number of interface points
    integer,             intent(in)  :: il
    CMPLX,               intent(out) :: green(:, :)

    integer             :: i
    CMPLX, allocatable  :: x(:,:), s(:,:), sub_o(:,:), sub_d(:,:), d(:), tmp(:,:)
    FLOAT               :: det

    call push_sub('td_trans_rti.green_lead')

    ALLOCATE(x(2*np,2*np),4*np**2)
    ALLOCATE(s(2*np,2*np),4*np**2)
    ALLOCATE(sub_o(np,np),np**2)
    ALLOCATE(sub_d(np,np),np**2)
    ALLOCATE(tmp(np,np),np**2)
    ALLOCATE(d(2*np),2*np)

    ! Algorithm taken from A. Umerski, Closed-form solutions to surface Green`s functions
    ! 1. create matrix x ( x = {{0,offdiag^(-1)},{-offdiag^T,(energy-diag)*offdiag^(-1)}} )
    ! 2. compute diagonalization matrix s, s^(-1)*x*s = d
    ! 3. extract submatrices ( s = {{o1,o2},{o3,o4}}; d = {{d1,0},{0,d2}} )
    ! 4. calculate g = o2*(offdiag*o2*d2)^(-1) for left lead or
    !              g = o1*d1*(o1*offdiag^(T))^(-1) for right lead

    ! 1. create matrix x ( x = {{0,offdiag^(-1)},{-offdiag^T,(energy-diag)*offdiag^(-1)}} )
    x(1:np,1:np) = M_z0
    ! use sub_o as tmp variable
    sub_o = offdiag
    if (il.eq.LEFT) then
      call lalg_invert_upper_triangular(np,  sub_o)
    else
      call lalg_invert_lower_triangular(np,  sub_o)
    end if
    x(1:np,np+1:2*np) = sub_o
    do i=1, np
      x(np+i,1:np) = -offdiag(1:np,i)
    end do
    ! use sub_d as tmp variable
    sub_d(:, :) = -diag(:, :)
    do i=1, np
      sub_d(i,i) = sub_d(i,i) + energy + CNST(1e-10)*M_zI
    end do
    if (il.eq.LEFT) then
      call lalg_trmm(np,np,'U','N','R',M_z1,sub_o,sub_d)
    else
      call lalg_trmm(np,np,'L','N','R',M_z1,sub_o,sub_d)
    end if
    x(np+1:2*np,np+1:2*np) = sub_d

    ! 2. compute diagonalization matrix s, s^(-1)*x*s = d
    s(:, :) = x(:, :)
    ! TODO: there is still a problem near the transversal eigenenergies,
    ! there the zgeev,zgeevx routine does not work
    ! a bigger imaginary eta does not solve this issue.
    call lalg_eigensolve_nonh(2*np, s, d)

    ! 3. extract submatrices ( S = {{o1,o2},{o3,o4}}; D = {{d1,0},{0,d2}} )
    sub_d = M_z0
    if (il.eq.LEFT) then ! left side --> o2, d2
      sub_o = s(1:np,np+1:2*np)
      do i=1, np
        sub_d(i,i) = d(np+i)
      end do
!write(*,*) 'd', d(np+1:2*np)
    else ! right side --> o1, d1
      sub_o = s(1:np,1:np)
      do i=1, np
        sub_d(i,i) = d(i)
      end do
!write(*,*) 'd', d(1:np)
    end if

    ! 4. calculate g
    if (il.eq.LEFT) then ! g = o2*(offdiag*o2*d2)^(-1)
      call zsymm('R','U',np,np,M_z1,sub_d,np,sub_o,np,M_z0,tmp,np)
      call lalg_trmm(np,np,'U','N','L',M_z1,offdiag,tmp)
      det = lalg_inverter(np, tmp, invert = .true.)
      call zgemm('N','N',np,np,np,M_z1,sub_o,np,tmp,np,M_z0,green,np)
    else ! g = o1*d1*(offdiag^T*o1)^(-1)
      call zsymm('R','U',np,np,M_z1,sub_d,np,sub_o,np,M_z0,tmp,np)
      call lalg_trmm(np,np,'L','T','L',M_z1,offdiag,sub_o)
      det = lalg_inverter(np, sub_o, invert = .true.)
      call zgemm('N','N',np,np,np,M_z1,tmp,np,sub_o,np,M_z0,green,np)
    end if
    ! now check if we have the correct solution: Im(Tr(g))<0
    ! if this is not the case just take the conjugate
    det = aimag(green(1,1))
    do i=2,np
      det = det + aimag(green(i,i))
    end do
    if (det.gt.M_ZERO) then
      green = conjg(green)
    end if
write(*,*) 'dos', -det

    deallocate(x, s, sub_o, sub_d, d, tmp)
    call pop_sub()
  end subroutine green_lead

  ! ---------------------------------------------------------
  ! check if the matrix is not just real valued
  logical function is_complex(matrix, np)
    CMPLX,          intent(in) :: matrix(np, np)
    integer,        intent(in) :: np

    integer :: i,j
    call push_sub('td_trans_rti.is_complex')

    is_complex = .false.
    mat_loop: do j=1,np
      do i=1,np
        if (abs(aimag(matrix(j,i))).gt.CNST(1e-10)) then
          is_complex = .true.
          exit mat_loop
        end if
      end do
    end do mat_loop
    call pop_sub()
  end function is_complex


  ! ---------------------------------------------------------
  ! compute the extended eigenstate(s) with imaginary part of the green function
  ! needs much memory: matrix do diagonalize has the size NP**2 
  subroutine ext_eigenst_imag_green(h, gr, np, diag, offdiag, order, energy, lead_pot, dx, st, phase, ik)
    type(hamiltonian_t), intent(in) :: h
    type(grid_t),   intent(in)      :: gr
    integer,        intent(in)      :: np                      ! number of interface points
    CMPLX,          intent(in)      :: diag(np, np, NLEADS)    ! lead hamiltonian (diagonal part)
    CMPLX,          intent(in)      :: offdiag(np, np, NLEADS) ! hopping operator V^T
    integer,        intent(in)      :: order                   ! discretization order
    FLOAT,          intent(in)      :: energy                  ! total energy
    FLOAT,          intent(in)      :: lead_pot(NLEADS)        ! lead potential at t=0
    FLOAT,          intent(in)      :: dx                      ! spacing
    type(states_t), intent(inout)   :: st                      ! states
    FLOAT,          intent(inout)   :: phase(:, :, :)          ! phaseshift at interface
    integer,        intent(in)      :: ik                      ! the index of the k-point

    FLOAT, allocatable :: der_coeff(:, :) ! the weights of the derivative
    CMPLX, allocatable :: green_l(:,:,:)
    CMPLX, allocatable :: g_cc(:,:)
    FLOAT, allocatable :: im_g_cc(:, :) ! imag of g_cc
    FLOAT, allocatable :: eigenvalues(:)
    integer            :: i
    FLOAT              :: eta, dos, q, fac

    call push_sub('td_trans_rti.ext_eigenst_imag_green')

    ALLOCATE(der_coeff(2*order+1, NLEADS),(2*order+1)*NLEADS)
    ALLOCATE(green_l(np, np, NLEADS),np**2*NLEADS)
    ALLOCATE(g_cc(NP, NP),NP**2)
    ALLOCATE(im_g_cc(NP, NP),NP**2)
    ALLOCATE(eigenvalues(NP),NP)

    ! 0. prepare some stuff for later
    call deriv_coeffs(order, dx, der_coeff(:,LEFT),  +1)
    call deriv_coeffs(order, dx, der_coeff(:,RIGHT), -1)
    !write(*,*) 'left stencil weights', der_coeff(:,LEFT)
    !write(*,*) 'right stencil weights', der_coeff(:,RIGHT)

    ! 1. calculate the extended eigenstate (unbounded)
    ! 1.1 calculate the retarded green function
    call green_lead(energy, diag(:,:,LEFT), offdiag(:,:,LEFT), np, LEFT, green_l(:,:,LEFT))
    ! TODO: apply_coupling
    call green_lead(energy, diag(:,:,RIGHT), offdiag(:,:,RIGHT), np, RIGHT, green_l(:,:,RIGHT))
    ! TODO: apply_coupling
    ! 1.2 calculate the green function for the central region
    g_cc = M_z0
    call nl_operator_op_to_matrix_cmplx(gr%f_der%der_discr%lapl, g_cc)
    g_cc = M_HALF*g_cc ! -H_cc
! FIXME: add all potentials
    do i=1, NP
      g_cc(i,i) = g_cc(i,i) -  h%vhxc(i, ik) - h%ep%vpsl(i)
    end do
    if (is_complex(green_l(:,:,LEFT),np).or.is_complex(green_l(:,:,RIGHT),np)) then
      eta = M_z0
    else
      eta = CNST(1E-10) ! check if not to small
      message(1) = 'Green function has no imaginary part!'
      message(2) = 'Energy probably to low or to high.'
      call write_warning(2)
    end if
    do i=1, NP
      g_cc(i,i) = energy + M_zI*eta
    end do
    g_cc(1:np,1:np) = g_cc(1:np,1:np) - green_l(1:np,1:np,LEFT)
    g_cc(NP-np+1:NP,NP-np+1:NP) = g_cc(NP-np+1:NP,NP-np+1:NP) - green_l(1:np,1:np,RIGHT)
    call lalg_sym_inverter('U', NP, g_cc)
    call make_symmetric(g_cc, NP)
    ! 1.3 take the imaginary part
    im_g_cc = aimag(g_cc)
    ! 1.4 calculate the density of states (without factor)
    dos = M_ZERO
    do i=1, NP
     dos = dos + im_g_cc(i,i)
    end do
    if (dos.eq.M_ZERO) then
      message(1) = 'Error in calculating the extended eigenstate:'
      message(2) = 'Density of states for Im(G_cc) is zero, can not continue!'
      call write_fatal(2)
    end if
    ! 1.5 diagonalize
    im_g_cc(:, :) = im_g_cc(:, :) / dos
    call lalg_eigensolve(NP, im_g_cc, im_g_cc, eigenvalues)
    ! 1.6 take the physical important eigenvalues and their corresponding eigenvectors = ext. eigenstates
    ! since the lead potentials must both be equal at t=0, the number of states is always 2 (or 0)
    if (abs(eigenvalues(NP)+eigenvalues(NP-1)-M_ONE).gt.CNST(1e-10)) then
      message(1) = 'Error in calculating the extended eigenstate:'
      message(2) = 'The two largest eigenvalues do not add up to 1!'
      call write_fatal(2)
    end if
    ! yet only tight binding
    st%zpsi(:,1,st%st_start,ik) = im_g_cc(:,NP)
    st%zpsi(:,1,st%st_start+1,ik) = im_g_cc(:,NP-1)

write(*,*) 'state 1', st%zpsi(:,1,st%st_start,ik)
write(*,*) 'state 2', st%zpsi(:,1,st%st_start+1,ik)
    ! 2. calculate phase shift of the plane waves in the leads (eq. (41) in paper)
    q = sqrt(M_TWO*(energy-lead_pot(LEFT)))
!    q = M_PI/(gr%sb%lsize(1)-(order-1)*dx)
    call calculate_phase_shift(order, NP, q, st, lead_pot, der_coeff, phase)

    ! 3. calculate the scaling to match to the lead eigenstates
    fac = scaling_factor(phase(st%st_start, 1, :), st%zpsi(:,1,st%st_start,ik), np, order-1, der_coeff, 2*order+1, q)
    st%zpsi(:,1,st%st_start,ik) = fac*st%zpsi(:,1,st%st_start,ik)
    fac = scaling_factor(phase(st%st_start+1, 1, :), st%zpsi(:,1,st%st_start,ik), np, order-1, der_coeff, 2*order+1, q)
    st%zpsi(:,1,st%st_start,ik) = fac*st%zpsi(:,1,st%st_start,ik)
    !write(6,*) fac

  ! 3.1 calculate the derivatives at the boundaries
!  write(*,*) 'vorher'
!      write(*,*) phase(LEFT), sin(phase(LEFT))**2, abs(st%zpsi(1, 1, ist, ik))**2
!      write(*,*) phase(RIGHT), sin(phase(RIGHT))**2, log(abs(st%zpsi(np, 1, ist, ik))**2)
!      c = M_ONE/sqrt(dx+abs(st%zpsi(1, 1, ist, ik))**2/(M_TWO*q)+abs(st%zpsi(np, 1, ist, ik))**2/(M_TWO*q))
!write(*,*) c
!      st%zpsi(:, :, ist, ik) = c*st%zpsi(:, :, ist, ik)

!  st_phase(:) = M_ZERO ! FIXME

    deallocate(der_coeff, green_l, g_cc, im_g_cc)
    call pop_sub()
  end subroutine ext_eigenst_imag_green

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
!    y(:) = x(:) ! no preconditioner

    deallocate(diag)

    call pop_sub()
  end subroutine preconditioner

  ! ---------------------------------------------------------
  ! compute the extended eigenstate(s) 
  subroutine ext_eigenstate_lip_sch(h, gr, np, diag, offdiag, order, energy, lead_pot, dx, st, phase, ik)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),   intent(inout)      :: gr
    integer,        intent(in)      :: np                      ! number of interface points
    CMPLX,          intent(in)      :: diag(np, np, NLEADS)    ! lead hamiltonian (diagonal part)
    CMPLX,          intent(in)      :: offdiag(np, np, NLEADS) ! hopping operator V^T
    integer,        intent(in)      :: order                   ! discretization order
    FLOAT, target,  intent(in)      :: energy                  ! total energy
    FLOAT,          intent(in)      :: lead_pot(NLEADS)        ! lead potential at t=0
    FLOAT,          intent(in)      :: dx                      ! spacing
    type(states_t), intent(inout)   :: st                      ! states
    FLOAT,          intent(inout)   :: phase(:, :, :)          ! phaseshift at interface
    integer,        intent(in)      :: ik                      ! the index of the k-point

    CMPLX, target, allocatable :: green_l(:,:,:)
    CMPLX, allocatable :: tmp(:, :)
    integer            :: i, id, iter, n(3)
    integer, pointer   :: lxyz(:,:)
    FLOAT              :: en, lsize(3), q(3)
    FLOAT              :: dres

    call push_sub('td_trans_rti.ext_eigenstate_lip_sch')

    ALLOCATE(green_l(np, np, NLEADS),np**2*NLEADS)
    ALLOCATE(tmp(NP_PART, 1), NP_PART)

    ! now generate the extended eigenstate with separated x-wavefunction
    ! yet only for a box psi=psi(x)*psi(y)*psi(z)
    ! e_tot = e_x + e_y + e_z, we need e_x for q (for box-sized leads)
    ! FIXME: (probably) run trough all possible linear combinations
    ! yet only the highest state of the separated function is multiplied

    ! 1. compute the unpertubed eigenstates
    ! for now set manually psi(x)*phi(y,z)
    energy_p => energy
    lxyz => gr%m%lxyz
    lsize = gr%sb%lsize
    en = energy
    ! TODO: run over all states
    ! subtract the transversal energy from the total energy to get the longitudinal part
    ! which is the energy used for the transport
    do id=2, gr%sb%dim ! FIXME: when the last gridpoint is not the border, recalculate transversal energy
      q(id) = M_PI/(M_TWO*(lsize(id)+gr%sb%h(id)))
      n(id) = sqrt(M_TWO*en)/q(id) ! use for now the largest transversal energy
      en = en - M_HALF*(n(id)*q(id))**2
    end do
!write(*,*) 'ny, ex, ey',n(2), en, M_HALF*(n(2)*q(2))**2
!write(*,*) 'ny, ex',n(2), en
    ! FIXME: tight binding approximation
    q(1) = sqrt(M_TWO*en)
    ! the rest is the energy for the transport direction
    if(en < M_ZERO) then
      write(message(1), '(a,f14.6,a)') "The input energy : '", energy, &
                  "' is smaller than the groundstate energy for the given system."
      call write_fatal(1)
    end if
    n(1) = 1
    do i=1,NP
      st%zpsi(i, 1, 1, ik) = exp(M_zI*n(1)*q(1)*lxyz(i,1)*gr%sb%h(1))
      do id=2, gr%sb%dim
        if (mod(n(id),2).eq.0) then
          st%zpsi(i, 1, 1, ik) = st%zpsi(i, 1, 1, ik)*sin(n(id)*q(id)*lxyz(i,id)*gr%sb%h(id))
        else
          st%zpsi(i, 1, 1, ik) = st%zpsi(i, 1, 1, ik)*cos(n(id)*q(id)*lxyz(i,id)*gr%sb%h(id))
        end if
      end do
    end do
    ! calculate right hand side (e-T-V(lead)-sum(a)[H_ca*g_a*H_ac]
    tmp(:,:) = M_z0
    call zkinetic(h, gr, st%zpsi(:, :, 1, ik), tmp(:,:))
    tmp(1:NP, :) = energy*st%zpsi(1:NP, :, 1, ik) - tmp(1:NP, :)
    ! TODO: the static potential of the lead
    call green_lead(energy, diag(:, :, LEFT), offdiag(:, :, LEFT), np, LEFT, green_l(:, :, LEFT))
    call apply_coupling(green_l(:, :, LEFT), offdiag(:, :, LEFT), green_l(:, :, LEFT), np, LEFT)
    call zsymv('U', np, -M_z1, green_l(:, :, LEFT), np, st%zpsi(1:np, 1, 1, 1), 1, M_z1, tmp(1:np, 1), 1)
    call green_lead(energy, diag(:, :, RIGHT), offdiag(:, :, RIGHT), np, RIGHT, green_l(:, :, RIGHT))
    call apply_coupling(green_l(:, :, RIGHT), offdiag(:, :, RIGHT), green_l(:, :, RIGHT), np, RIGHT)

    call zsymv('U', np, -M_z1, green_l(:, :, RIGHT), np, &
      st%zpsi(NP - np + 1:NP, 1, 1, 1), 1, M_z1, tmp(NP - np + 1:NP, 1), 1)

    ! now solve the equation
    green_l_p => green_l
    ! now solve the linear system to get the extended eigenstate
    iter = 5000
    ! zconjugate_gradients fails in some cases (wrong solution) and takes longer
    ! therefore take the symmetric quasi-minimal residual solver (QMR)
    call zqmr_sym(NP, st%zpsi(:, 1, 1, ik), tmp(:, 1), h_eff_lip_sch, preconditioner, &
                    iter, residue=dres, threshold=cg_tol, showprogress = .true.)
    write(*,*) 'iter =',iter, 'residue =', dres

    deallocate(green_l, tmp)
    call pop_sub()
  end subroutine ext_eigenstate_lip_sch

  ! ---------------------------------------------------------
  ! compute the extended eigenstate(s)
  subroutine calculate_ext_eigenstate(h, gr, np, diag, offdiag, order, energy, lead_pot, dx, st, phase, ik)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),   intent(inout)      :: gr
    integer,        intent(in)      :: np                      ! number of interface points
    CMPLX,          intent(in)      :: diag(np, np, NLEADS)    ! lead hamiltonian (diagonal part)
    CMPLX,          intent(in)      :: offdiag(np, np, NLEADS) ! hopping operator V^T
    integer,        intent(in)      :: order                   ! discretization order
    FLOAT,          intent(in)      :: energy                  ! total energy
    FLOAT,          intent(in)      :: lead_pot(NLEADS)        ! lead potential at t=0
    FLOAT,          intent(in)      :: dx                      ! spacing
    type(states_t), intent(inout)   :: st                      ! states
    FLOAT,          intent(out)     :: phase(:, :, :)          ! phaseshift at interface
    integer,        intent(in)      :: ik                      ! the index of the k-point

    integer  :: eigenstate_type

    call push_sub('td_trans_rti.calculate_extended_eigenstate')

    eigenstate_type = 2

    select case(eigenstate_type)
    case(1) ! reference implementation, only feasible in 1D (maybe 2D)
      call ext_eigenst_imag_green(h, gr, np, diag, offdiag, order, energy, lead_pot, dx, st, phase, ik)
    case(2) ! method by Florian Lorenzen and Heiko Appel
      call ext_eigenstate_lip_sch(h, gr, np, diag, offdiag, order, energy, lead_pot, dx, st, phase, ik)
    case(3) ! method by Gianluca Stefannuci
     !call ext_eigenstate_stefanucci(h, gr, np, diag, offdiag, order, energy, lead_pot, dx, st, phase, ik)
    end select

    call pop_sub()
  end subroutine calculate_ext_eigenstate

  ! ---------------------------------------------------------
  ! Crank-Nicholson timestep with source and memory.
  ! Only non-interacting electrons for the moment, so no
  ! predictor-corrector scheme.
  subroutine cn_src_mem_dt(intf, mem_type, mem, sp_mem, st_intf, ks, h, st, gr, dt,  &
    t, max_iter, sm_u, u, timestep, src_factor, diag, offdiag, st_sincos, st_psi0, st_phase, mem_s, &
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
    CMPLX, target,               intent(in)    :: diag(intf%np, intf%np, NLEADS)
    CMPLX, target,               intent(in)    :: offdiag(intf%np, intf%np, NLEADS)   ! hopping operator V^T
    CMPLX, target,               intent(in)    :: st_sincos(:, :, :) ! sin & cos vectors
    CMPLX, target,               intent(in)    :: st_psi0(:, :, :, :, :, :)
    FLOAT, target,               intent(inout) :: st_phase(:, :, :)   ! phaseshift at interface
    CMPLX, target,               intent(in)    :: mem_s(intf%np, intf%np, 2, NLEADS)
    integer,                     intent(in)    :: additional_terms
    integer, target,             intent(in)    :: mapping(:)   ! the mapping

    integer            :: il, it, m, cg_iter, j, order
    integer, target    :: ist, ik
    FLOAT              :: energy
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

!if (m.eq.0) then
! FIXME: this does NOT belong here
  energy = (timestep/M_TEN+3.825)*(M_PI/(gr%sb%lsize(2)+gr%sb%h(2)))**2/M_EIGHT 
!  energy = (timestep/M_TEN+3.4)*(M_PI/(gr%sb%lsize(2)+gr%sb%h(2)))**2/M_EIGHT 
!  energy = (timestep/M_TEN+1.175)*(M_PI/(gr%sb%lsize(2)+gr%sb%h(2)))**2/M_EIGHT
!  energy = (timestep/M_TEN+0.925)*(M_PI/(gr%sb%lsize(2)+gr%sb%h(2)))**2/M_EIGHT
!  energy = (timestep/M_TEN/M_TWO)
!  energy = (timestep/M_TEN**5+3.96750)*(M_PI/(gr%sb%lsize(2)+gr%sb%h(2)))**2/M_EIGHT 
!  energy = ((timestep+1)**2-0.1)*(M_PI/(gr%sb%lsize(2)+gr%sb%h(2)))**2/M_EIGHT 
!write(*,*) 'E/Ey', (timestep/M_TEN**5+3.96750)
!  energy = 0.1
  do ik = 1, st%d%nik
   ! call calculate_ext_eigenstate(h, gr, intf%np, diag, offdiag, order, energy, u(0,:), gr%sb%h(1), st, st_phase, ik)
  end do
!end if
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
       !   call apply_h_eff(h, gr, mem(:, :, 0, :), intf, -M_ONE, dt, t, ist, ik, st%zpsi(:, :, ist, ik))
        else
        !  call apply_h_eff_sp(h, gr, sp_mem(:, 0, :), intf, -M_ONE, dt, t, ist, ik, &
         !                     st%zpsi(:, :, ist, ik), mem_s(:, :, :, :), mapping)
        end if

        do il = 1, NLEADS
          ! 2. Add source term
          if(iand(additional_terms, src_term_flag).ne.0) then
            if (mem_type.eq.1) then
              call calc_source_wf(max_iter, m, il, st_phase(ist, ik, il), offdiag(:, :, il), src_factor(:), &
                mem(:, :, :, il), dt, intf%np, st_sincos(:, :, il), tmp_wf(:))
              !call calc_source_wf_new(max_iter, m, il, offdiag(:, :, il), src_factor(:), &
              !  mem(:, :, :, il), dt, intf%np, st_psi0(:,:,1,ist,ik,il), tmp_wf(:))
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
    write(c_np,'(i10)') np*2
    do j = 1, np
      write(*,'(a,i3,a,'//trim(c_np)//'e20.9)') "m(",j,",i)=",matr(j,1:np)
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
