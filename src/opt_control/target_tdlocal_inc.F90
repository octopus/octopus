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


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_tdlocal(gr, namespace, tg, td)
    type(grid_t),      intent(in)    :: gr
    type(namespace_t), intent(in)    :: namespace
    type(target_t),    intent(inout) :: tg
    type(td_t),        intent(in)    :: td

    type(block_t)       :: blk
    PUSH_SUB(target_init_tdlocal)

    tg%move_ions = ion_dynamics_ions_move(td%ions_dyn)
    tg%dt = td%dt

    !%Variable OCTTdTarget
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% (Experimental) If <tt>OCTTargetOperator = oct_tg_td_local</tt>, then you must supply
    !% a OCTTdTarget block. The block should only contain one element, a string cotaining the
    !% definition of the time-dependent local target, <i>i.e.</i> a function of x,y,z and t that 
    !% is to be maximized along the evolution.
    !%End
    if(parse_block(namespace, 'OCTTdTarget', blk)==0) then
      call parse_block_string(blk, 0, 0, tg%td_local_target)
      call conv_to_C_string(tg%td_local_target)
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))
      call parse_block_end(blk)
    else
      message(1) = 'If OCTTargetOperator = oct_tg_td_local, you must supply a OCTTdTarget block.'
      call messages_fatal(1)
    end if
    SAFE_ALLOCATE(tg%td_fitness(0:td%max_iter))
    tg%td_fitness = M_ZERO
    call target_build_tdlocal(tg, gr, M_ZERO)

    POP_SUB(target_init_tdlocal)
  end subroutine target_init_tdlocal


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_tdlocal(tg)
    type(target_t),   intent(inout) :: tg

    PUSH_SUB(target_end_tdlocal)

    SAFE_DEALLOCATE_A(tg%rho)
    SAFE_DEALLOCATE_A(tg%td_fitness)

    POP_SUB(target_end_tdlocal)
  end subroutine target_end_tdlocal


  ! ----------------------------------------------------------------------
  subroutine target_output_tdlocal(tg, namespace, space, gr, dir, ions, outp)
    type(target_t),    intent(inout) :: tg
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space
    type(grid_t),      intent(in)    :: gr
    character(len=*),  intent(in)    :: dir
    type(ions_t),      intent(in)    :: ions
    type(output_t),    intent(in)    :: outp

    integer :: ierr
    PUSH_SUB(target_output_tdlocal)
    
    call io_mkdir(trim(dir), namespace)
    call target_build_tdlocal(tg, gr, M_ZERO)
    call dio_function_output(outp%how(0), trim(dir), 'td_local_target', namespace, space, gr%mesh, &
      tg%rho, units_out%length**(-space%dim), ierr, ions = ions)

    POP_SUB(target_output_tdlocal)
  end subroutine target_output_tdlocal
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_tdlocal(tg) result(j1)
    type(target_t), intent(in) :: tg

    integer :: maxiter
    PUSH_SUB(target_j1_tdlocal)

    maxiter = size(tg%td_fitness) - 1
    j1 = M_HALF * tg%dt * tg%td_fitness(0) + & 
         M_HALF * tg%dt * tg%td_fitness(maxiter) + & 
         tg%dt * sum(tg%td_fitness(1:maxiter-1))


    POP_SUB(target_j1_tdlocal)
  end function target_j1_tdlocal


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_tdlocal(chi_out)
    type(states_elec_t), intent(inout) :: chi_out

    integer :: ik, ib
    PUSH_SUB(target_chi_tdlocal)

    !We assume that there is no time-independent operator.
    
    do ik = chi_out%d%kpt%start, chi_out%d%kpt%end
      do ib = chi_out%group%block_start, chi_out%group%block_end
        call batch_set_zero(chi_out%group%psib(ib, ik))
      end do
    end do

    POP_SUB(target_chi_tdlocal)
  end subroutine target_chi_tdlocal


  ! ---------------------------------------------------------
  !> 
  !!
  subroutine target_tdcalc_tdlocal(tg, gr, psi, time)
    type(target_t),      intent(inout) :: tg
    type(grid_t),        intent(in)    :: gr
    type(states_elec_t), intent(in)    :: psi
    integer,             intent(in)    :: time

    CMPLX, allocatable :: opsi(:, :), zpsi(:, :)
    integer :: ist, ip
    PUSH_SUB(target_tdcalc_tdlocal)

    tg%td_fitness(time) = M_ZERO

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:psi%d%dim))
    
    !!!! WARNING Here one should build the time-dependent target.
    select case(psi%d%ispin)
    case(UNPOLARIZED)
      ASSERT(psi%d%nik  ==  1)
      SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))
      opsi = M_z0
      do ist  = psi%st_start, psi%st_end

        call states_elec_get_state(psi, gr%mesh, ist, 1, zpsi)
        
        do ip = 1, gr%mesh%np
          opsi(ip, 1) = tg%rho(ip)*zpsi(ip, 1)
        end do
        
        tg%td_fitness(time) = tg%td_fitness(time) + psi%occ(ist, 1)*TOFLOAT(zmf_dotp(gr%mesh, psi%d%dim, zpsi, opsi))

      end do
      SAFE_DEALLOCATE_A(opsi)
    case(SPIN_POLARIZED)
      message(1) = 'Error in target.target_tdcalc: spin_polarized.'
      call messages_fatal(1)
    case(SPINORS)
      message(1) = 'Error in target.target_tdcalc: spinors.'
      call messages_fatal(1)
    end select

    SAFE_DEALLOCATE_A(zpsi)

    POP_SUB(target_tdcalc_tdlocal)
  end subroutine target_tdcalc_tdlocal
  ! ----------------------------------------------------------------------


  !----------------------------------------------------------
  subroutine target_build_tdlocal(tg, gr, time)
    type(target_t), intent(inout) :: tg
    type(grid_t),   intent(in)    :: gr
    FLOAT,          intent(in)    :: time

    integer :: ip
    FLOAT :: xx(gr%sb%dim), rr, re, im

    PUSH_SUB(target_build_tdlocal)

    do ip = 1, gr%mesh%np
      call mesh_r(gr%mesh, ip, rr, coords = xx)
      call parse_expression(re, im, gr%sb%dim, xx, rr, time, tg%td_local_target)
      tg%rho(ip) = re
    end do

    POP_SUB(target_build_tdlocal)
  end subroutine target_build_tdlocal
  !----------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
