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
!! $Id: target_tdlocal_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_tdlocal(gr, tg, td)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg
    type(td_t),       intent(in)    :: td

    type(block_t)       :: blk
    PUSH_SUB(target_init_tdlocal)

    if(parse_block(datasets_check('OCTTdTarget'), blk)==0) then
      call parse_block_string(blk, 0, 0, tg%td_local_target)
      call conv_to_C_string(tg%td_local_target)
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))
      call parse_block_end(blk)
    else
      message(1) = 'If OCTTargetMode = oct_targetmode_td, you must suppy a OCTTDTarget block.'
      call messages_fatal(1)
    end if
    tg%dt = td%dt
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
    SAFE_DEALLOCATE_P(tg%rho)
    SAFE_DEALLOCATE_P(tg%td_fitness)
    POP_SUB(target_end_tdlocal)
  end subroutine target_end_tdlocal


  ! ----------------------------------------------------------------------
  subroutine target_output_tdlocal(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    integer :: ierr
    PUSH_SUB(target_output_tdlocal)
    
    call loct_mkdir(trim(dir))
    call target_build_tdlocal(tg, gr, M_ZERO)
    if(outp%how /= 0) then
      call dio_function_output(outp%how, trim(dir), 'td_local_target', gr%mesh, &
        tg%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
    end if

    POP_SUB(target_output_tdlocal)
  end subroutine target_output_tdlocal
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_tdlocal(tg) result(j1)
    type(target_t),   intent(inout) :: tg

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
  subroutine target_chi_tdlocal(gr, chi_out)
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: chi_out

    integer :: ik, idim, ist, ip
    PUSH_SUB(target_chi_tdlocal)

    !We assume that there is no time-independent operator.
    forall(ik = 1:chi_out%d%nik, ist = chi_out%st_start:chi_out%st_end, idim = 1:chi_out%d%dim, ip = 1:gr%mesh%np)
      chi_out%zpsi(ip, idim, ist, ik) = M_z0
    end forall


    POP_SUB(target_chi_tdlocal)
  end subroutine target_chi_tdlocal


  ! ---------------------------------------------------------
  !> 
  !!
  subroutine target_tdcalc_tdlocal(tg, gr, psi, time)
    type(target_t),      intent(inout) :: tg
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: psi
    integer,             intent(in)    :: time

    CMPLX, allocatable :: opsi(:, :)
    integer :: ist, ip
    PUSH_SUB(target_tdcalc_tdlocal)

    tg%td_fitness(time) = M_ZERO

    !!!! WARNING Here one should build the time-dependent target.
    select case(psi%d%ispin)
    case(UNPOLARIZED)
      ASSERT(psi%d%nik  ==  1)
      SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))
      opsi = M_z0
      do ist  = psi%st_start, psi%st_end
        do ip = 1, gr%mesh%np
          opsi(ip, 1) = tg%rho(ip) * psi%zpsi(ip, 1, ist, 1)
        end do
        tg%td_fitness(time) = &
          tg%td_fitness(time) + psi%occ(ist, 1) * &
          real(zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, 1), opsi(:, :)), REAL_PRECISION)
      end do
      SAFE_DEALLOCATE_A(opsi)
    case(SPIN_POLARIZED)
      message(1) = 'Error in target.target_tdcalc: spin_polarized.'
      call messages_fatal(1)
    case(SPINORS)
      message(1) = 'Error in target.target_tdcalc: spinors.'
      call messages_fatal(1)
    end select

    POP_SUB(target_tdcalc_tdlocal)
  end subroutine target_tdcalc_tdlocal
  ! ----------------------------------------------------------------------


  !----------------------------------------------------------
  subroutine target_build_tdlocal(tg, gr, time)
    type(target_t), intent(inout) :: tg
    type(grid_t),   intent(in)    :: gr
    FLOAT,          intent(in)    :: time

    integer :: ip
    FLOAT :: xx(MAX_DIM), rr, re, im

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
