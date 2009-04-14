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
!! $Id: td_trans_src.F90 3030 2008-01-03 10:30:05Z nitsche $

! Calculation of the source term for the modified Crank-Nicholson
! propagator.

#include "global.h"

module ob_src_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use io_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use loct_parser_m
  use math_m
  use messages_m
  use nl_operator_m
  use ob_mem_m
  use ob_terms_m
  use profiling_m
  use simul_box_m
  use states_m
  use system_m
  use varinfo_m

  implicit none

  private
  public ::             &
    ob_src_init,        &
    ob_src_end,         &
    calc_source_wf,     &
    calc_source_wf_sp

contains
  
  ! ---------------------------------------------------------
  ! Allocate memory to calculate source term.
  subroutine ob_src_init(ob, st, np)
    type(ob_terms_t), intent(inout) :: ob
    type(states_t),   intent(in)    :: st
    integer,          intent(in)    :: np

    integer :: size
    call push_sub('ob_src.ob_src_init')

    ! FIXME: spinor index is ignored here.
    size = np*st%lnst*st%d%kpt%nlocal*NLEADS
    ALLOCATE(ob%src_prev(np, 1, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end, NLEADS), size)

    ob%src_prev = M_z0

    call pop_sub()
  end subroutine ob_src_init


  ! ---------------------------------------------------------
  ! Calculates the source-wavefunction for the source-term (recursive version)
  ! S(m) = f*u(m)*u(m-1)*S(m-1) + dt**2/2*lambda(m,0)/u(m)*f0*(Q(m)+Q(m-1))*psi_c(0)
  subroutine calc_source_wf(maxiter, m, np, il, offdiag, mem, dt, psi0, u, f0, factor, lambda, src)
    integer, intent(in)    :: maxiter
    integer, intent(in)    :: m             ! m-th timestep.
    integer, intent(in)    :: np            ! intf%np, the size of the wave functions.
    integer, intent(in)    :: il            ! Which lead.
    CMPLX,   intent(in)    :: offdiag(:, :) ! Matrix V^T.
    CMPLX,   intent(in)    :: mem(np, np)   ! the effective memory coefficient
    FLOAT,   intent(in)    :: dt            ! Timestep.
    CMPLX,   intent(in)    :: psi0(:, :)    ! (np, INNER/OUTER).
    CMPLX,   intent(in)    :: u(0:maxiter)
    CMPLX,   intent(in)    :: f0
    CMPLX,   intent(in)    :: factor
    CMPLX,   intent(in)    :: lambda
    CMPLX,   intent(inout) :: src(np)        ! Old wave function in, the new one out.
    
    CMPLX   :: tmp, alpha

    call push_sub('ob_src.calc_source_wf')

    if(m.eq.0) then
      src(1:np) = psi0(1:np, OUTER)
      if (il.eq.LEFT) then
        call ztrmv('U', 'N', 'N', np, offdiag, np, src, 1)
      else
        call ztrmv('L', 'N', 'N', np, offdiag, np, src, 1)
      end if
      tmp = -M_zI*dt*u(0)*f0
      call zsymv('U', np, tmp*M_zI*M_HALF*dt, mem, np, psi0(1:np, INNER), 1, tmp,  src, 1)
    else
      tmp   = factor*u(m)*u(m-1)
      alpha = M_HALF*dt**2*f0*lambda/u(m)
      call zsymv('U', np, alpha, mem, np, psi0(1:np, INNER), 1, tmp, src, 1)
    end if

    call pop_sub()
  end subroutine calc_source_wf


  ! ---------------------------------------------------------
  ! Calculates the source-wavefunction for the source-term (for sparse
  ! mem-coefficients).
  subroutine calc_source_wf_sp(maxiter, m, np, il, offdiag, sp_mem, dt, order, dim, psi0, &
    mem_s, mapping, u, f0, factor, lambda, src)
    integer, intent(in)    :: maxiter        ! Maximum timestep.
    integer, intent(in)    :: m              ! m-th timestep
    integer, intent(in)    :: np             ! intf%np, the size of the wave functions.
    integer, intent(in)    :: order
    integer, intent(in)    :: dim
    integer, intent(in)    :: il             ! Which lead.
    CMPLX,   intent(in)    :: offdiag(:, :)  ! matrix V^T.
    CMPLX,   intent(in)    :: sp_mem(:)
    FLOAT,   intent(in)    :: dt             ! Timestep.
    CMPLX,   intent(in)    :: psi0(:, :)     ! (np, INNER/OUTER)
    CMPLX,   intent(in)    :: mem_s(:, :, :)
    integer, intent(in)    :: mapping(:)
    CMPLX,   intent(in)    :: u(0:maxiter)
    CMPLX,   intent(in)    :: f0
    CMPLX,   intent(in)    :: factor
    CMPLX,   intent(in)    :: lambda
    CMPLX,   intent(inout) :: src(:)         ! Old wave function in, the new one out.
    
    CMPLX,   allocatable :: tmem(:, :)

    call push_sub('ob_src.calc_source_wf_sp')

    ALLOCATE(tmem(np, np), np**2)
    call make_full_matrix(np, order, dim, sp_mem, mem_s, tmem, mapping)
    call calc_source_wf(maxiter, m, np, il, offdiag, tmem, dt, psi0, u, f0, factor, lambda, src)
    deallocate(tmem)

    call pop_sub()
  end subroutine calc_source_wf_sp


  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine ob_src_end(ob)
    type(ob_terms_t), intent(inout) :: ob

    call push_sub('ob_src.ob_src_end')

    DEALLOCATE(ob%src_prev)

    call pop_sub()
  end subroutine ob_src_end
end module ob_src_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
