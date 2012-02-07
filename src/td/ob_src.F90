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

!> Calculation of the source term for the modified Crank-Nicholson
!! propagator.

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
  use parser_m
  use math_m
  use messages_m
  use nl_operator_m
  use ob_interface_m
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
  subroutine ob_src_init(ob, st, intf)
    type(ob_terms_t),  intent(inout) :: ob
    type(states_t),    intent(in)    :: st
    type(interface_t), intent(in)    :: intf(:)

    integer :: il, np
    PUSH_SUB(ob_src_init)

    ! FIXME: spinor index is ignored here.
    do il=1, NLEADS
      np = intf(il)%np_intf
      SAFE_ALLOCATE(ob%lead(il)%src_prev(1:np, 1, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
      ob%lead(il)%src_prev = M_z0
    end do

    POP_SUB(ob_src_init)
  end subroutine ob_src_init


  ! ---------------------------------------------------------
  ! Calculates the source-wavefunction for the source-term (recursive version)
  ! S(m) = f*u(m)*u(m-1)*S(m-1) + dt**2/2*lambda(m,0)/u(m)*f0*(Q(m)+Q(m-1))*psi_c(0)
  subroutine calc_source_wf(maxiter, m, np, il, offdiag, mem, dt, psi0, u, f0, factor, lambda, src)
    integer, intent(in)    :: maxiter
    integer, intent(in)    :: m             ! m-th timestep.
    integer, intent(in)    :: np            ! intf%np, the size of the wavefunctions.
    integer, intent(in)    :: il            ! Which lead.
    CMPLX,   intent(in)    :: offdiag(:, :) ! Matrix V^T.
    CMPLX,   intent(in)    :: mem(np, np)   ! the effective memory coefficient
    FLOAT,   intent(in)    :: dt            ! Timestep.
    CMPLX,   intent(in)    :: psi0(:)       ! (np).
    CMPLX,   intent(in)    :: u(0:)
    CMPLX,   intent(in)    :: f0
    CMPLX,   intent(in)    :: factor
    CMPLX,   intent(in)    :: lambda
    CMPLX,   intent(inout) :: src(np)        ! Old wavefunction in, the new one out.
    
    CMPLX   :: tmp, alpha

    PUSH_SUB(calc_source_wf)
    if(m.eq.0) then
      ! initial src is V*psi0(outer) (precalculated in gs mode)
      tmp = -M_zI*dt*f0*u(0)
      !call lalg_gemv(np, np, tmp*M_zI*M_HALF*dt, mem, psi0, tmp, src)
      call lalg_symv(np, tmp*M_zI*M_HALF*dt, mem, psi0, tmp, src)
    else
      tmp   = factor*u(m)*u(m-1)
      if(m.gt.maxiter) then
        src(1:np) = tmp*src(1:np)
      else
        alpha = M_HALF*dt**2*f0*lambda/u(m)
!        call lalg_gemv(np, np, alpha, mem, psi0, tmp, src)
        call lalg_symv(np, alpha, mem, psi0, tmp, src)
      end if
    end if

    POP_SUB(calc_source_wf)
  end subroutine calc_source_wf


  ! ---------------------------------------------------------
  ! Calculates the source-wavefunction for the source-term (for sparse
  ! mem-coefficients).
  subroutine calc_source_wf_sp(maxiter, m, np, il, offdiag, sp_mem, dt, order, dim, psi0, &
    mem_s, mapping, u, f0, factor, lambda, src)
    integer, intent(in)    :: maxiter        ! Maximum timestep.
    integer, intent(in)    :: m              ! m-th timestep
    integer, intent(in)    :: np             ! intf%np, the size of the wavefunctions.
    integer, intent(in)    :: order
    integer, intent(in)    :: dim
    integer, intent(in)    :: il             ! Which lead.
    CMPLX,   intent(in)    :: offdiag(:, :)  ! matrix V^T.
    CMPLX,   intent(in)    :: sp_mem(:)
    FLOAT,   intent(in)    :: dt             ! Timestep.
    CMPLX,   intent(in)    :: psi0(:)        ! (np)
    CMPLX,   intent(in)    :: mem_s(:, :, :)
    integer, intent(in)    :: mapping(:)
    CMPLX,   intent(in)    :: u(0:maxiter)
    CMPLX,   intent(in)    :: f0
    CMPLX,   intent(in)    :: factor
    CMPLX,   intent(in)    :: lambda
    CMPLX,   intent(inout) :: src(:)         ! Old wavefunction in, the new one out.
    
    CMPLX,   allocatable :: tmem(:, :)

    PUSH_SUB(calc_source_wf_sp)

    SAFE_ALLOCATE(tmem(1:np, 1:np))
    call make_full_matrix(np, order, dim, sp_mem, mem_s, tmem, mapping)
    call calc_source_wf(maxiter, m, np, il, offdiag, tmem, dt, psi0, u, f0, factor, lambda, src)
    SAFE_DEALLOCATE_A(tmem)

    POP_SUB(calc_source_wf_sp)
  end subroutine calc_source_wf_sp


  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine ob_src_end(ob)
    type(ob_terms_t), intent(inout) :: ob

    integer :: il

    PUSH_SUB(ob_src_end)

    do il=1, NLEADS
      SAFE_DEALLOCATE_P(ob%lead(il)%src_prev)
    end do

    POP_SUB(ob_src_end)
  end subroutine ob_src_end
end module ob_src_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
