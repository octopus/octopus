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
  use global_m

  private
  public :: &
    type(ob_terms_t)

  integer, parameter :: &
    mem_term_flag  = 1, &
    src_term_flag  = 2, &
    save_cpu_time  = 1, &
    save_ram_usage = 2

  ! Additional source and memory terms.
  type ob_terms_t
    integer          :: mem_type                  ! 1: fast/lots of memory, 2: slow/little memory.
    integer          :: sp_length                 ! Length of the sparse array (as 1d array).
    integer          :: additional_terms          ! Shall we add source and memory term?
    CMPLX, pointer   :: mem_coeff(:, :, :, :)     ! Memory coefficients, for mem_type=1. (i, j, t, il) 
    CMPLX, pointer   :: mem_sp_coeff(:, :, :)     ! Memory coefficients, for mem_type=2, (sp_length, t, il)
    integer, pointer :: sp2full_map(:)            ! The mapping indices from sparse to full matrices.
    CMPLX, pointer   :: src_mem_u(:, :)           ! Td-bias. (t, il)
    CMPLX, pointer   :: src_prev(:, :, :, :, :)   ! Source term of previous iteration. (np_intf, ndim, nst, nik, il)
    CMPLX, pointer   :: mem_s(:, :, :, :)         ! the matrices to diagonalize coeff0
    CMPLX, pointer   :: st_intface(:, :, :, :, :) ! The mempry. (np_intf, nst, nik, nleads, max_iter)
  end type ob_terms_t
  
contains


end module ob_rti_m
