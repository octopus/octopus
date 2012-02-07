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

#include "global.h"

module ob_terms_m
  use simul_box_m

  implicit none

  private
  public ::     &
    ob_terms_t, &
    ob_memsrc_t

  integer, parameter, public :: &
    MEM_TERM_FLAG  = 1,         &
    SRC_TERM_FLAG  = 2

  !> source and mem-coeffs
  type ob_memsrc_t
    integer          :: sp_length              !< Length of the sparse array (as 1d array).
    CMPLX, pointer   :: q(:, :, :)             !< Memory coefficients, for mem_type=1. (i, j, t)
    CMPLX, pointer   :: p(:, :, :)             !< Memory coefficients, for mem_type=1. (i, j, t)
    CMPLX, pointer   :: q_sp(:, :)             !< Memory coefficients, for mem_type=2. (sp_length, t)
    CMPLX, pointer   :: q_s(:, :, :)           !< The matrices to diagonalize coeff0.
    integer, pointer :: sp2full_map(:)         !< The mapping indices from sparse to full matrices.
    CMPLX, pointer   :: src_prev(:, :, :, :)   !< Source term of previous iteration. (ip_intf, idim, ist, ik)
    CMPLX, pointer   :: st_intface(:, :, :, :) !< The memory. (ip_intf, ist, ik, t)
  end type ob_memsrc_t

  !> Additional source and memory terms.
  type ob_terms_t
    integer           :: mem_type              !< SAVE_CPU_TIME: fast/lots of memory, SAVE_RAM_USAGE: slow/little memory.
    integer           :: additional_terms      !< Shall we add source and memory term?
    integer           :: max_mem_coeffs        !< How many memory coefficients? (for open system)
    type(ob_memsrc_t) :: lead(2*MAX_DIM)       !< source and memory coefficients for each lead
    CMPLX, pointer    :: src_mem_u(:, :)       !< Prefactor for source and memory term. (t, il)
    FLOAT             :: td_pot0(2*MAX_DIM)    !< bias at t=0
  end type ob_terms_t

end module ob_terms_m
