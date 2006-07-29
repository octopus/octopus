!! Copyright (C) 2006 X. Andrade, M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module preconditioners_m
  use global_m
  use grid_m
  use nl_operator_m
  use stencil_star_m

  implicit none
  private
  
  integer, public, parameter ::     &
    PRECONDITIONER_SMOOTHING = 1001
  
  public ::                         &
    init_preconditioner,            &
    end_preconditioner,             &
    dapply_preconditioner,          &
    zapply_preconditioner,          &
    preconditioner_smoothing_t
  
  interface init_preconditioner
    module procedure init_preconditioner_smoothing
  end interface
  
  interface end_preconditioner
    module procedure end_preconditioner_smoothing
  end interface

  interface dapply_preconditioner
    module procedure dapply_precond_smoothing
    module procedure dapply_precond_smoothing_wfs
  end interface

  interface zapply_preconditioner
    module procedure zapply_precond_smoothing
    module procedure zapply_precond_smoothing_wfs
  end interface

  type preconditioner_smoothing_t
    type(nl_operator_t) :: op
  end type preconditioner_smoothing_t
  
contains

  ! --------------------------------------------------------- 
  subroutine init_preconditioner_smoothing(this, gr)
    type(preconditioner_smoothing_t), intent(inout) :: this 
    type(grid_t),        intent(in)    :: gr
    
    FLOAT, parameter :: alpha = M_HALF
    
    ! the smoothing has a star stencil like the laplacian
    call nl_operator_init(this%op, 2*NDIM + 1)
    call stencil_star_get_lapl(NDIM, 1, this%op%stencil)
    call nl_operator_build(gr%m, this%op, NP, const_w = .true.)
    
    this%op%w_re(1, 1) = alpha
    this%op%w_re(2:,1) = M_HALF*(M_ONE-alpha)/NDIM

  end subroutine init_preconditioner_smoothing

  subroutine end_preconditioner_smoothing(this)
    type(preconditioner_smoothing_t), intent(inout) :: this 

    call nl_operator_end(this%op)
  end subroutine end_preconditioner_smoothing

#include "undef.F90"
#include "complex.F90"

#include "preconditioners_inc.F90"

#include "undef.F90"
#include "real.F90"

#include "preconditioners_inc.F90"

end module preconditioners_m
