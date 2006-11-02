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
  use datasets_m
  use functions_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lib_basic_alg_m
  use lib_oct_parser_m
  use messages_m
  use nl_operator_m
  use stencil_star_m
  use varinfo_m

  implicit none
  private
  
  integer, public, parameter ::     &
    PRECONDITIONER_NONE      = 0,   &
    PRECONDITIONER_SMOOTHING = 1,   &
    PRECONDITIONER_JACOBI    = 2
  
  public ::                         &
    preconditioner_t,               &
    preconditioner_init,            &
    preconditioner_end,             &
    dpreconditioner_apply,          &
    zpreconditioner_apply
  
  type preconditioner_t
    integer :: which

    type(nl_operator_t) :: op
    FLOAT, pointer      :: diag_lapl(:) ! diagonal of the laplacian
  end type preconditioner_t
  
contains

  ! --------------------------------------------------------- 
  subroutine preconditioner_init(this, gr)
    type(preconditioner_t), intent(out)   :: this 
    type(grid_t),           intent(inout) :: gr
    
    FLOAT, parameter :: alpha = M_HALF
    
    !%Variable Preconditioner
    !%Type integer
    !%Default filter
    !%Section SCF::EigenSolver
    !%Description
    !% Which preconditioner to use in order to solve the Kohn-Sham equations or
    !% the linear-response equations.
    !%Option no 0
    !% Do not apply preconditioner.
    !%Option filter 1
    !% Filter preconditioner.
    !%Option jacobi 2
    !% Jacobi preconditioner. only the local part of the pseudopotential is used.
    !%End
    call loct_parse_int(check_inp('Preconditioner'), PRECONDITIONER_SMOOTHING, this%which)
    if(.not.varinfo_valid_option('Preconditioner', this%which)) call input_error('Preconditioner')
    call messages_print_var_option(stdout, "Preconditioner", this%which)
    

    select case(this%which)
    case(PRECONDITIONER_SMOOTHING)
      ! the smoothing has a star stencil like the laplacian
      call nl_operator_init(this%op, 2*NDIM + 1)
      call stencil_star_get_lapl(NDIM, 1, this%op%stencil)
      call nl_operator_build(gr%m, this%op, NP, const_w = .true.)
      
      this%op%w_re(1, 1) = alpha
      this%op%w_re(2:,1) = M_HALF*(M_ONE-alpha)/NDIM

    case(PRECONDITIONER_JACOBI)
      ALLOCATE(this%diag_lapl(NP), NP)
      call f_laplacian_diag (gr%sb, gr%f_der, this%diag_lapl)
      call lalg_scal(NP, -M_HALF, this%diag_lapl(:))
    end select

  end subroutine preconditioner_init

  
  ! ---------------------------------------------------------
  subroutine preconditioner_end(this)
    type(preconditioner_t), intent(inout) :: this 

    select case(this%which)
    case(PRECONDITIONER_SMOOTHING)
      call nl_operator_end(this%op)

    case(PRECONDITIONER_JACOBI)
      deallocate(this%diag_lapl)
    end select

  end subroutine preconditioner_end


#include "undef.F90"
#include "complex.F90"
#include "preconditioners_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "preconditioners_inc.F90"

end module preconditioners_m
