!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: em_resp.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"

module sternheimer_m
  use datasets_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lib_basic_alg_m
  use lib_oct_parser_m
  use linear_solver_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mix_m
  use output_m
  use poisson_m
  use restart_m
  use scf_tol_m
  use states_m
  use string_m
  use system_m
  use units_m

  implicit none

  private
  public :: &
       sternheimer_t,      &
       sternheimer_init,   &
       sternheimer_end,    &
       dsternheimer_solve, & 
       zsternheimer_solve
  
  type sternheimer_t
     private
     type(linear_solver_t) :: solver
     type(mix_t) :: mixer
     type(scf_tol_t) :: scftol
     logical :: add_fxc
     logical :: add_hartree
     logical :: from_scratch
     logical :: ok
     logical :: hermitian 
     logical  :: orth_response
 end type sternheimer_t

contains

  subroutine sternheimer_init(this, gr, prefix, hermitian)
    type(sternheimer_t), intent(out) :: this
    type(grid_t),      intent(inout) :: gr
    character(len=*),  intent(in)    :: prefix
    logical, optional, intent(in)    :: hermitian
    
    this%add_fxc = .true.
    this%add_hartree = .true.
    this%orth_response = .true. 
    this%from_scratch = .false.

    if(present(hermitian)) then 
      if(.not. hermitian) then 
        call linear_solver_init(this%solver, gr, prefix, def_solver=LS_BICGSTAB)
      else
        call linear_solver_init(this%solver, gr, prefix, def_solver=LS_CG)
      end if
    else
      call linear_solver_init(this%solver, gr, prefix)
    end if

    call scf_tol_init(this%scftol, prefix)

  end subroutine sternheimer_init

  subroutine sternheimer_end(this)
    type(sternheimer_t) :: this

    call linear_solver_end(this%solver)
    call scf_tol_end(this%scftol)

  end subroutine sternheimer_end

#include "complex.F90"
#include "sternheimer_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "sternheimer_inc.F90"

end module sternheimer_m
