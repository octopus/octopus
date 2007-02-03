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
       zsternheimer_solve, &
       sternheimer_add_fxc, &
       sternheimer_add_hartree
  
  type sternheimer_t
     private
     type(linear_solver_t) :: solver
     type(mix_t) :: mixer
     type(scf_tol_t) :: scftol
     logical :: add_fxc
     logical :: add_hartree
     logical :: ok
     logical :: hermitian 
     logical :: orth_response
 end type sternheimer_t

contains

  subroutine sternheimer_init(this, h, gr, prefix, hermitian)
    type(sternheimer_t), intent(out) :: this
    type(hamiltonian_t), intent(in)  :: h
    type(grid_t),      intent(inout) :: gr
    character(len=*),  intent(in)    :: prefix
    logical, optional, intent(in)    :: hermitian

    integer :: ham_var

    !%Variable OrthResponse
    !%Type logical
    !%Default true
    !%Section Linear Response::Sternheimer
    !%Description
    !% Wheter variations of the wavefunctios should be orthogonalized
    !% or not against the occupied states.
    !%End

    call loct_parse_logical(check_inp('PolOrthResponse'), .true., this%orth_response)

    !%Variable HamiltonianVariation
    !%Type integer
    !%Default hartree+fxc
    !%Section Linear Response::Sternheimer
    !%Description
    !% The terms are considered in the variation of the
    !% hamiltonian. V_ext is always considered. The default is to include
    !% also the exchange, correlation and hartree terms. If you want
    !% to choose the exchange and correlation kernel use the variable
    !% XCKernel.
    !%Option hartree 1
    !% The variation of the hartree potential.
    !%Option fxc 2
    !% The exchange and correlation kernel, the variation of the
    !% exchange and correlation potential.
    !%End

    if(.not. h%ip_app) then 
      call loct_parse_int(check_inp('PolHamiltonianVariation'), 3, ham_var)    
      this%add_fxc = ((ham_var/2) == 1)
      this%add_hartree = (mod(ham_var, 2) == 1)
    else
      this%add_fxc = .false. 
      this%add_hartree = .false.
    end if
    
    message(1) = "Variation of the hamiltonian in Sternheimer equation: V_ext"
    if(this%add_hartree) write(message(1), '(2a)') trim(message(1)), ' + hartree'
    if(this%add_fxc)     write(message(1), '(2a)') trim(message(1)), ' + fxc'
    call write_info(1)

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

  logical function sternheimer_add_fxc(this) result(r)
    type(sternheimer_t) :: this
    r = this%add_fxc
  end function sternheimer_add_fxc

  logical function sternheimer_add_hartree(this) result(r)
    type(sternheimer_t) :: this
    r = this%add_hartree
  end function sternheimer_add_hartree

#include "complex.F90"
#include "sternheimer_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "sternheimer_inc.F90"

end module sternheimer_m
