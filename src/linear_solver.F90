!! Copyright (C) 2004 E.S. Kadantsev, M. Marques
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
!! $Id: linear_response.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"

module linear_solver_m
  use datasets_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lib_basic_alg_m
  use lib_oct_m
  use lib_oct_parser_m
  use mesh_m
  use messages_m
  use preconditioners_m
  use states_m

  implicit none

  private

  integer, public, parameter ::&
       LS_HX_FIXED = 11, & 
       LS_HX = 10,       &
       LS_CG = 5,        &
       LS_BICGSTAB = 3

  public :: &
       linear_solver_t, &
       linear_solver_init, &
       linear_solver_end, &
       dsolve_HXeY, & 
       zsolve_HXeY, &
       linear_solver_ops_per_iter
  

  type linear_solver_t
     integer :: solver         
     type(preconditioner_t) :: pre
     FLOAT   :: abs_psi
     FLOAT   :: initial_tol
     FLOAT   :: final_tol
     FLOAT   :: tol
     integer :: iter
     integer :: max_iter
  end type linear_solver_t

contains

  ! ---------------------------------------------------------
  subroutine linear_solver_init(this, gr, prefix, def_solver)
    type(linear_solver_t),        intent(out)   :: this
    type(grid_t),      intent(inout) :: gr
    character(len=*),  intent(in)    :: prefix
    integer, optional, intent(in)    :: def_solver

    integer :: fsolver

    call push_sub('linear_solver.linear_solver_init')

    !%Variable LinearSolver
    !%Type integer
    !%Default cg
    !%Section Linear Response::Solver
    !%Description
    !% To calculate response using density functional perturbation
    !% theory is necessary to solve the sterheimer equation, this is a self
    !% consistent linear equation where the operator is the Kohn-Sham
    !% hamiltonian with a complex shift. This variable selects which
    !% method to use in order to solve this linear equation.
    !%Option hx 10
    !% This solver aproximates the solution of the non-hermitian
    !% equation by a power series in terms of the inverse of the
    !% hamiltonian. Each term implies solving a hermitian linear equation by
    !% conjugated gradients.
    !% Altough this might sound inefficient, the hermitian equations
    !% to solve are better conditioned than the original, so this can be more
    !% efficient than a non-hermitian solver.
    !% This version of the solvers apply the number of terms in the
    !% series as needed to be under the tolerance required.
    !%Option hx_fixed 11
    !% This is the same prevoius solver, but only two steps of the
    !% series are applied.
    !%Option cg 5
    !% Conjugated gradients. This is the fastest solver but does not
    !% work when an imaginary shift is added. This is the default.
    !%Option bicgstab 3
    !% Biconjugated gradients stabilized. This is an improved version
    !% of bcg that is faster and more stable. This is the default when
    !% complex polarizabilities are calculated.
    !%End

    if(present(def_solver)) then
      call loct_parse_int  (check_inp(trim(prefix)//"LinearSolver"), def_solver, fsolver)
    else 
      call loct_parse_int  (check_inp(trim(prefix)//"LinearSolver"), LS_CG, fsolver)
    end if

    !the last 2 digits select the linear solver
    this%solver = mod(fsolver, 100)

    call preconditioner_init(this%pre, gr)

    !%Variable LinearSolverMaxIter
    !%Type integer
    !%Default 1000
    !%Section Linear Response::Solver
    !%Description
    !% Maximum number of iterations the linear solver does, even if
    !% convergency is not achieved.
    !%End

    call loct_parse_int  (check_inp(trim(prefix)//"LinearSolverMaxIter"), 1000, this%max_iter)

    !%Variable LinearSolverInitTol
    !%Type float
    !%Default 1e-2
    !%Section Linear Response::Solver
    !%Description
    !% This is the tolerance to determine that the linear solver has converged.
    !%End
    
    call loct_parse_float(check_inp(trim(prefix)//"LinearSolverInitTol"), &
        CNST(1e-2), this%initial_tol)

    !%Variable LinearSolverTol
    !%Type float
    !%Default 1e-6
    !%Section Linear Response::Solver
    !%Description
    !% This is the tolerance to determine that the linear solver has converged.
    !%End
    
    call loct_parse_float(check_inp(trim(prefix)//"LinearSolverTol"), &
        CNST(1e-6), this%final_tol)
 
    this%tol = this%final_tol
    
    !WRITE INFO
    
    write(message(1),'(a)') 'Linear Solver'
    call messages_print_stress(stdout, trim(message(1)))
    
    
    ! solver 
    select case(this%solver)
      case(LS_CG)
        message(1)='Linear Solver: Conjugated Gradients'

      case(LS_BICGSTAB)
        message(1)='Linear Solver: Biconjugated Gradients Stabilized'

      case(LS_HX_FIXED)
        message(1)='Linear Solver: Fixed Hermitian Expansion'

      case(LS_HX)
        message(1)='Linear Solver: Hermitian Expansion'

    end select

    call write_info(1)
    
    call messages_print_stress(stdout)

    call pop_sub()

  end subroutine linear_solver_init

  ! ---------------------------------------------------------
  subroutine linear_solver_end(this)
    type(linear_solver_t), intent(inout) :: this
    this%solver = -1
  end subroutine linear_solver_end

  integer function linear_solver_ops_per_iter(this) result(n)
    type(linear_solver_t), intent(inout) :: this
    
    select case(this%solver)
      
    case(LS_CG)
      n = 1
    case(LS_HX_FIXED)
      n = 1
    case(LS_HX)
      n = 1
    case(LS_BICGSTAB)
      n = 2
    case default 
      message(1)="Unknown linear response solver"
      call write_fatal(1)
    end select
  
  end function linear_solver_ops_per_iter

#include "undef.F90"

#include "real.F90"

#include "linear_solver_inc.F90"

#include "undef.F90"

#include "complex.F90"
#include "linear_solver_inc.F90"

end module linear_solver_m
