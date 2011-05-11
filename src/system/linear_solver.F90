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
!! $Id: linear_response.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"

module linear_solver_m
  use datasets_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_basic_m
  use linear_response_m
  use loct_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use preconditioners_m
  use scf_tol_m
  use smear_m
  use solvers_m
  use states_m

  implicit none

  private

  integer, public, parameter ::&
       LS_CG              = 5,  &
       LS_BICGSTAB        = 4,  &
       LS_MULTIGRID       = 7,  &
       LS_QMR_SYMMETRIC   = 81, &
       LS_QMR_SYMMETRIZED = 82, &
       LS_QMR_DOTP        = 83, &
       LS_QMR_GENERAL     = 84, &
       LS_SOS             = 9

  public :: &
       linear_solver_t,       &
       linear_solver_init,    &
       linear_solver_end,     &
       dsolve_HXeY,           & 
       zsolve_HXeY,           &
       linear_solver_ops_per_iter,       &
       linear_solver_obsolete_variables

  type linear_solver_t
     integer                :: solver         
     type(preconditioner_t) :: pre
     FLOAT                  :: abs_psi
     integer                :: iter
     integer                :: max_iter
  end type linear_solver_t

  type(profile_t), save :: prof

  type linear_solver_args_t
    type(linear_solver_t), pointer :: ls
    type(hamiltonian_t),   pointer :: hm
    type(grid_t),          pointer :: gr
    type(states_t),        pointer :: st
    integer                        :: ist
    integer                        :: ik
    FLOAT                          :: dshift
    CMPLX                          :: zshift
  end type linear_solver_args_t

  type(linear_solver_args_t) :: args

contains

  ! ---------------------------------------------------------
  subroutine linear_solver_init(this, gr, prefix, def_solver, tol_scheme)
    type(linear_solver_t),  intent(out)   :: this
    type(grid_t),           intent(inout) :: gr
    character(len=*),       intent(in)    :: prefix
    integer, optional,      intent(in)    :: def_solver
    integer, optional,      intent(in)    :: tol_scheme

    integer :: fsolver
    integer :: defsolver_ 
    PUSH_SUB(linear_solver_init)

    !%Variable LinearSolver
    !%Type integer
    !%Default qmr_symmetric
    !%Section Linear Response::Solver
    !%Description
    !% Method for solving linear equations, which occur for Sternheimer linear
    !% response and OEP. The solvers vary in speed, reliability (ability to
    !% converge), and domain of applicability. QMR solvers are most reliable.
    !%Option bicgstab 4
    !% Biconjugate gradients stabilized. Slower than <tt>cg</tt>, but more reliable.
    !% General matrices.
    !%Option cg 5
    !% Conjugate gradients. Fast but unreliable. Hermitian matrices only
    !% (no eta in Sternheimer).
    !%Option multigrid 7
    !% Multigrid solver, currently only Gauss-Jacobi (experimental).
    !% Slow, but fairly reliable. General matrices.
    !%Option qmr_symmetric 81
    !% Quasi-minimal residual solver, for (complex) symmetric matrices. [Real symmetric
    !% is equivalent to Hermitian.] Slightly slower than <tt>bicgstab</tt> but more reliable.
    !% For Sternheimer, must be real wavefunctions, but can have eta.
    !%Option qmr_symmetrized 82
    !% Quasi-minimal residual solver, using the symmetrized form A^T A x = A^T y instead of
    !% A x = y. Reliable but very slow. General matrices.
    !%Option qmr_dotp 83
    !% Quasi-minimal residual solver, for Hermitian matrices, using the
    !% symmetric algorithm with conjugated dot product (experimental). Slightly slower than <tt>bicgstab</tt>
    !% but more reliable. Can always be used in Sternheimer.
    !%Option qmr_general 84
    !% Quasi-minimal residual solver, for general matrices, using the
    !% most general form of the algorithm. Slow and unreliable.
    !%Option sos 9
    !% Sum over states: the Sternheimer equation is solved by using
    !% the explicit solution in terms of the ground-state
    !% wavefunctions. You need unoccupied states to use this method.
    !% Unlike the other methods, may not give the correct answer.
    !%End

    defsolver_ = LS_QMR_SYMMETRIZED
    if(present(def_solver)) defsolver_ = def_solver

    if (parse_isdef(datasets_check(trim(prefix)//"LinearSolver")) /= 0 ) then 
      call parse_integer  (datasets_check(trim(prefix)//"LinearSolver"), defsolver_, fsolver)
    else
      call parse_integer  (datasets_check("LinearSolver"), defsolver_, fsolver)
    end if

    !the last 2 digits select the linear solver
    this%solver = mod(fsolver, 100)

    call preconditioner_init(this%pre, gr, prefix)

    !%Variable LinearSolverMaxIter
    !%Type integer
    !%Default 1000
    !%Section Linear Response::Solver
    !%Description
    !% Maximum number of iterations the linear solver does, even if
    !% convergence is not achieved.
    !%End
    if (parse_isdef(datasets_check(trim(prefix)//"LinearSolverMaxIter")) /= 0) then 
      call parse_integer  (datasets_check(trim(prefix)//"LinearSolverMaxIter"), 1000, this%max_iter)
    else
      call parse_integer  (datasets_check("LinearSolverMaxIter"), 1000, this%max_iter)
    end if

    write(message(1),'(a)') 'Linear Solver'
    call messages_print_stress(stdout, trim(message(1)))
    
    ! solver 
    select case(this%solver)
      case(LS_CG)
        message(1)='Linear Solver: Conjugate Gradients'

      case(LS_BICGSTAB)
        message(1)='Linear Solver: Biconjugate Gradients Stabilized'

      case(LS_MULTIGRID)
        message(1)='Multigrid (currently only Gauss-Jacobi - EXPERIMENTAL)'

      case(LS_QMR_SYMMETRIC)
        message(1)='Linear Solver: Quasi-Minimal Residual, for symmetric matrix'

      case(LS_QMR_SYMMETRIZED)
        message(1)='Linear Solver: Quasi-Minimal Residual, for symmetrized matrix'

      case(LS_QMR_DOTP)
        message(1)='Linear Solver: Quasi-Minimal Residual, symmetric with conjugated dot product'

      case(LS_QMR_GENERAL)
        message(1)='Linear Solver: Quasi-Minimal Residual, general algorithm'

      case(LS_SOS)
        message(1)='Linear Solver: Sum-over-States'
    end select

    call messages_info(1)
    
    call messages_print_stress(stdout)

    if(this%solver == LS_MULTIGRID) call messages_experimental("Multigrid linear solver")
    if(this%solver == LS_QMR_DOTP)  call messages_experimental("QMR solver (symmetric with conjugated dot product)")

    POP_SUB(linear_solver_init)

  end subroutine linear_solver_init

  ! ---------------------------------------------------------
  subroutine linear_solver_end(this)
    type(linear_solver_t), intent(inout) :: this
    this%solver = -1

    call preconditioner_end(this%pre)

  end subroutine linear_solver_end


  ! ---------------------------------------------------------
  integer function linear_solver_ops_per_iter(this) result(n)
    type(linear_solver_t), intent(inout) :: this
    
    select case(this%solver)
    case(LS_BICGSTAB)
      n = 2
    case default ! LS_CG, LS_MULTIGRID, LS_QMR, LS_SOS
      n = 1
    end select
  
  end function linear_solver_ops_per_iter

  ! ----------------------------------------------------------
  
  subroutine linear_solver_obsolete_variables(old_prefix, new_prefix)
    character(len=*),    intent(in)    :: old_prefix
    character(len=*),    intent(in)    :: new_prefix
    
    call messages_obsolete_variable(trim(old_prefix)//"LinearSolver", trim(new_prefix)//"LinearSolver")
    call messages_obsolete_variable(trim(old_prefix)//"LinearSolverMaxIter", trim(new_prefix)//"LinearSolverMaxIter")

    call preconditioner_obsolete_variables(old_prefix, new_prefix)

  end subroutine linear_solver_obsolete_variables

#include "undef.F90"

#include "real.F90"

#include "linear_solver_inc.F90"

#include "undef.F90"

#include "complex.F90"
#include "linear_solver_inc.F90"

end module linear_solver_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
