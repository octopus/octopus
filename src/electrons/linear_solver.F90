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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module linear_solver_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use lalg_basic_oct_m
  use linear_response_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use preconditioners_oct_m
  use smear_oct_m
  use solvers_oct_m
  use space_oct_m
  use states_elec_oct_m
  use wfs_elec_oct_m

  implicit none

  private

  public ::                              &
       linear_solver_t,                  &
       linear_solver_init,               &
       linear_solver_end,                &
       dlinear_solver_solve_HXeY,        & 
       zlinear_solver_solve_HXeY,        &
       dlinear_solver_solve_HXeY_batch,  & 
       zlinear_solver_solve_HXeY_batch,  &
       linear_solver_ops_per_iter,       &
       linear_solver_obsolete_variables

  type linear_solver_t
    private
    integer,                public :: solver
    type(preconditioner_t), public :: pre
    integer                        :: max_iter
    type(multigrid_t), allocatable :: mgrid
  end type linear_solver_t

  type(profile_t), save :: prof, prof_batch

  type linear_solver_args_t
    private
    type(namespace_t),        pointer :: namespace
    type(linear_solver_t),    pointer :: ls
    type(hamiltonian_elec_t), pointer :: hm
    type(mesh_t),             pointer :: mesh
    type(states_elec_t),      pointer :: st
    integer                           :: ist
    integer                           :: ik
    FLOAT                             :: dshift
    CMPLX                             :: zshift
  end type linear_solver_args_t

  type(linear_solver_args_t) :: args

contains

  ! ---------------------------------------------------------
  subroutine linear_solver_init(this, namespace, gr, states_are_real, mc, space)
    type(linear_solver_t),  intent(out)   :: this
    type(namespace_t),      intent(in)    :: namespace
    type(grid_t),           intent(inout) :: gr
    logical,                intent(in)    :: states_are_real !< for choosing solver
    type(multicomm_t),      intent(in)    :: mc
    type(space_t),          intent(in)    :: space

    integer :: fsolver
    integer :: defsolver

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
    !% Quasi-minimal residual solver, using the symmetrized form <math>A^\dagger A x = A^\dagger y</math> instead of
    !% <math>A x = y</math>. Reliable but very slow. General matrices.
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
    !%Option idrs 11
    !% This is the "Induced Dimension Reduction", IDR(s) (for s=4). IDR(s) is a robust and efficient short recurrence 
    !% Krylov subspace method for solving large nonsymmetric systems of linear equations. It is described in 
    !% [Peter Sonneveld and Martin B. van Gijzen, SIAM J. Sci. Comput. 31, 1035 (2008)]. We have adapted the code
    !% released by M. B. van Gizjen [http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html].
    !%End

    if(conf%devel_version) then
      defsolver = OPTION__LINEARSOLVER__QMR_DOTP
    else
      if(states_are_real) then
        defsolver = OPTION__LINEARSOLVER__QMR_SYMMETRIC
        ! in this case, it is equivalent to LS_QMR_DOTP
      else
        defsolver = OPTION__LINEARSOLVER__QMR_SYMMETRIZED
      end if
    end if

    call parse_variable(namespace, "LinearSolver", defsolver, fsolver)

    ! set up pointer for dot product and norm in QMR solvers
    call mesh_init_mesh_aux(gr%mesh)

    !the last 2 digits select the linear solver
    this%solver = mod(fsolver, 100)

    call preconditioner_init(this%pre, namespace, gr, mc, space)

    !%Variable LinearSolverMaxIter
    !%Type integer
    !%Default 1000
    !%Section Linear Response::Solver
    !%Description
    !% Maximum number of iterations the linear solver does, even if
    !% convergence is not achieved.
    !%End
    call parse_variable(namespace, "LinearSolverMaxIter", 1000, this%max_iter)

    write(message(1),'(a)') 'Linear Solver'
    call messages_print_stress(stdout, trim(message(1)), namespace=namespace)
    
    ! solver 
    select case(this%solver)
      case(OPTION__LINEARSOLVER__CG)
        message(1)='Linear Solver: Conjugate Gradients'

      case(OPTION__LINEARSOLVER__BICGSTAB)
        message(1)='Linear Solver: Biconjugate Gradients Stabilized'

      case(OPTION__LINEARSOLVER__IDRS)
        message(1)='Linear Solver: IDRS'

      case(OPTION__LINEARSOLVER__MULTIGRID)
        message(1)='Multigrid (currently only Gauss-Jacobi - EXPERIMENTAL)'

      case(OPTION__LINEARSOLVER__QMR_SYMMETRIC)
        message(1)='Linear Solver: Quasi-Minimal Residual, for symmetric matrix'

      case(OPTION__LINEARSOLVER__QMR_SYMMETRIZED)
        message(1)='Linear Solver: Quasi-Minimal Residual, for symmetrized matrix'

      case(OPTION__LINEARSOLVER__QMR_DOTP)
        message(1)='Linear Solver: Quasi-Minimal Residual, symmetric with conjugated dot product'

      case(OPTION__LINEARSOLVER__QMR_GENERAL)
        message(1)='Linear Solver: Quasi-Minimal Residual, general algorithm'

      case(OPTION__LINEARSOLVER__SOS)
        message(1)='Linear Solver: Sum-over-States'
    end select

    call messages_info(1)
    
    call messages_print_stress(stdout, namespace=namespace)

    if (this%solver == OPTION__LINEARSOLVER__MULTIGRID) then
      call messages_experimental("Multigrid linear solver")

      SAFE_ALLOCATE(this%mgrid)
      call multigrid_init(this%mgrid, namespace, space, gr%cv, gr%mesh, gr%der, gr%stencil, mc)
    end if

    if (this%solver == OPTION__LINEARSOLVER__QMR_DOTP) then
      call messages_experimental("QMR solver (symmetric with conjugated dot product)")
    end if

    POP_SUB(linear_solver_init)
  end subroutine linear_solver_init

  ! ---------------------------------------------------------
  subroutine linear_solver_end(this)
    type(linear_solver_t), intent(inout) :: this
    this%solver = -1

    call preconditioner_end(this%pre)
    if (allocated(this%mgrid)) then
      call multigrid_end(this%mgrid)
      SAFE_DEALLOCATE_A(this%mgrid)
    end if

  end subroutine linear_solver_end


  ! ---------------------------------------------------------
  integer function linear_solver_ops_per_iter(this) result(n)
    type(linear_solver_t), intent(inout) :: this
    
    select case(this%solver)
    case(OPTION__LINEARSOLVER__BICGSTAB)
      n = 2
    case default ! LS_CG, LS_MULTIGRID, LS_QMR, LS_SOS
      n = 1
    end select
  
  end function linear_solver_ops_per_iter

  ! ----------------------------------------------------------
  
  subroutine linear_solver_obsolete_variables(namespace, old_prefix, new_prefix)
    type(namespace_t),   intent(in)    :: namespace
    character(len=*),    intent(in)    :: old_prefix
    character(len=*),    intent(in)    :: new_prefix
    
    call messages_obsolete_variable(namespace, trim(old_prefix)//"LinearSolver", trim(new_prefix)//"LinearSolver")
    call messages_obsolete_variable(namespace, trim(old_prefix)//"LinearSolverMaxIter", trim(new_prefix)//"LinearSolverMaxIter")

    call preconditioner_obsolete_variables(namespace, old_prefix, new_prefix)

  end subroutine linear_solver_obsolete_variables

#include "undef.F90"

#include "real.F90"

#include "linear_solver_inc.F90"

#include "undef.F90"

#include "complex.F90"
#include "linear_solver_inc.F90"



end module linear_solver_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
