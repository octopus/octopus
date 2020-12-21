!! Copyright (C) 2005-2006 Heiko Appel
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

module root_solver_oct_m
  use global_oct_m
  use lalg_adv_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                               &
    root_solver_t,                        &
    root_solver_init,                     &
    root_solver_read,                     &
    droot_solver_run

  integer, public, parameter ::           &
    ROOT_NEWTON      =  3

  type root_solver_t
    private
    integer :: solver_type    !< what solver to use (see ROOT_* variables above)_oct_m
    integer :: dim            !< dimensionality of the problem
    integer :: maxiter        !< maximal number of iterations
    integer :: usediter       !< number of iterations actually performed
    FLOAT   :: abs_tolerance
    FLOAT   :: rel_tolerance
  end type root_solver_t

contains

  ! ---------------------------------------------------------
  subroutine root_solver_init(rs, namespace, dimensionality, solver_type, maxiter, &
    rel_tolerance, abs_tolerance)
    type(root_solver_t), intent(out) :: rs
    type(namespace_t),   intent(in)  :: namespace
    integer,             intent(in)  :: dimensionality
    integer, optional,   intent(in)  :: solver_type, maxiter
    FLOAT, optional,     intent(in)  :: rel_tolerance, abs_tolerance

    ! no push_sub, called too often

    ! Set dimension
    rs%dim = dimensionality

    ! Read config parameters for the solver
    call root_solver_read(rs, namespace)

    if(present(solver_type))     rs%solver_type     = solver_type
    if(present(maxiter))         rs%maxiter         = maxiter
    if(present(rel_tolerance))   rs%rel_tolerance   = rel_tolerance
    if(present(abs_tolerance))   rs%abs_tolerance   = abs_tolerance

  end subroutine root_solver_init


  ! ---------------------------------------------------------
  subroutine root_solver_read(rs, namespace)
    type(root_solver_t),   intent(inout) :: rs
    type(namespace_t),     intent(in)  :: namespace

    PUSH_SUB(root_solver_read)

    !%Variable RootSolver
    !%Type integer
    !%Default root_newton
    !%Section Math::RootSolver
    !%Description
    !% Specifies what kind of root solver will be used.
    !%Option root_newton 3
    !% Newton method.
    !%End
    call parse_variable(namespace, 'RootSolver', ROOT_NEWTON, rs%solver_type)
    if (.not. varinfo_valid_option('RootSolver', rs%solver_type)) then
      call messages_input_error(namespace, 'RootSolver')
    end if

    !%Variable RootSolverMaxIter
    !%Type integer
    !%Default 100
    !%Section Math::RootSolver
    !%Description
    !% In case of an iterative root solver, this variable determines the maximum number
    !% of iteration steps.
    !%End
    call parse_variable(namespace, 'RootSolverMaxIter', 500, rs%maxiter)

    !%Variable RootSolverRelTolerance
    !%Type float
    !%Default 1e-8
    !%Section Math::RootSolver
    !%Description
    !% Relative tolerance for the root-finding process.
    !%End
    call parse_variable(namespace, 'RootSolverRelTolerance', CNST(1e-10), rs%rel_tolerance)

    !%Variable RootSolverAbsTolerance
    !%Type float
    !%Default 1e-8
    !%Section Math::RootSolver
    !%Description
    !% Relative tolerance for the root-finding process.
    !%End
    call parse_variable(namespace, 'RootSolverAbsTolerance', CNST(1e-10), rs%abs_tolerance)

    POP_SUB(root_solver_read)
  end subroutine root_solver_read

  ! ---------------------------------------------------------
  subroutine droot_solver_run(rs, func, root, success, startval)
    type(root_solver_t), intent(in)  :: rs
    FLOAT,               intent(out) :: root(:)        !< roots we are searching
    logical,             intent(out) :: success
    FLOAT, optional,     intent(in)  :: startval(:)    !< start value for the search
    interface
      subroutine func(z, f, jf)
        implicit none
        FLOAT, intent(in)  :: z(:)
        FLOAT, intent(out) :: f(:), jf(:, :)
      end subroutine func
    end interface

    ! no push_sub, called too often

    ! Initializations
    root = M_ZERO
    success = .false.

    select case(rs%solver_type)
    case(ROOT_NEWTON)
      call droot_newton(rs, func, root, startval, success)
    case default
      write(message(1), '(a,i4,a)') "Error in droot_solver_run: '", rs%solver_type, "' is not a valid root solver"
      call messages_fatal(1)
    end select

  end subroutine

  ! ---------------------------------------------------------
  !> Newton-Raphson scheme can only be used in the real case.
  subroutine droot_newton(rs, func, root, startval, success)
    type(root_solver_t), intent(in)  :: rs
    FLOAT,               intent(out) :: root(:)        !< root we are searching
    FLOAT,               intent(in)  :: startval(:)    !< start value for the search
    logical,             intent(out) :: success
    interface
      subroutine func(z, f, jf)
        implicit none
        FLOAT, intent(in)  :: z(:)
        FLOAT, intent(out) :: f(:), jf(:, :)
      end subroutine func
    end interface

    integer :: iter
    FLOAT   :: err
    FLOAT, allocatable :: f(:), jf(:, :), delta(:, :), rhs(:, :)

    ! no push_sub, called too often

    SAFE_ALLOCATE(    f(1:rs%dim))
    SAFE_ALLOCATE(   jf(1:rs%dim, 1:rs%dim))
    SAFE_ALLOCATE(delta(1:rs%dim, 1))
    SAFE_ALLOCATE(  rhs(1:rs%dim, 1))

    root = startval
    call func(root, f, jf)
    err = sum(f(1:rs%dim)**2)

    success = .true.
    iter = 0
    do while(err > rs%abs_tolerance)
      rhs(1:rs%dim, 1) = -f(1:rs%dim)
      call lalg_linsyssolve(rs%dim, 1, jf, rhs, delta)
      root(1:rs%dim) = root(1:rs%dim) + delta(1:rs%dim, 1)
      iter = iter + 1
      if(iter > rs%maxiter) then
        success = .false.
        exit
      end if
      call func(root, f, jf)
      err = sum(f(1:rs%dim)**2)
    end do

    SAFE_DEALLOCATE_A(f)
    SAFE_DEALLOCATE_A(jf)
    SAFE_DEALLOCATE_A(delta)
    SAFE_DEALLOCATE_A(rhs)

  end subroutine droot_newton

end module root_solver_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
