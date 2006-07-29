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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module eigen_solver_m
  use global_m
  use messages_m
  use datasets_m
  use profiling_m
  use lib_oct_parser_m
  use lib_basic_alg_m
  use lib_adv_alg_m
  use math_m, only : sort
  use io_m
  use nl_operator_m
  use stencil_star_m
  use mesh_function_m
  use mesh_m
  use grid_m
  use functions_m
  use states_m
  use hamiltonian_m
  use preconditioners_m

  implicit none

  private
  public ::             &
    eigen_solver_t,     &
    eigen_solver_init,  &
    eigen_solver_end,   &
    eigen_solver_run

  private ::            &
    dmv,                &
    zmv

  type eigen_solver_t
    integer :: es_type    ! which eigen solver to use

    FLOAT :: init_tol
    FLOAT :: final_tol
    integer  :: final_tol_iter
    integer  :: es_maxiter

    integer :: arnoldi_vectors
    FLOAT :: imag_time

    logical :: do_subspace_diag

    ! Stores information about how well it performed.
    FLOAT, pointer :: diff(:, :)
    integer           :: matvec
    integer           :: converged
  end type eigen_solver_t


#ifdef HAVE_TRLAN
  integer, parameter :: RS_LANCZOS = 1
#endif
  integer, parameter :: RS_PLAN    = 11
#if defined(HAVE_ARPACK)
  integer, parameter :: ARPACK     = 3
#endif
#if defined(HAVE_JDQZ)
  integer, parameter :: JDQZ       = 4
#endif
  integer, parameter :: RS_CG      = 5
  integer, parameter :: RS_CG_NEW  = 6
  integer, parameter :: EVOLUTION  = 7

  ! For the TRLan package
#ifdef HAVE_TRLAN
  integer                      :: ik_trlan
  type(hamiltonian_t), pointer :: h_trlan
  type(mesh_t),        pointer :: m_trlan
  type(states_t),      pointer :: st_trlan
#endif

  type(hamiltonian_t), pointer :: h_
  type(grid_t),        pointer :: gr_
  integer                         :: ik_

  type(preconditioner_smoothing_t) :: filter

contains

  ! ---------------------------------------------------------
  subroutine eigen_solver_init(gr, eigens, st, max_iter_default)
    type(grid_t),         intent(in)  :: gr
    type(eigen_solver_t), intent(out) :: eigens
    type(states_t),       intent(in)  :: st
    integer,              intent(in)  :: max_iter_default

    call push_sub('eigen.eigen_solver_init')

    !%Variable EigenSolver
    !%Type integer
    !%Default cg
    !%Section SCF::EigenSolver
    !%Description
    !% Decides the eigensolver that obtains the lowest eigenvalues
    !% and eigenfunctions of the Kohn-Sham Hamiltonian.
    !%Option cg 5
    !% Conjugate-gradients algorithm
    !%Option trlan 1
    !% Lanczos scheme. Requires the TRLan package.
    !%Option plan 11
    !% Preconditioned Lanczos scheme.
    !%Option arpack 3
    !% Implicitly Restarted Arnoldi Method. Requires the ARPACK package
    !%Option jdqz 5
    !% Jacobi-Davidson scheme. Requires the JDQZ package
    !%Option cg_new 6
    !% A rewritting of the cg option, that will eventually substitute it.
    !%Option evolution 7
    !% Propagation in imaginary time
    !%End
    call loct_parse_int(check_inp('EigenSolver'), RS_CG, eigens%es_type)

    select case(eigens%es_type)
    case(RS_CG_NEW)
    case(RS_CG)
#ifdef HAVE_TRLAN
    case(RS_LANCZOS)
#endif
    case(RS_PLAN)
#if defined(HAVE_ARPACK)
    case(ARPACK)

      !%Variable EigenSolverArnoldiVectors
      !%Type integer
      !%Default 20
      !%Section SCF::EigenSolver
      !%Description
      !% This indicates how many Arnoldi vectors are generated
      !% It must satisfy EigenSolverArnoldiVectors - Number Of Eigenvectors >= 2.
      !% See the ARPACK documentation for more details.
      !%End
      call loct_parse_int(check_inp('EigenSolverArnoldiVectors'), 20, eigens%arnoldi_vectors)
      if(eigens%arnoldi_vectors-st%nst < 2) call input_error('EigenSolverArnoldiVectors')
#endif
#if defined(HAVE_JDQZ)
    case(JDQZ)
#endif
    case(EVOLUTION)
      !%Variable EigenSolverImaginaryTime
      !%Type float
      !%Default 1.0
      !%Section SCF::EigenSolver
      !%Description
      !% The imaginary-time step that is used in the imaginary-time evolution
      !% method to obtain the lowest eigenvalues/eigenvectors.
      !% It must satisfy EigenSolverImaginaryTime > 0.
      !%End
      call loct_parse_float(check_inp('EigenSolverImaginaryTime'), M_ONE, eigens%imag_time)
      if(eigens%imag_time <= M_ZERO) call input_error('EigenSolverImaginaryTime')
    case default
      call input_error('EigenSolver')
    end select
    call messages_print_var_option(stdout, "EigenSolver", eigens%es_type)

    !%Variable EigenSolverInitTolerance
    !%Type float
    !%Default 1.0e-6
    !%Section SCF::EigenSolver
    !%Description
    !% This is the initial tolerance for the eigenvectors.
    !%End
    call loct_parse_float(check_inp('EigenSolverInitTolerance'), CNST(1.0e-6), eigens%init_tol)
    if(eigens%init_tol < 0) call input_error('EigenSolverInitTolerance')

    !%Variable EigenSolverFinalTolerance
    !%Type float
    !%Default 1.0e-6
    !%Section SCF::EigenSolver
    !%Description
    !% This is the final tolerance for the eigenvectors. Must be smaller than <tt>EigenSolverInitTolerance</tt>.
    !%End
    call loct_parse_float(check_inp('EigenSolverFinalTolerance'), CNST(1.0e-6), eigens%final_tol)
    if(eigens%final_tol < 0 .or. eigens%final_tol > eigens%init_tol) call input_error('EigenSolverFinalTolerance')

    !%Variable EigenSolverFinalToleranceIteration
    !%Type integer
    !%Default 7
    !%Section SCF::EigenSolver
    !%Description
    !% Determines how many interactions are needed 
    !% to go from <tt>EigenSolverInitTolerance</tt> to <tt>EigenSolverFinalTolerance</tt>.
    !% Must be larger than 1.
    !%End
    call loct_parse_int(check_inp('EigenSolverFinalToleranceIteration'), 7, eigens%final_tol_iter)
    if(eigens%final_tol_iter <= 1) call input_error('EigenSolverFinalToleranceIteration')

    !%Variable EigenSolverMaxIter
    !%Type integer
    !%Default 25
    !%Section SCF::EigenSolver
    !%Description
    !% Determines the maximum number of iterations 
    !% for the eigensolver (per state) -- that is, if this number is reached, the diagonalization
    !% is stopped even if the desired tolerance was not achieved. Must be larger or equal than 1.
    !%End
    call loct_parse_int(check_inp('EigenSolverMaxIter'), max_iter_default, eigens%es_maxiter)
    if(eigens%es_maxiter < 1) call input_error('EigenSolverMaxIter')

    !%Variable EigenSolverSubspaceDiag
    !%Type logical
    !%Default false
    !%Section SCF::EigenSolver
    !%Description
    !% When true, a Rayleigh-Ritz projection is performed after the eigensolver routine. 
    !% This is a diagonalization of the Hamiltonian in the subspace spanned by the eigenvectors 
    !% that separates eigenvectors that have close eigenvalues.
    !%End
    call loct_parse_logical(check_inp('EigenSolverSubspaceDiag'), .false., eigens%do_subspace_diag)
    
    select case(eigens%es_type)
    case(RS_PLAN)
      call init_preconditioner(filter, gr)
    end select

    nullify(eigens%diff)
    ALLOCATE(eigens%diff(st%nst, st%d%nik), st%nst*st%d%nik)

    eigens%converged = 0
    eigens%matvec    = 0

    call pop_sub()

  end subroutine eigen_solver_init


  ! ---------------------------------------------------------
  subroutine eigen_solver_end(eigens)
    type(eigen_solver_t), intent(inout) :: eigens

    select case(eigens%es_type)
    case(RS_PLAN)
      call end_preconditioner(filter)
    end select

    deallocate(eigens%diff)
    nullify(eigens%diff)

  end subroutine eigen_solver_end


  ! ---------------------------------------------------------
  subroutine eigen_solver_run(eigens, gr, st, h, iter, conv, verbose)
    type(eigen_solver_t), intent(inout) :: eigens
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(inout) :: h
    integer,              intent(in)    :: iter
    logical,    optional, intent(inout) :: conv
    logical,    optional, intent(in)    :: verbose

    logical :: verbose_
    integer :: maxiter
    FLOAT :: tol

    call profiling_in(C_PROFILING_EIGEN_SOLVER)
    call push_sub('eigen.eigen_solver_run')

    verbose_ = .false.; if(present(verbose)) verbose_ = verbose

    if(iter < eigens%final_tol_iter) then
      tol = log(eigens%final_tol/eigens%init_tol)/(eigens%final_tol_iter - 1)*(iter - 1) + &
        log(eigens%init_tol)
      tol = exp(tol)
    else
      tol = eigens%final_tol
    end if

    if(present(conv)) conv = .false.
    maxiter = eigens%es_maxiter

    if (st%d%wfs_type == M_REAL) then 
      select case(eigens%es_type)
      case(RS_CG_NEW)
        call deigen_solver_cg2_new(gr, st, h, tol, maxiter, &
             eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_CG)
        call deigen_solver_cg2(gr, st, h, tol, maxiter, &
             eigens%converged, eigens%diff, verbose = verbose_)
#ifdef HAVE_TRLAN
      case(RS_LANCZOS)
        call eigen_solver_cg3(gr%m, st, h, tol, maxiter, &
             eigens%converged, eigens%diff)
#endif
      case(RS_PLAN)
        call deigen_solver_plan(gr, st, h, tol, maxiter, eigens%converged, eigens%diff)
#if defined(HAVE_ARPACK)
      case(ARPACK)
        call deigen_solver_arpack(gr, st, h, tol, maxiter, eigens%arnoldi_vectors, &
             eigens%converged, eigens%diff)
#endif
      case(EVOLUTION)
        call deigen_solver_evolution(gr, st, h, tol, maxiter, eigens%converged, eigens%diff, &
             tau = eigens%imag_time)
      end select

      if(eigens%do_subspace_diag) call deigen_diagon_subspace(gr, st, h)
    else
      select case(eigens%es_type)
      case(RS_CG_NEW)
        call zeigen_solver_cg2_new(gr, st, h, tol, maxiter, &
             eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_CG)
        call zeigen_solver_cg2(gr, st, h, tol, maxiter, &
             eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_PLAN)
        call zeigen_solver_plan(gr, st, h, tol, maxiter, eigens%converged, eigens%diff)
      case(EVOLUTION)
        call zeigen_solver_evolution(gr, st, h, tol, maxiter, eigens%converged, eigens%diff, &
             tau = eigens%imag_time)
      end select

      if(eigens%do_subspace_diag) call zeigen_diagon_subspace(gr, st, h)
    end if

    eigens%matvec = maxiter
    if(present(conv).and. eigens%converged == st%nst*st%d%nik) conv = .true.

    call pop_sub()
    call profiling_out(C_PROFILING_EIGEN_SOLVER)
  end subroutine eigen_solver_run


#if defined(HAVE_ARPACK)
#include "undef.F90"
#include "real.F90"
#include "eigen_arpack_inc.F90"
#include "undef.F90"
#include "complex.F90"
#include "eigen_arpack_inc.F90"
#endif
#ifdef HAVE_TRLAN
#include "eigen_trlan.F90"
#endif


#include "undef.F90"
#include "real.F90"
#include "eigen_inc.F90"
#include "eigen_cg_inc.F90"
#include "eigen_plan_inc.F90"
#include "eigen_evolution_inc.F90"


#include "undef.F90"
#include "complex.F90"
#include "eigen_inc.F90"
#include "eigen_cg_inc.F90"
#include "eigen_plan_inc.F90"
#include "eigen_evolution_inc.F90"

end module eigen_solver_m
