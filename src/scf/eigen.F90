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
!! $Id$

#include "global.h"

module eigen_solver_m
  use datasets_m
  use eigen_cg_m
  use eigen_lobpcg_m
  use eigen_rmmdiis_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_parser_m
  use varinfo_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use mpi_lib_m
  use preconditioners_m
  use profiling_m
  use states_m
  use states_block_m
  use subspace_m
  use units_m
  use td_exp_m

  implicit none

  private
  public ::             &
    eigen_solver_t,     &
    eigen_solver_init,  &
    eigen_solver_end,   &
    eigen_solver_run

  type eigen_solver_t
    integer :: es_type    ! which eigensolver to use
    logical :: verbose    ! If true, the solver prints additional information.

    FLOAT   :: init_tol
    FLOAT   :: final_tol
    integer :: final_tol_iter
    integer :: es_maxiter

    integer :: arnoldi_vectors
    FLOAT   :: imag_time

    ! Stores information about how well it performed.
    FLOAT, pointer :: diff(:, :)
    integer        :: matvec
    integer, pointer :: converged(:)

    ! Stores information about the preconditioning.
    type(preconditioner_t) :: pre
  end type eigen_solver_t


  integer, public, parameter :: &
       RS_PLAN    = 11,         &
       RS_CG      =  5,         &
       RS_MG      =  7,         &
       RS_CG_NEW  =  6,         &
       RS_EVO     =  9,         &
       RS_LOBPCG  =  8,         &
       RS_RMMDIIS = 10

contains

  ! ---------------------------------------------------------
  subroutine eigen_solver_init(gr, eigens, st)
    type(grid_t),         intent(inout) :: gr
    type(eigen_solver_t), intent(out)   :: eigens
    type(states_t),       intent(in)    :: st

    call push_sub('eigen.eigen_solver_init')

    !%Variable Eigensolver
    !%Type integer
    !%Default cg
    !%Section SCF::Eigensolver
    !%Description
    !% Decides the eigensolver that obtains the lowest eigenvalues and
    !% eigenfunctions of the Kohn-Sham Hamiltonian. The default is
    !% conjugate gradients (cg), when parallelization in states is
    !% enabled the default is lobpcg.
    !%Option cg 5
    !% Conjugate-gradients algorithm.
    !%Option plan 11
    !% Preconditioned Lanczos scheme.
    !%Option cg_new 6
    !% An alternative conjugate-gradients eigensolver, faster for
    !% larger systems but less mature.
    !%Option evolution 9
    !% Propagation in imaginary time. WARNING: Sometimes it misbehaves. Use with 
    !% caution.
    !%Option lobpcg 8
    !% Locally optimal block-preconditioned conjugate-gradient algorithm
    !% (only available if DevelVersion = yes),
    !% see: A. Knyazev. Toward the Optimal Preconditioned Eigensolver: Locally
    !% Optimal Block Preconditioned Conjugate Gradient Method. SIAM
    !% Journal on Scientific Computing, 23(2):517­541, 2001.
    !%Option rmmdiis 10
    !% Residual minimization scheme, direct inversion in the iterative subspace.
    !%Option multigrid 7
    !% Multigrid eigensolver (experimental).
    !%End
    call loct_parse_int(check_inp('EigenSolver'), RS_CG, eigens%es_type)

    if(st%parallel_in_states .and. .not. eigen_solver_parallel_in_states(eigens)) then
      message(1) = "The selected eigensolver is not parallel in states."
      message(2) = "Please use the lobpcg or rmmdiis eigensolvers."
      call write_fatal(2)
    end if

    !%Variable EigenSolverVerbose
    !%Type logical
    !%Default no
    !%Section SCF::EigenSolver
    !%Description
    !% If enabled the eigensolver prints additional information.
    !%End
    call loct_parse_logical(check_inp('EigenSolverVerbose'), .false., eigens%verbose)
    
    select case(eigens%es_type)
    case(RS_CG_NEW)
    case(RS_MG)
    case(RS_CG)
    case(RS_PLAN)
    case(RS_EVO)
      !%Variable EigenSolverImaginaryTime
      !%Type float
      !%Default 10.0
      !%Section SCF::EigenSolver
      !%Description
      !% The imaginary-time step that is used in the imaginary-time evolution
      !% method to obtain the lowest eigenvalues/eigenvectors.
      !% It must satisfy EigenSolverImaginaryTime > 0.
      !%End
      call loct_parse_float(check_inp('EigenSolverImaginaryTime'), CNST(10.0), eigens%imag_time)
      if(eigens%imag_time <= M_ZERO) call input_error('EigenSolverImaginaryTime')
    case(RS_LOBPCG)
    case(RS_RMMDIIS)
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
    !% Determines the maximum number of iterations that the
    !% eigensolver will perform if the desired tolerance is not
    !% achieved. The default is 25 iterations.
    !%End
    call loct_parse_int(check_inp('EigenSolverMaxIter'), 25, eigens%es_maxiter)
    if(eigens%es_maxiter < 1) call input_error('EigenSolverMaxIter')

    select case(eigens%es_type)
    case(RS_PLAN, RS_CG, RS_LOBPCG, RS_RMMDIIS)
      call preconditioner_init(eigens%pre, gr)
    end select

    nullify(eigens%diff)
    ALLOCATE(eigens%diff(st%nst, st%d%nik), st%nst*st%d%nik)

    ALLOCATE(eigens%converged(st%d%nik), st%d%nik)
    eigens%converged(1:st%d%nik) = 0
    eigens%matvec    = 0

    call pop_sub()

  end subroutine eigen_solver_init


  ! ---------------------------------------------------------
  subroutine eigen_solver_end(eigens)
    type(eigen_solver_t), intent(inout) :: eigens

    select case(eigens%es_type)
    case(RS_PLAN, RS_CG, RS_LOBPCG, RS_RMMDIIS)
      call preconditioner_end(eigens%pre)
    end select

    deallocate(eigens%converged)
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

    verbose_ = eigens%verbose; if(present(verbose)) verbose_ = verbose

    if(iter < eigens%final_tol_iter) then
      tol = log(eigens%final_tol/eigens%init_tol)/(eigens%final_tol_iter - 1)*(iter - 1) + &
        log(eigens%init_tol)
      tol = exp(tol)
    else
      tol = eigens%final_tol
    end if

    if(present(conv)) conv = .false.
    maxiter = eigens%es_maxiter

    if (st%wfs_type == M_REAL) then 
      select case(eigens%es_type)
      case(RS_CG_NEW)
        call deigen_solver_cg2_new(gr, st, h, tol, maxiter, &
             eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_CG)
        call deigen_solver_cg2(gr, st, h, eigens%pre, tol, maxiter, &
             eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_PLAN)
        call deigen_solver_plan(gr, st, h, eigens%pre, tol, maxiter, eigens%converged, eigens%diff)
      case(RS_EVO)
        call deigen_solver_evolution(gr, st, h, tol, maxiter, eigens%converged, eigens%diff, &
             tau = eigens%imag_time)
      case(RS_LOBPCG)
        if(conf%devel_version) then
          call deigen_solver_lobpcg(gr, st, h, eigens%pre, tol, maxiter, eigens%converged, &
            eigens%diff, verbose = verbose_)
        else
          message(1) = 'LOBPCG is still under development. Put'
          message(2) = ''
          message(3) = '  DevelVersion = yes'
          message(4) = ''
          message(5) = 'in your input file if you really want to use it. Be warned.'
          call write_fatal(5)
        end if
      case(RS_MG)
        call deigen_solver_mg(gr, st, h, tol, maxiter, eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_RMMDIIS)
        call deigen_solver_rmmdiis(gr, st, h, eigens%pre, tol, maxiter, eigens%converged, eigens%diff, verbose = verbose_)
      end select

      call dsubspace_diag(gr, st, h, eigens%diff)

    else
      select case(eigens%es_type)
      case(RS_CG_NEW)
        call zeigen_solver_cg2_new(gr, st, h, tol, maxiter, &
             eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_CG)
        call zeigen_solver_cg2(gr, st, h, eigens%pre, tol, maxiter, &
             eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_PLAN)
        call zeigen_solver_plan(gr, st, h, eigens%pre, tol, maxiter, eigens%converged, eigens%diff)
      case(RS_EVO)
        call zeigen_solver_evolution(gr, st, h, tol, maxiter, eigens%converged, eigens%diff, &
             tau = eigens%imag_time)
      case(RS_LOBPCG)
        if(conf%devel_version) then
          call zeigen_solver_lobpcg(gr, st, h, eigens%pre, tol, maxiter, eigens%converged, &
            eigens%diff, verbose = verbose_)
        else
          message(1) = 'LOBPCG is still under development. Put'
          message(2) = ''
          message(3) = '  DevelVersion = yes'
          message(4) = ''
          message(5) = 'in your input file if you really want to use it. Be warned.'
          call write_fatal(5)
        end if
      case(RS_MG)
        call zeigen_solver_mg(gr, st, h, tol, maxiter, eigens%converged, eigens%diff, verbose = verbose_)
      case(RS_RMMDIIS)
        call zeigen_solver_rmmdiis(gr, st, h, eigens%pre, tol, maxiter, eigens%converged, eigens%diff, verbose = verbose_)
      end select

      call zsubspace_diag(gr, st, h)

    end if

    eigens%matvec = maxiter
    if(present(conv).and. sum(eigens%converged(1:st%d%nik)) == st%nst*st%d%nik) conv = .true.

    call pop_sub()
    call profiling_out(C_PROFILING_EIGEN_SOLVER)
  end subroutine eigen_solver_run

  logical function eigen_solver_parallel_in_states(this) result(par_stat)
    type(eigen_solver_t), intent(in) :: this
    
    par_stat = .false.

    select case(this%es_type)
    case(RS_RMMDIIS, RS_LOBPCG)
      par_stat = .true.
    end select
    
  end function eigen_solver_parallel_in_states
    
#include "undef.F90"
#include "real.F90"
#include "eigen_mg_inc.F90"
#include "eigen_plan_inc.F90"
#include "eigen_evolution_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "eigen_mg_inc.F90"
#include "eigen_plan_inc.F90"
#include "eigen_evolution_inc.F90"

end module eigen_solver_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
