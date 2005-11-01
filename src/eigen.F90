!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module eigen_solver
  use global
  use messages
  use syslabels
  use lib_oct_parser
  use lib_basic_alg
  use lib_adv_alg
  use math, only : sort
  use io
  use nl_operator
  use stencil_star
  use mesh_function
  use mesh
  use grid
  use functions
  use states
  use hamiltonian

  implicit none

  private
  public ::             &
    eigen_solver_type,  &
    eigen_solver_init,  &
    eigen_solver_end,   &
    eigen_solver_run

  type eigen_solver_type
    integer :: es_type    ! which eigen solver to use

    FLOAT :: init_tol
    FLOAT :: final_tol
    integer  :: final_tol_iter
    integer  :: es_maxiter

    integer :: arnoldi_vectors
    FLOAT :: imag_time

    ! Stores information about how well it performed.
    FLOAT, pointer :: diff(:, :)
    integer           :: matvec
    integer           :: converged
  end type eigen_solver_type


#ifdef HAVE_TRLAN
  integer, parameter :: RS_LANCZOS = 1
#endif
  integer, parameter :: RS_PLAN    = 2
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
  integer                         :: ik_trlan
  type(hamiltonian_type), pointer :: h_trlan
  type(mesh_type),        pointer :: m_trlan
  type(states_type),      pointer :: st_trlan
#endif

  type(hamiltonian_type), pointer :: h_
  type(grid_type),        pointer :: gr_
  integer                         :: ik_

  type(nl_operator_type) :: filter

contains

  ! ---------------------------------------------------------
  subroutine eigen_solver_init(gr, eigens, st, max_iter_default)
    type(grid_type),         intent(in)  :: gr
    type(eigen_solver_type), intent(out) :: eigens
    type(states_type),       intent(in)  :: st
    integer,                 intent(in)  :: max_iter_default

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
    !% Lanczos scheme, as implemented in the TRLan package.
    !%Option plan 2
    !% Preconditioned Lanczos scheme.
    !%Option arpack 3
    !% Implicitly Restarted Arnoldi Method, as implemented in the ARPACK package
    !%Option jdqz 5
    !% Jacobi-Davidson scheme as implemented in the JDQZ package
    !%Option cg_new 6
    !% A rewritting of the cg option, that will eventually substitute it.
    !%Option evolution 7
    !% Propagation in imaginary time
    !%End
    call loct_parse_int(check_inp('EigenSolver'), RS_CG, eigens%es_type)
    select case(eigens%es_type)
    case(RS_CG_NEW)
      message(1) = 'Info: Eigensolver type: Real-space conjugate gradients (new)'
    case(RS_CG)
      message(1) = 'Info: Eigensolver type: Real-space conjugate gradients'
#ifdef HAVE_TRLAN
    case(RS_LANCZOS)
      message(1) = 'Info: Eigensolver type: Lanczos algorithm (TRLan package)'
#endif
    case(RS_PLAN)
      message(1) = 'Info: Eigensolver type: Preconditioned Lanczos'
#if defined(HAVE_ARPACK)
    case(ARPACK)
      message(1) = 'Info: Eigensolver type: ARPACK Arnoldi method'

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
      message(1) = 'Info: Eigensolver type: Jacobi-Davidson through the JDQZ package'
#endif
    case(EVOLUTION)
      message(1) = 'Info: Eigensolver type: Imaginary time evolution'

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
      write(message(1), '(a,i4,a)') "Input: '", eigens%es_type, &
        "' is not a valid EigenSolver"
      message(2) = '( EigenSolver =  cg | trlan | plan | arpack | cg_new )'
      call write_fatal(2)
    end select
    call write_info(1)

    call loct_parse_float(check_inp('EigenSolverInitTolerance'), CNST(1.0e-6), eigens%init_tol)
    if(eigens%init_tol < 0) then
      write(message(1), '(a,e14.4,a)') "Input: '", eigens%init_tol, &
        "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverInitTolerance > 0)'
      call write_fatal(2)
    end if

    call loct_parse_float(check_inp('EigenSolverFinalTolerance'), CNST(1.0e-6), eigens%final_tol)
    if(eigens%final_tol < 0 .or. eigens%final_tol > eigens%init_tol) then
      write(message(1),'(a,e14.4,a)') "Input: '", eigens%init_tol, &
        "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverFinalTolerance > 0 and '
      message(3) = ' EigenSolverFinalTolerance < EigenSolverInitTolerance)'
      call write_fatal(3)
    end if

    call loct_parse_int(check_inp('EigenSolverFinalToleranceIteration'), 7, eigens%final_tol_iter)
    if(eigens%final_tol_iter <= 1) then
      write(message(1),'(a,i5,a)') "Input: '", eigens%final_tol_iter, &
        "' is not a valid EigenSolverFinalToleranceIter"
      message(2) = '(EigenSolverFinalToleranceIter > 1)'
      call write_fatal(2)
    end if

    call loct_parse_int(check_inp('EigenSolverMaxIter'), max_iter_default, eigens%es_maxiter)
    if(eigens%es_maxiter < 1) then
      write(message(1),'(a,i5,a)') "Input: '", eigens%es_maxiter, &
        "' is not a valid EigenSolverMaxIter"
      message(3) = '(EigenSolverMaxIter >=1 )'
      call write_fatal(2)
    end if

    select case(eigens%es_type)
    case(RS_PLAN)
      call init_filter()
    end select

    nullify(eigens%diff)
    allocate(eigens%diff(st%nst, st%d%nik))

    eigens%converged = 0
    eigens%matvec    = 0

    call pop_sub()

  contains
    subroutine init_filter()
      FLOAT, parameter :: alpha = M_HALF

      ! the filter has a star stencil like the laplacian
      call nl_operator_init(filter, 2*NDIM + 1)
      call stencil_star_get_lapl(NDIM, 1, filter%stencil)
      call nl_operator_build(gr%m, filter, NP, const_w = .true.)

      filter%w_re(1, 1) = alpha
      filter%w_re(2:,1) = M_HALF*(M_ONE-alpha)/NDIM

    end subroutine init_filter
  end subroutine eigen_solver_init


  ! ---------------------------------------------------------
  subroutine eigen_solver_end(eigens)
    type(eigen_solver_type), intent(inout) :: eigens

    select case(eigens%es_type)
    case(RS_PLAN)
      call nl_operator_end(filter)
    end select

    deallocate(eigens%diff); nullify(eigens%diff)

  end subroutine eigen_solver_end


  ! ---------------------------------------------------------
  subroutine eigen_solver_run(eigens, gr, st, h, iter, conv, verbose)
    type(eigen_solver_type), intent(inout) :: eigens
    type(grid_type),         intent(inout) :: gr
    type(states_type),       intent(inout) :: st
    type(hamiltonian_type),  intent(IN)    :: h
    integer,                 intent(in)    :: iter
    logical,       optional, intent(inout) :: conv
    logical,       optional, intent(in)    :: verbose

    logical :: verbose_
    integer :: maxiter
    FLOAT :: tol

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

    select case(eigens%es_type)
    case(RS_CG_NEW)
      call eigen_solver_cg2_new(gr, st, h, tol, maxiter, &
        eigens%converged, eigens%diff, verbose = verbose_)
    case(RS_CG)
      call eigen_solver_cg2(gr, st, h, tol, maxiter, &
        eigens%converged, eigens%diff, verbose = verbose_)
#ifdef HAVE_TRLAN
    case(RS_LANCZOS)
      call eigen_solver_cg3(gr%m, st, h, tol, maxiter, &
        eigens%converged, eigens%diff)
#endif
    case(RS_PLAN)
      call eigen_solver_plan(gr, st, h, tol, maxiter, eigens%converged, eigens%diff)
#if defined(HAVE_ARPACK)
    case(ARPACK)
      call eigen_solver_arpack(gr, st, h, tol, maxiter, eigens%arnoldi_vectors, &
        eigens%converged, eigens%diff)
#endif
    case(EVOLUTION)
      call eigen_solver_evolution(gr, st, h, tol, maxiter, eigens%converged, eigens%diff, &
        tau = eigens%imag_time)
    end select

    eigens%matvec = maxiter
    if(present(conv).and. eigens%converged == st%nst*st%d%nik) conv = .true.

    call pop_sub()
  end subroutine eigen_solver_run


!!! This routine in principle diagonalises the hamiltonian in the
!!! basis defined by st. It has not been tested, and it is not used now
  subroutine eigen_diagon_subspace(gr, st, h)
    type(grid_type),        intent(inout) :: gr
    type(states_type),      intent(inout) :: st
    type(hamiltonian_type), intent(IN)    :: h

    R_TYPE, allocatable :: h_subspace(:,:), vec(:,:), f(:,:,:)
    integer :: ik, i, j

    allocate(h_subspace(st%nst, st%nst), vec(st%nst, st%nst))
    allocate(f(NP, st%d%dim, st%nst))

    ik_loop: do ik = 1, st%d%nik
      f = st%X(psi)(:,:,:, ik)

      eigenfunction_loop : do i = 1, st%nst
        call X(Hpsi)(h, gr, st%X(psi)(:,:, i, ik) , f(:,:, 1), ik)
        h_subspace(i, i) = st%eigenval(i, ik)
        do j = i, st%nst
          h_subspace(i, j) = X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:,:, j, ik), f(:,:, 1))
          h_subspace(j, i) = R_CONJ(h_subspace(i, j))
        end do
      end do eigenfunction_loop

      call lalg_eigensolve(st%nst, h_subspace, vec, st%eigenval(:, ik))

      do i = 1, st%nst
        ! build new state
        st%X(psi)(:,:, i, ik) = vec(i, i)*st%X(psi)(:,:, i, ik)
        do j = 1, st%nst
          if(i.ne.j) st%X(psi)(:,:,i, ik) = st%X(psi)(:,:,i, ik) + vec(i, j)*f(:,:,j)
        end do

        ! renormalize
        st%X(psi)(:,:, i, ik) = st%X(psi)(:,:, i, ik)/X(states_nrm2)(gr%m, st%d%dim, st%X(psi)(:,:, i, ik))
      end do
    end do ik_loop

    deallocate(f, h_subspace, vec)

  end subroutine eigen_diagon_subspace


#include "eigen_cg.F90"
#ifdef HAVE_TRLAN
#include "eigen_trlan.F90"
#endif
#include "eigen_plan.F90"
#if defined(HAVE_ARPACK)
#include "eigen_arpack.F90"
#endif
#include "eigen_evolution.F90"

end module eigen_solver
