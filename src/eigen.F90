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

#include "global.h"

module eigen_solver
  use global
  use oct_parser
  use linalg
  use io
  use states
  use mesh
  use hamiltonian
  use system

  implicit none

  type eigen_solver_type
    integer :: es_type    ! which eigen solver to use

    ! for new conjugated gradients (es_type = 0)
    FLOAT :: init_tol
    FLOAT :: final_tol
    integer  :: final_tol_iter
    integer  :: es_maxiter  

    ! Stores information about how well it performed.
    FLOAT, pointer :: diff(:, :)
    integer           :: matvec
    integer           :: converged
  end type eigen_solver_type

  integer, parameter :: RS_CG      = 0
#ifdef HAVE_TRLAN
  integer, parameter :: RS_LANCZOS = 1
#endif
  integer, parameter :: RS_PLAN    = 2

  ! For the TRLan package
#ifdef HAVE_TRLAN
  integer                         :: ik_trlan
  type(hamiltonian_type), pointer :: h_trlan
  type(mesh_type), pointer        :: m_trlan
  type(states_type), pointer      :: st_trlan
  type(system_type), pointer      :: sys_trlan
#endif

  type(derivatives_type) :: filter

contains

subroutine eigen_solver_init(eigens, st, m)
  type(eigen_solver_type), intent(out) :: eigens
  type(states_type), intent(in)        :: st
  type(mesh_type), intent(in)          :: m
  
  call push_sub('eigen_solver_init')

  call oct_parse_int("EigenSolver", RS_CG, eigens%es_type)
  select case(eigens%es_type)
  case(RS_CG)
    message(1) = 'Info: Eigensolver type: Real-space conjugate gradients'
#ifdef HAVE_TRLAN
  case(RS_LANCZOS)
    message(1) = 'Info: Eigensolver type: Lanczos algorithm (TRLan package)'
#endif
  case(RS_PLAN)
    message(1) = 'Info: Eigensolver type: Preconditioned Lanczos'
    call derivatives_init_filter(m, 1, filter)
  case default
    write(message(1), '(a,i4,a)') "Input: '", eigens%es_type, &
         "' is not a valid EigenSolver"
    message(2) = '(0 <= EigenSolver <= 0)'
    call write_fatal(2)
  end select
  call write_info(1)

  call oct_parse_float("EigenSolverInitTolerance", CNST(1.0e-10), eigens%init_tol)
  if(eigens%init_tol <= 0) then
      write(message(1), '(a,e14.4)') "Input: '", eigens%init_tol, &
           "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverInitTolerance >= 0)'
      call write_fatal(2)
  endif

  call oct_parse_float("EigenSolverFinalTolerance", CNST(1.0e-14), eigens%final_tol)
  if(eigens%final_tol <= 0 .or. eigens%final_tol > eigens%init_tol) then
      write(message(1),'(a,e14.4)') "Input: '", eigens%init_tol, &
           "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverInitTolerance >= 0 and '
      message(3) = ' EigenSolverFinalTolerance < EigenSolverInitTolerance)'
      call write_fatal(3)
  endif

  call oct_parse_int("EigenSolverFinalToleranceIteration", 7, eigens%final_tol_iter)
  if(eigens%final_tol_iter <= 1) then
      write(message(1),'(a,i5,a)') "Input: '", eigens%final_tol_iter, &
           "' is not a valid EigenSolverFinalToleranceIter"
      message(2) = '(EigenSolverFinalToleranceIter > 1)'
      call write_fatal(2)
  endif

  call oct_parse_int("EigenSolverMaxIter", 25, eigens%es_maxiter)
  if(eigens%es_maxiter < 1) then
      write(message(1),'(a,i5,a)') "Input: '", eigens%es_maxiter, &
           "' is not a valid EigenSolverMaxIter"
      message(3) = '(EigenSolverMaxIter >=1 )'
      call write_fatal(2)
  endif

  nullify(eigens%diff)
  allocate(eigens%diff(st%nst, st%nik))  

  eigens%converged = 0
  eigens%matvec    = 0

  call pop_sub(); return
end subroutine eigen_solver_init

subroutine eigen_solver_end(eigens)
  type(eigen_solver_type), intent(inout) :: eigens
  select case(eigens%es_type)
  case(RS_PLAN)
    call derivatives_end(filter)
  end select
  nullify(eigens%diff)
end subroutine eigen_solver_end

subroutine eigen_solver_run(eigens, st, sys, h, iter, conv)
  type(eigen_solver_type), intent(inout) :: eigens
  type(states_type), intent(inout) :: st
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(IN) :: h
  integer, intent(in) :: iter
  logical, intent(inout), optional :: conv

  integer :: ik, maxiter, errorflag
  FLOAT :: tol

  call push_sub('eigen_solver_run')
  
  if(iter < eigens%final_tol_iter) then
      tol = log(eigens%final_tol/eigens%init_tol)/(eigens%final_tol_iter - 1)*(iter - 1) + &
            log(eigens%init_tol)
      tol = exp(tol)
  else
      tol = eigens%final_tol
  end if

  if(present(conv)) conv = .false.
  maxiter = eigens%es_maxiter
  eigens%converged = 0

  select case(eigens%es_type)
  case(RS_CG)
    call eigen_solver_cg2(st, sys, h, &
         tol, maxiter, eigens%converged, errorflag, eigens%diff)
#ifdef HAVE_TRLAN
  case(RS_LANCZOS)
    call eigen_solver_cg3(st, sys, h, &
         tol, maxiter, eigens%converged, errorflag, eigens%diff)
#endif
  case(RS_PLAN)
    call eigen_solver_plan(st, sys, h, tol, maxiter, eigens%converged, eigens%diff)
  end select

  eigens%matvec = maxiter
  if(present(conv).and. eigens%converged == st%nst*st%nik) conv = .true.

  call pop_sub()
end subroutine eigen_solver_run

!!! This routine in principle diagonalises the hamiltonian in the
!!! basis defined by st. It has not been tested, and it is not used now
subroutine eigen_diagon_subspace(st, sys, h)
  type(states_type), intent(inout)   :: st
  type(system_type), intent(in)      :: sys
  type(hamiltonian_type), intent(in) :: h

  R_TYPE, allocatable :: h_subspace(:,:), vec(:,:), f(:,:,:)
  integer :: ik, i, j

  allocate(h_subspace(st%nst, st%nst), vec(st%nst, st%nst))
  allocate(f(sys%m%np, st%dim, st%nst))

  ik_loop: do ik = 1, st%nik
    f = st%X(psi)(:,:,:, ik)

    eigenfunction_loop : do i = 1, st%nst
      call X(Hpsi)(h, sys%m, st%X(psi)(:,:, i, ik) , f(:,:, 1), sys, ik)
      h_subspace(i, i) = st%eigenval(i, ik)
      do j = i, st%nst
        h_subspace(i, j) = X(states_dotp) (sys%m, st%dim, st%X(psi)(:,:, j, ik), f(:,:, 1))
        h_subspace(j, i) = R_CONJ(h_subspace(i, j))
      end do
    end do eigenfunction_loop

    call X(eigensolve) (st%nst, h_subspace, vec, st%eigenval(:, ik))

    do i = 1, st%nst
      ! build new state
      st%X(psi)(:,:, i, ik) = vec(i, i)*st%X(psi)(:,:, i, ik)
      do j = 1, st%nst
        if(i.ne.j) st%X(psi)(:,:,i, ik) = st%X(psi)(:,:,i, ik) + vec(i, j)*f(:,:,j)
      end do
      
      ! renormalize
      st%X(psi)(:,:, i, ik) = st%X(psi)(:,:, i, ik)/X(states_nrm2)(sys%m, st%dim, st%X(psi)(:,:, i, ik))
    end do
  end do ik_loop

  deallocate(f, h_subspace, vec)

end subroutine eigen_diagon_subspace


#include "eigen_cg2.F90"
#ifdef HAVE_TRLAN
#include "eigen_cg3.F90"
#endif
#include "eigen_plan.F90"

end module eigen_solver
