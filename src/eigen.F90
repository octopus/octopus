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

#include "config_F90.h"

module eigen_solver
use global
use oct_parser
use hamiltonian
use states

implicit none

type eigen_solver_type
  integer :: es_type    ! which eigen solver to use

  ! for new conjugated gradients (es_type = 0)
  real(r8) :: init_tol
  real(r8) :: final_tol
  integer  :: final_tol_iter
  integer  :: es_maxiter  

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

contains

subroutine eigen_solver_init(eigens)
  type(eigen_solver_type), intent(out) :: eigens
  
  sub_name = 'eigen_solver_init'; call push_sub()

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
  case default
    write(message(1), '(a,i4,a)') "Input: '", eigens%es_type, &
         "' is not a valid EigenSolver"
    message(2) = '(0 <= EigenSolver <= 0)'
    call write_fatal(2)
  end select
  call write_info(1)

  call oct_parse_double("EigenSolverInitTolerance", 1.0e-10_r8, eigens%init_tol)
  if(eigens%init_tol <= 0) then
      write(message(1), '(a,e14.4)') "Input: '", eigens%init_tol, &
           "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverInitTolerance >= 0)'
      call write_fatal(2)
  endif

  call oct_parse_double("EigenSolverFinalTolerance", 1.0e-14_r8, eigens%final_tol)
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
  
  call pop_sub(); return
end subroutine eigen_solver_init

subroutine eigen_solver_run(eigens, st, sys, h, iter, conv)
  type(eigen_solver_type), intent(IN) :: eigens
  type(states_type), intent(inout) :: st
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(IN) :: h
  integer, intent(in) :: iter
  logical, intent(inout), optional :: conv

  integer :: ik, maxiter, converged, errorflag
  real(r8) :: tol
  real(r8), allocatable :: diff(:, :)

  sub_name = 'eigen_solver_run'; call push_sub()

  allocate(diff(st%nst, st%nik)); diff = 1.0_r8
  
  if(iter < eigens%final_tol_iter) then
      tol = (eigens%final_tol - eigens%init_tol)/(eigens%final_tol_iter - 1)*(iter - 1) + &
           eigens%init_tol
  else
      tol = eigens%final_tol
  end if

  if(present(conv)) conv = .false.
  maxiter = eigens%es_maxiter
  converged = 0


  select case(eigens%es_type)
  case(RS_CG)
    call eigen_solver_cg2(st, sys, h, &
         tol, maxiter, converged, errorflag, diff)
#ifdef HAVE_TRLAN
  case(RS_LANCZOS)
    call eigen_solver_cg3(st, sys, h, &
         tol, maxiter, converged, errorflag, diff)
#endif
  case(RS_PLAN)
    call eigen_solver_plan(st, sys, h, tol, maxiter, converged, diff)
  end select
  write(message(1),'(a,i5)') 'Info: Converged = ',converged
  write(message(2),'(a,i8)') 'Info: Matrix-Vector multiplications = ', maxiter
  call write_info(2)
  call states_write_eigenvalues(stdout, st%nst, st, diff)

  if(present(conv).and. converged == st%nst*st%nik) conv = .true.

  deallocate(diff)
  call pop_sub(); return
end subroutine eigen_solver_run

#include "eigen_cg2.F90"
#ifdef HAVE_TRLAN
#include "eigen_cg3.F90"
#endif
#include "eigen_plan.F90"

end module eigen_solver
