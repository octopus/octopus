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

  ! for old conjugated gradients (es_type = 1)
  integer :: no_cg    ! # of conjugated gradients steps
end type eigen_solver_type

integer, parameter :: ES_NEW_CG = 0, ES_OLD_CG = 1

contains

subroutine eigen_solver_init(eigens)
  type(eigen_solver_type), intent(out) :: eigens
  
  sub_name = 'eigen_solver_init'; call push_sub()

  call oct_parse_int(C_string("EigenSolver"), 0, eigens%es_type)
  if(eigens%es_type < 0 .or. eigens%es_type > 1) then
    write(message(1), '(a,i4,a)') "Input: '", eigens%es_type, &
         "' is not a valid EigenSolver"
    message(2) = '(0 <= EigenSolver <= 1)'
    call write_fatal(2)
  endif
  
  select case (eigens%es_type)
  case(0)
    call oct_parse_double(C_string("EigenSolverInitTolerance"), 1.0e-10_r8, eigens%init_tol)
    if(eigens%init_tol <= 0) then
      write(message(1), '(a,e14.4)') "Input: '", eigens%init_tol, &
           "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverInitTolerance >= 0)'
      call write_fatal(2)
    endif

    call oct_parse_double(C_string("EigenSolverFinalTolerance"), 1.0e-14_r8, eigens%final_tol)
    if(eigens%final_tol <= 0 .or. eigens%final_tol > eigens%init_tol) then
      write(message(1),'(a,e14.4)') "Input: '", eigens%init_tol, &
           "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverInitTolerance >= 0 and '
      message(3) = ' EigenSolverFinalTolerance < EigenSolverInitTolerance)'
      call write_fatal(3)
    endif

    call oct_parse_int(C_string("EigenSolverFinalToleranceIteration"), 7, eigens%final_tol_iter)
    if(eigens%final_tol_iter <= 1) then
      write(message(1),'(a,i5,a)') "Input: '", eigens%final_tol_iter, &
           "' is not a valid EigenSolverFinalToleranceIter"
      message(2) = '(EigenSolverFinalToleranceIter > 1)'
      call write_fatal(2)
    endif

    call oct_parse_int(C_string("EigenSolverMaxIter"), 25, eigens%es_maxiter)
    if(eigens%es_maxiter < 1) then
      write(message(1),'(a,i5,a)') "Input: '", eigens%es_maxiter, &
           "' is not a valid EigenSolverMaxIter"
      message(3) = '(EigenSolverMaxIter >=1 )'
      call write_fatal(2)
    endif

  case(1)
    call oct_parse_int(C_string("NumberCG"), 3, eigens%no_cg)
    if(eigens%no_cg <= 0) then
      write(message(1), '(a,i4,a)') "Input: '", eigens%no_cg, &
           "' is not a valid NumberCG"
      message(2) = '(0 < NumberCG)'
      call write_fatal(2)
    endif

  end select
  
  call pop_sub()
end subroutine eigen_solver_init

subroutine eigen_solver_run(eigens, sys, h, iter, diff)
  type(eigen_solver_type), intent(IN) :: eigens
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(IN) :: h
  integer, intent(in) :: iter
  real(r8), intent(out) :: diff(sys%st%nst, sys%st%nik)

  integer :: ik, maxiter, converged, errorflag
  real(r8) :: tol

  sub_name = 'eigen_solver_run'; call push_sub()

  if(eigens%es_type .ne. ES_OLD_CG) then
    if(iter < eigens%final_tol_iter) then
      tol = (eigens%final_tol - eigens%init_tol)/(eigens%final_tol_iter - 1)*(iter - 1) + &
           eigens%init_tol
    else
      tol = eigens%final_tol
    end if
  end if

  select case(eigens%es_type)
  case(ES_NEW_CG)
    maxiter = eigens%es_maxiter
    converged = 0

    call eigen_solver_cg2(sys, h, sys%st, &
         tol, maxiter, converged, errorflag, diff)
    write(message(1),'(a,i5)') 'Info: Converged = ',converged
    call write_info(1)
  case(ES_OLD_CG)
    call eigen_solver_cg1(eigens%no_cg, sys, h, sys%st, diff)
    do ik = 1, sys%st%nik
      call R_FUNC(states_gram_schmidt) (sys%st%nst, sys%m, sys%st%dim, &
           sys%st%R_FUNC(psi)(:,:,:, ik))
    end do
  end select

  call pop_sub()
end subroutine eigen_solver_run

#include "eigen_cg1.F90"
#include "eigen_cg2.F90"

end module eigen_solver
