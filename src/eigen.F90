#include "config.h"

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

contains

subroutine eigen_solver_init(eigens)
  type(eigen_solver_type), intent(out) :: eigens
  
  sub_name = 'eigen_solver_init'; call push_sub()

  eigens%es_type = fdf_integer("EigenSolver", 0)
  if(eigens%es_type < 0 .or. eigens%es_type > 1) then
    write(message(1), '(a,i4,a)') "Input: '", eigens%es_type, &
         "' is not a valid EigenSolver"
    message(2) = '(0 <= EigenSolver <= 1)'
    call write_fatal(2)
  endif
  
  select case (eigens%es_type)
  case(0)
    eigens%init_tol = fdf_double("EigenSolverInitTolerance", 1.0e-10_r8)
    if(eigens%init_tol <= 0) then
      write(message(1), '(a,e14.4)') "Input: '", eigens%init_tol, &
           "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverInitTolerance >= 0)'
      call write_fatal(2)
    endif

    eigens%final_tol = fdf_double("EigenSolverFinalTolerance", 1.0e-14_r8)
    if(eigens%final_tol <= 0 .or. eigens%final_tol > eigens%init_tol) then
      write(message(1),'(a,e14.4)') "Input: '", eigens%init_tol, &
           "' is not a valid EigenSolverInitTolerance"
      message(2) = '(EigenSolverInitTolerance >= 0 and '
      message(3) = ' EigenSolverFinalTolerance < EigenSolverInitTolerance)'
      call write_fatal(3)
    endif

    eigens%final_tol_iter = fdf_integer("EigenSolverFinalToleranceIteration", 20)
    if(eigens%final_tol_iter <= 1) then
      write(message(1),'(a,i5,a)') "Input: '", eigens%final_tol_iter, &
           "' is not a valid EigenSolverFinalToleranceIter"
      message(2) = '(EigenSolverFinalToleranceIter > 1)'
      call write_fatal(2)
    endif

    eigens%es_maxiter = fdf_integer("EigenSolverMaxIter",25)
    if(eigens%es_maxiter < 1) then
      write(message(1),'(a,i5,a)') "Input: '", eigens%es_maxiter, &
           "' is not a valid EigenSolverMaxIter"
      message(3) = '(EigenSolverMaxIter >=1 )'
      call write_fatal(2)
    endif

  case(1)
    eigens%no_cg = fdf_integer("NumberCG", 3);
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

  if(eigens%es_type .ne. 0) then
    if(iter < eigens%final_tol_iter) then
      tol = (eigens%final_tol - eigens%init_tol)/(eigens%final_tol_iter - 1)*(iter - 1) + &
           eigens%init_tol
    else
      tol = eigens%final_tol
    end if
  end if

  select case(eigens%es_type)
  case(0)
    maxiter = eigens%es_maxiter
    converged = 0

    call eigen_solver_cg2(sys, h, sys%st, &
         tol, maxiter, converged, errorflag, diff)
  case(1)
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
