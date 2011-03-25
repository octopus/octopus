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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$


! ---------------------------------------------------------
subroutine X(ode_solver_init)(os)
  type(ode_solver_t), intent(out) :: os

  PUSH_SUB(X(ode_solver_init))

  !%Variable ODESolver
  !%Type integer
  !%Default ode_rk4
  !%Section Math::General
  !%Description
  !% Specifies what kind of ODE solver will be used.
  !%Option ode_rk4 1
  !% Standard Runge-Kutta, 4th order.
  !%Option ode_fb78 2
  !% Fehlberg solver.
  !%Option ode_vr89 3
  !% Verner solver.
  !%Option ode_pd89 4
  !% Prince-Dormand solver.
  !%End
  call parse_integer(datasets_check('ODESolver'),       ODE_RK4, os%solver_type)
  if( os%solver_type .lt. ODE_MINVAL .or. os%solver_type .gt. ODE_MAXVAL ) then
    call input_error(datasets_check('ODESolver'))
  end if

  !%Variable ODESolverNSteps
  !%Type integer
  !%Default 100
  !%Section Math::General
  !%Description
  !% Number of steps which the chosen ODE solver should perform
  !% in the integration interval [a,b] of the ODE.
  !%End
  call parse_integer(datasets_check('ODESolverNSteps'),     100, os%nsteps)

  call X(ode_solver_create)(os)

  POP_SUB(X(ode_solver_init))
end subroutine X(ode_solver_init)


! ---------------------------------------------------------
subroutine X(ode_solver_create)(os)
  type(ode_solver_t), intent(inout) :: os

  PUSH_SUB(X(ode_solver_create))

  select case(os%solver_type)
  case(ODE_RK4)
    os%vsize = 4
    message(1) = 'Info: ode_solver: Using Runge-Kutta, 4th order.'
    call messages_info(1)
  case(ODE_FB78)
    os%vsize = 13
    message(1) = 'Info: ode_solver: Using Fehlberg, 7th/8th order.'
    call messages_info(1)
  case(ODE_VR89)
    os%vsize = 16
    message(1) = 'Info: ode_solver: Using Verner, 8th/9th order.'
    call messages_info(1)
  case(ODE_PD89)
    os%vsize = 13
    message(1) = 'Info: ode_solver: Using Prince-Dormand, 8th/9th order.'
    call messages_info(1)
  end select

  SAFE_ALLOCATE(os%a(1:os%vsize, 1:os%vsize))
  SAFE_ALLOCATE(os%b(1:os%vsize))
  SAFE_ALLOCATE(os%c(1:os%vsize))
  SAFE_ALLOCATE(os%e(1:os%vsize))

  POP_SUB(X(ode_solver_create))
end subroutine X(ode_solver_create)


! ---------------------------------------------------------
subroutine X(ode_solver_run)(os, func, startval, solutionp, solutionvec)
  type(ode_solver_t), intent(inout) :: os
  R_TYPE,             intent(in)    :: startval(:)
  ! values of the solution only at the endpoint of the interval
  R_TYPE,             intent(out)   :: solutionp(:)
  R_TYPE, optional,   intent(out)   :: solutionvec(:,:) ! full solution for all t

  interface
    subroutine func(size, t, z, res)
      integer, intent(in)  :: size
      FLOAT,   intent(in)  :: t
      R_TYPE,  intent(in)  :: z(:)
      R_TYPE,  intent(out) :: res(:)

    end subroutine func
  end interface

  PUSH_SUB(X(ode_solver_run))

  ! initialize array
  solutionp = M_ZERO

  ! setup coefficients
  select case(os%solver_type)
  case(ODE_RK4)
    call ode_rk4_coeff(os)
  case(ODE_FB78)
    call ode_fb78_coeff(os)
  case(ODE_VR89)
    call ode_vr89_coeff(os)
  case(ODE_PD89)
    call ode_pd89_coeff(os)
  case default
    write(message(1), '(a,i4,a)') "Input: '", os%solver_type, &
      "' is not a valid ODE solver"
    message(2) = '( ODE solver =  ode_rk4 | ode_fb7 | ode_vr8 | ode_pd8 )'
    call messages_fatal(2)
  end select

  ! start stepping
  if(present(solutionvec)) then
    call X(ode_step)(os, func, startval, solutionp, solutionvec)
  else
    call X(ode_step)(os, func, startval, solutionp)
  end if

  POP_SUB(X(ode_solver_run))
end subroutine X(ode_solver_run)


! ---------------------------------------------------------
! ODE stepping for equally spaced steps
!
subroutine X(ode_step)(os, func, startval, solutionp, solutionvec)
  type(ode_solver_t), intent(in)  :: os
  R_TYPE,                intent(in)  :: startval(:)
  R_TYPE,                intent(out) :: solutionp(:)
  R_TYPE, optional,      intent(out) :: solutionvec(:,:) ! full solution for all t

  interface
    subroutine func(size, t, z, res)
      integer, intent(in)  :: size
      FLOAT,   intent(in)  :: t
      R_TYPE,  intent(in)  :: z(:)
      R_TYPE,  intent(out) :: res(:)

    end subroutine func
  end interface

  R_TYPE, allocatable :: yn(:), kv(:,:), y0(:)
  FLOAT   :: tn, dh
  integer, parameter :: order = 4
  integer :: ii, jj, kk

  PUSH_SUB(X(ode_step))

  SAFE_ALLOCATE(kv(1:os%nsize, 1:os%vsize))
  SAFE_ALLOCATE(yn(1:os%nsize))
  SAFE_ALLOCATE(y0(1:os%nsize))

  dh = (os%tmax-os%tmin)/real(os%nsteps, REAL_PRECISION)
  tn = os%tmin
  yn = startval

  ! actual ODE stepping
  do ii = 1, os%nsteps
    if (os%full_solution) then
      solutionvec(ii, :) = yn
    end if

    ! calculate auxilary vectors
    do jj = 1, os%vsize
      y0 = yn
      do kk = 1, jj-1
        y0 = y0 + dh*os%a(jj, kk)*kv(:, kk)
      end do
      call func(os%nsize, tn + os%c(jj)*dh, y0, kv(:, jj))
    end do
    ! step forward
    tn = tn + dh
    do jj = 1, os%vsize
      yn = yn + dh*os%b(jj)*kv(:, jj)
    end do

  end do

  solutionp = yn

  SAFE_DEALLOCATE_A(kv)
  SAFE_DEALLOCATE_A(yn)
  SAFE_DEALLOCATE_A(y0)

  POP_SUB(X(ode_step))
end subroutine X(ode_step)


! ---------------------------------------------------------
subroutine X(ode_solver_end)(os)
  type(ode_solver_t), intent(inout) :: os

  PUSH_SUB(X(ode_solver_end))

  ! cleanup
  SAFE_DEALLOCATE_P(os%a)
  SAFE_DEALLOCATE_P(os%b)
  SAFE_DEALLOCATE_P(os%c)
  SAFE_DEALLOCATE_P(os%e)

  POP_SUB(X(ode_solver_end))
end subroutine X(ode_solver_end)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
