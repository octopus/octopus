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
!! $Id$

! ---------------------------------------------------------
subroutine X(ode_solver_run)(os, func, startval, solutionp, solutionvec)
  type(ode_solver_t), intent(inout) :: os
  R_TYPE,             intent(in)    :: startval(:)
  !> values of the solution only at the endpoint of the interval
  R_TYPE,             intent(out)   :: solutionp(:)
  R_TYPE, optional,   intent(out)   :: solutionvec(:,:) !< full solution for all t

  interface
    subroutine func(size, t, z, res)
      implicit none
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
!> ODE stepping for equally spaced steps
subroutine X(ode_step)(os, func, startval, solutionp, solutionvec)
  type(ode_solver_t), intent(in)  :: os
  R_TYPE,                intent(in)  :: startval(:)
  R_TYPE,                intent(out) :: solutionp(:)
  R_TYPE, optional,      intent(out) :: solutionvec(:,:) !< full solution for all t

  interface
    subroutine func(size, t, z, res)
      implicit none
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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
