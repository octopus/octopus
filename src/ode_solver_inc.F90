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


! ---------------------------------------------------------
subroutine X(ode_solver_init)(os)
  type(ode_solver_type), intent(out) :: os

  call push_sub('ode_solver_inc.ode_solver_init')

  !%Variable ODESolver
  !%Type integer
  !%Section 1 Generalities
  !%Description
  !% Specifies what kind of root solver will be used
  !%Option ode_rk4 1
  !% Standard Runge-Kutta 4th order
  !%Option ode_fb78 2
  !% Fehlberg solver
  !%Option ode_vr89 3
  !% Verner solver
  !%Option ode_pd89 4
  !% Prince-Dormand solver
  !%End
  call loct_parse_int(check_inp('ODESolver'),       ODE_RK4, os%solver_type)
  call loct_parse_int(check_inp('ODESolverNSteps'),     100, os%nsteps)
  if (os%solver_type .gt. ODE_PD89) then
     message(1) = 'Error: Unknown ODE solver type'
     call write_fatal(1)
  endif

  call X(ode_solver_create)(os)

  call pop_sub()
end subroutine X(ode_solver_init)


! ---------------------------------------------------------
subroutine X(ode_solver_create)(os)
  type(ode_solver_type), intent(inout) :: os

  call push_sub('ode_solver_inc.ode_solver_create')

  select case(os%solver_type)
  case(ODE_RK4)
     os%vsize = 4
     message(1) = 'Info: ode_solver: Using Runge-Kutta 4th order.'
     call write_info(1)
  case(ODE_FB78)
     os%vsize = 13
     message(1) = 'Info: ode_solver: Using Fehlberg 7th/8th order.'
     call write_info(1)
  case(ODE_VR89)
     os%vsize = 16
     message(1) = 'Info: ode_solver: Using Verner 8th/9th order.'
     call write_info(1)
  case(ODE_PD89)
     os%vsize = 13
     message(1) = 'Info: ode_solver: Using Prince-Dormand 8th/9th order.'
     call write_info(1)
  end select

  allocate(os%a(os%vsize, os%vsize), os%b(os%vsize), os%c(os%vsize), os%e(os%vsize))

  call pop_sub()
end subroutine X(ode_solver_create)


! ---------------------------------------------------------
subroutine X(ode_solver_run)(os, func, startval, solutionp, solutionvec, filename)
  type(ode_solver_type),      intent(inout) :: os
  R_TYPE,                     intent(in)    :: startval(:)
  ! values of the solution only at the endpoint of the interval
  R_TYPE,                     intent(out)   :: solutionp(:)
  R_TYPE, optional,           intent(out)   :: solutionvec(:,:) ! full solution for all t
  character(len=*), optional, intent(in)    :: filename

  character(len=128) :: filename_

  interface
     subroutine func(size, t, z, res)
       integer, intent(in)  :: size
       FLOAT,   intent(in)  :: t
       R_TYPE,  intent(in)  :: z(:)
       R_TYPE,  intent(out) :: res(:)

     end subroutine func
  end interface


  call push_sub('ode_solver_inc.ode_solver_run')

  if (present(filename)) then
     filename_ = trim(filename)
  else
     filename_ = ''
  endif

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
     call write_fatal(2)
  end select

  ! start stepping
  if(present(solutionvec)) then
     call X(ode_step)(os, func, startval, solutionp, solutionvec)
  else
     call X(ode_step)(os, func, startval, solutionp)
  endif


  call pop_sub()
end subroutine X(ode_solver_run)


! ---------------------------------------------------------
! ODE stepping for equally spaced steps
!
subroutine X(ode_step)(os, func, startval, solutionp, solutionvec)
  type(ode_solver_type), intent(in)  :: os
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
  integer :: i, j, k

  call push_sub('ode_solver_inc.ode_step')

  allocate(kv(os%nsize, os%vsize), yn(os%nsize), y0(os%nsize))

  dh = (os%tmax-os%tmin)/real(os%nsteps, PRECISION)
  tn = os%tmin
  yn = startval

  ! actual ODE stepping
  do i = 1, os%nsteps
     if (os%full_solution) then
        solutionvec(i, :) = yn
     endif

     ! calculate auxilary vectors
     do j = 1, os%vsize
        y0 = yn
        do k = 1, j-1
           y0 = y0 + dh*os%a(j,k)*kv(:,k)
        enddo
        call func(os%nsize, tn + os%c(j)*dh, y0, kv(:,j))
     enddo
     ! step forward
     tn = tn + dh
     do j = 1, os%vsize
        yn = yn + dh*os%b(j)*kv(:,j)
     enddo

  enddo

  solutionp = yn

  deallocate(kv, yn, y0)

  call pop_sub()
end subroutine X(ode_step)


! ---------------------------------------------------------
subroutine X(ode_solver_end)(os)
  type(ode_solver_type), intent(in)  :: os

  call push_sub('ode_solver_inc.ode_solver_end')

  ! cleanup
  deallocate(os%a, os%b, os%c, os%e)

  call pop_sub()
end subroutine X(ode_solver_end)

