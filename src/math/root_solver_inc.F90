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
subroutine X(root_solver_run)(rs, func, root, success, startval, coeff)
  type(root_solver_t), intent(inout) :: rs
  R_TYPE,                intent(out)  :: root(:)        !< roots we are searching
  logical,               intent(out)  :: success
  R_TYPE, optional,      intent(in)   :: startval(:)    !< start value for the search
  R_TYPE, optional,      intent(in)   :: coeff(:)       !< polynomial coefficients
  interface
    subroutine func(z, f, jf)
      implicit none
      R_TYPE, intent(in)  :: z(:)
      R_TYPE, intent(out) :: f(:), jf(:, :)
    end subroutine func
  end interface

  ! no push_sub, called too often

  ! Initializations
  root = M_ZERO
  success = .false.

  select case(rs%solver_type)
#if defined(R_TREAL)
  case(ROOT_NEWTON)
    call droot_newton(rs, func, root, startval, success)
#endif
  case(ROOT_WATTERSTROM)
#ifdef R_TREAL
    message(1) = 'root_solver: Watterstrom method not defined for pure real arithmetic'
    call messages_fatal(1)
#endif
#ifdef R_TCOMPLEX
    if(present(coeff)) then
      message(1) = 'Info: root_solver: Using Watterstrom method.'
      call messages_info(1)
      call zroot_watterstrom(rs, root, coeff)
    else
      message(1) = 'root_solver: Watterstrom method only valid for polynomials.'
      call messages_fatal(1)
    end if
#endif

  case default
    write(message(1), '(a,i4,a)') "Input: '", rs%solver_type, &
      "' is not a valid root solver"
    message(2) = '( root solver = root_newton | root_watterstrom )'
    call messages_fatal(2)
  end select

end subroutine X(root_solver_run)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
