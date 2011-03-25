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
!! $Id: util.F90 3013 2007-06-21 15:53:17Z xavier $

#include "global.h"

! This module is intended to contain simple general-purpose utility functions
! and procedures.

module utils_m
  use global_m
  use messages_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::              &
    get_divisors,        &
    index2axis,          &
    output_tensor,       &
    output_dipole


contains

  ! ---------------------------------------------------------
  subroutine get_divisors(nn, n_divisors, divisors)
    integer, intent(in)    :: nn
    integer, intent(inout) :: n_divisors
    integer, intent(out)   :: divisors(:)

    integer :: ii, max_d

    PUSH_SUB(get_divisors)

    ASSERT(n_divisors > 1)
    max_d = n_divisors

    n_divisors = 1
    divisors(n_divisors) = 1
    do ii = 2, nn / 2
      if(mod(nn, ii)==0) then
        n_divisors = n_divisors + 1

        if(n_divisors > max_d - 1) then
          message(1) = "Internal error in get_divisors. Please increase n_divisors"
          call messages_fatal(1)
        end if

        divisors(n_divisors) = ii
      end if
    end do
    n_divisors = n_divisors + 1
    divisors(n_divisors) = nn

    POP_SUB(get_divisors)
  end subroutine get_divisors


  ! ---------------------------------------------------------
  character function index2axis(idir) result(ch)
    integer, intent(in) :: idir
    
    PUSH_SUB(index2axis)

    select case(idir)
      case(1)
        ch = 'x'
      case(2)
        ch = 'y'
      case(3)
        ch = 'z'
      case(4)
        ch = 'w'
      case default
        write(ch,'(i1)') idir
    end select

    POP_SUB(index2axis)
  end function index2axis


  ! ---------------------------------------------------------
  subroutine output_tensor(iunit, tensor, ndim, unit, write_average)
    integer,           intent(in) :: iunit
    FLOAT,             intent(in) :: tensor(:,:)
    integer,           intent(in) :: ndim
    type(unit_t),      intent(in) :: unit
    logical, optional, intent(in) :: write_average
    
    FLOAT :: trace
    integer :: jj, kk
    logical :: write_average_

    PUSH_SUB(output_tensor)

    write_average_ = .true.
    if(present(write_average)) write_average_ = write_average

    trace = M_z0
    do jj = 1, ndim
      write(iunit, '(3f20.6)') (units_from_atomic(unit, tensor(jj, kk)), kk=1,ndim)
      trace = trace + tensor(jj, jj)
    end do

    trace = units_from_atomic(unit, trace/TOFLOAT(ndim))

    if(write_average_) write(iunit, '(a, f20.6)')  'Isotropic average', trace

    POP_SUB(output_tensor)
  end subroutine output_tensor


  ! ---------------------------------------------------------
  subroutine output_dipole(iunit, dipole, ndim) 
    integer, intent(in) :: iunit
    FLOAT,   intent(in) :: dipole(:)
    integer, intent(in) :: ndim
    
    integer :: idir

    PUSH_SUB(output_dipole)

    write(iunit, '(a,a20,a17)') 'Dipole:', '[' // trim(units_abbrev(units_out%length)) // ']', &
          '[' // trim(units_abbrev(unit_debye)) // ']'
    do idir = 1, ndim
      write(iunit, '(6x,3a,es14.5,3x,2es14.5)') '<', index2axis(idir), '> = ', &
        units_from_atomic(units_out%length, dipole(idir)), units_from_atomic(unit_debye, dipole(idir))
    end do

    POP_SUB(output_dipole)
  end subroutine output_dipole

end module utils_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
