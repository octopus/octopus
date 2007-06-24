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

! This module is intended to contain simple general purpose utility functions
! and procedures.

module utils_m
  use global_m
  use messages_m

  implicit none

  private
  public ::                        &
       get_divisors

contains

  ! ---------------------------------------------------------
  subroutine get_divisors(n, n_divisors, divisors)
    integer, intent(in)    :: n
    integer, intent(inout) :: n_divisors
    integer, intent(out)   :: divisors(:)

    integer :: i, max_d

    ASSERT(n_divisors > 1)
    max_d = n_divisors

    n_divisors = 1
    divisors(n_divisors) = 1
    do i = 2, n/2
      if(mod(n, i)==0) then
        n_divisors = n_divisors + 1

        if(n_divisors > max_d - 1) then
          message(1) = "Internal error in get_divisors. Please increase n_divisors"
          call write_fatal(1)
        end if

        divisors(n_divisors) = i        
      end if
    end do
    n_divisors = n_divisors + 1
    divisors(n_divisors) = n
  end subroutine get_divisors

end module utils_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
