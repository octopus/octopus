!! Copyright (C) 2009 X. Andrade
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
!! $Id: lookup.F90 3613 2007-11-29 16:47:41Z xavier $

#include "global.h"

module lookup_m
  use global_m
  use loct_pointer_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::                &
    lookup_t,              &
    lookup_init,           &
    lookup_end,            &
    lookup_copy,           &
    lookup_get_list
    
  type lookup_t
    private
    integer        :: nobjs
    integer        :: dim
    FLOAT, pointer :: pos(:, :)
  end type lookup_t
  
contains

  subroutine lookup_init(this, dim, nobjs, pos)
    type(lookup_t), intent(out) :: this
    integer,        intent(in)  :: dim
    integer,        intent(in)  :: nobjs
    FLOAT,          intent(in)  :: pos(:, :)
    
    PUSH_SUB(lookup_init)

    this%nobjs = nobjs
    this%dim = dim
    SAFE_ALLOCATE(this%pos(1:this%dim, 1:this%nobjs))
    
    this%pos(1:this%dim, 1:this%nobjs) = pos(1:this%dim, 1:this%nobjs)

    POP_SUB(lookup_init)
  end subroutine lookup_init

  ! -----------------------------------------
  subroutine lookup_end(this)
    type(lookup_t), intent(inout) :: this

    PUSH_SUB(lookup_end)
    SAFE_DEALLOCATE_P(this%pos)

    POP_SUB(lookup_end)
  end subroutine lookup_end

  ! -----------------------------------------

  subroutine lookup_copy(cin, cout)
    type(lookup_t), intent(in) :: cin
    type(lookup_t), intent(out) :: cout

    PUSH_SUB(lookup_copy)

    cout%nobjs = cin%nobjs
    cout%dim = cin%dim
    call loct_pointer_copy(cout%pos, cin%pos)
   
    POP_SUB(lookup_copy)
  end subroutine lookup_copy

  ! ------------------------------------------

  subroutine lookup_get_list(this, npoint, points, radius, nlist, list)
    type(lookup_t), intent(in)  :: this
    integer,        intent(in)  :: npoint
    FLOAT,          intent(in)  :: points(:, :)
    FLOAT,          intent(in)  :: radius
    integer,        intent(out) :: nlist(:)
    integer, optional, pointer  :: list(:, :)

    FLOAT :: r2
    integer :: ii, ipoint

    ! No PUSH SUB, called too often.

    if(present(list)) then
      SAFE_ALLOCATE(list(1:this%nobjs, 1:npoint))
    end if

    nlist(1:npoint) = 0    

    do ipoint = 1, npoint
      do ii = 1, this%nobjs
        r2 = sum((this%pos(1:this%dim, ii) - points(1:this%dim, ipoint))**2)
        if(r2 < radius**2) then
          nlist(ipoint) = nlist(ipoint) + 1
!This is a PGI pragma to force the optimization level of this file to -O0.
!-O2 or below is needed for 10.5. -O1 or below is needed for 10.8.   
!The line after the pragma causes a segmentation fault otherwise.
!pgi$r opt=0
          if(present(list)) list(nlist(ipoint), ipoint) = ii
        end if
      end do
    end do

  end subroutine lookup_get_list

end module lookup_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
