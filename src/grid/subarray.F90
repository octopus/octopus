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
!! $Id: subarray.F90 4396 2008-07-21 16:14:17Z xavier $

#include "global.h"

module subarray_m
  use batch_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use global_m
  use math_m
  use messages_m
  use opencl_m
  use types_m
  use profiling_m


  implicit none

  private

  public                         &
       subarray_t,               &
       subarray_init,            &
       subarray_end,             &
       subarray_size,            &
       isubarray_gather,         &
       dsubarray_gather,         &
       zsubarray_gather,         &
       dsubarray_gather_batch,   &
       zsubarray_gather_batch,   &
       get_blocks
  
  type subarray_t
    private
    integer, pointer :: offsets(:)
    integer, pointer :: blength(:)
    integer          :: nblocks
    integer          :: npoints
  end type subarray_t
  
contains

  subroutine get_blocks(npos, positions, nblocks, blocklengths, offsets)
    integer, intent(in)  :: npos
    integer, intent(in)  :: positions(:)
    integer, intent(out) :: nblocks
    integer, intent(out) :: blocklengths(:)
    integer, intent(out) :: offsets(:)

    integer :: ip, ib, ii

    blocklengths = 1

    nblocks = 1
    offsets(1) = positions(1)
    blocklengths(1) = 1

    do ip = 2, npos
      if(positions(ip) == positions(ip - 1) + 1) then
        blocklengths(nblocks) = blocklengths(nblocks) + 1
      else
        nblocks = nblocks + 1
        offsets(nblocks) = positions(ip)
      end if
    end do

    ! All the rest is only for checking that the function works
    ! correctly.

    ASSERT( npos == sum(blocklengths(1:nblocks)) )

    ip = 1 
    do ib = 1, nblocks
      do ii = 1, blocklengths(ib)
        ASSERT(positions(ip) == offsets(ib) + ii - 1)
        ip = ip + 1
      end do
    end do
    
    ASSERT(npos == ip - 1)

  end subroutine get_blocks

  subroutine subarray_init(this, total, positions)
    type(subarray_t),    intent(out) :: this
    integer,             intent(in)  :: total
    integer,             intent(in)  :: positions(:)

    integer, allocatable :: bltmp(:), ostmp(:)

    PUSH_SUB(subarray_init)

    this%npoints = total
    
    SAFE_ALLOCATE(bltmp(1:total))
    SAFE_ALLOCATE(ostmp(1:total))

    call get_blocks(total, positions, this%nblocks, bltmp, ostmp)
    
    SAFE_ALLOCATE(this%blength(1:this%nblocks))
    SAFE_ALLOCATE(this%offsets(1:this%nblocks))

    this%blength(1:this%nblocks) = bltmp(1:this%nblocks)
    this%offsets(1:this%nblocks) = ostmp(1:this%nblocks)

    SAFE_DEALLOCATE_A(bltmp)
    SAFE_DEALLOCATE_A(ostmp)

    POP_SUB(subarray_init)
  end subroutine subarray_init

  subroutine subarray_end(this)
    type(subarray_t), intent(inout) :: this
    
    PUSH_SUB(subarray_end)

    SAFE_DEALLOCATE_P(this%offsets)
    SAFE_DEALLOCATE_P(this%blength)

    POP_SUB(subarray_end)
  end subroutine subarray_end

  integer pure function subarray_size(this) result(size)
    type(subarray_t), intent(in) :: this

    size = this%npoints
  end function subarray_size

#include "undef.F90"
#include "integer.F90"
#include "subarray_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "subarray_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "subarray_inc.F90"

end module subarray_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
