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
!! $Id: subarray_inc.F90 3030 2007-06-25 16:45:05Z marques $

subroutine X(subarray_gather)(this, array, subarray)
    type(subarray_t),    intent(in)  :: this
    R_TYPE,              intent(in)  :: array(:)
    R_TYPE,              intent(out) :: subarray(:)

    type(profile_t), save :: prof
    integer :: iblock, ii, isa
    
    call profiling_in(prof, "SUBARRAY_GATHER")

    isa = 0
    do iblock = 1, this%nblocks
      do ii = 1, this%blength(iblock)
        subarray(isa + ii) = array(this%offsets(iblock) + ii - 1)
      end do
      isa = isa + this%blength(iblock)
    end do

    call profiling_out(prof)
end subroutine X(subarray_gather)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
