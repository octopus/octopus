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

#include "global.h"

module mesh_lib_m
  use global_m

  implicit none

  private
  public :: mesh_index


contains
  ! this function takes care of the boundary conditions
  ! for a given x,y,z it returns the true index of the point
  integer function mesh_index(dim, periodic_dim, nr, Lxyz_inv, ix_) &
    result(index)
    integer, intent(in) :: dim               ! Number of dimensions.
    integer, intent(in) :: periodic_dim      ! Number of periodic dimensions.
    integer, intent(in) :: nr(:, :)          ! Dimensions of the box.
                                             ! (x, y, z) to point-no.
    integer, intent(in) :: &
      Lxyz_inv(            &
      nr(1,1):nr(2,1),     &
      nr(1,2):nr(2,2),     &
      nr(1,3):nr(2,3))
    integer, intent(in) :: ix_(:)            ! Coodinates of requested point.

    integer :: i, ix(3)    ! ix has to go until 3, not sb%dim

    ix = 0
    ix(1:dim) = ix_(1:dim) ! make a local copy that we can change

    index = 1
    do i = 1, dim
      if(ix(i) < nr(1, i)) then    ! first look left
        if(i <= periodic_dim) then ! fold point
          ix(i) = ix(i) + abs(nr(2,i) - nr(1,i) + 1)
        else
          ix(i) = nr(1, i)
          index = 0
        end if
      else if(ix(i) > nr(2, i)) then  ! the same, but on the right
        if(i <= periodic_dim) then
          ix(i) = ix(i) - abs(nr(2,i) - nr(1,i) + 1)
        else
          ix(i) = nr(2, i)
          index = 0
        end if
      end if
    end do

    if(index.ne.0) index = Lxyz_inv(ix(1), ix(2), ix(3))

  end function mesh_index

end module mesh_lib_m
