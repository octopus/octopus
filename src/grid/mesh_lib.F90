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
!! $Id$

#include "global.h"

module mesh_lib_m

  implicit none

  private
  public :: mesh_index,    &
            index_t
  
  type index_t
    integer, pointer :: Lxyz(:,:)       ! return x, y and z for each point
    integer, pointer :: Lxyz_inv(:,:,:) ! return points # for each xyz
    integer, pointer :: Lxyz_tmp(:,:,:) ! init_1 and init_2
    ! some other vars
    integer :: nr(2, MAX_DIM)              ! dimensions of the box where the points are contained
    integer :: ll(MAX_DIM)                  ! literally n(2,:) - n(1,:) + 1 - 2*enlarge(:)
  end type index_t

contains
  ! this function takes care of the boundary conditions
  ! for a given x,y,z it returns the true index of the point
  integer function mesh_index(dim, nr, Lxyz_inv, ix_) result(index)
    integer, intent(in) :: dim               ! Number of dimensions.
    integer, intent(in) :: nr(:, :)          ! Dimensions of the box.
                                             ! (x, y, z) to point-no.
    integer, intent(in) :: &
      Lxyz_inv(            &
      nr(1,1):nr(2,1),     &
      nr(1,2):nr(2,2),     &
      nr(1,3):nr(2,3))
    integer, intent(in) :: ix_(:)            ! Coodinates of requested point.


    integer :: i, ix(MAX_DIM)      ! ix has to go until 3, not sb%dim

    ix = 0
    ix(1:dim) = ix_(1:dim) ! make a local copy that we can change

    index = 1
    do i = 1, dim
      if(ix(i) < nr(1, i)) then    ! first look left
        ix(i) = nr(1, i)
        index = 0
      else if(ix(i) > nr(2, i)) then  ! the same, but on the right
        ix(i) = nr(2, i)
        index = 0
      end if
    end do
    
    if(index.ne.0) index = Lxyz_inv(ix(1), ix(2), ix(3))

  end function mesh_index

end module mesh_lib_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
