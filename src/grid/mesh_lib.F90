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
  public ::                  &
       index_t,              &
       index_from_coords,    &
       index_to_coords
  
  type index_t
    integer          :: nr(2, MAX_DIM)              ! dimensions of the box where the points are contained
    integer          :: ll(MAX_DIM)                 ! literally n(2,:) - n(1,:) + 1 - 2*enlarge(:)
    integer, pointer :: Lxyz(:,:)       ! return x, y and z for each point
    integer, pointer :: Lxyz_inv(:,:,:) ! return points # for each xyz
    integer, pointer :: Lxyz_tmp(:,:,:) ! init_1 and init_2
  end type index_t

contains

  ! this function takes care of the boundary conditions
  ! for a given x,y,z it returns the true index of the point
  integer pure function index_from_coords(idx, dim, ix) result(index)
    type(index_t),      intent(in)    :: idx
    integer,            intent(in)    :: dim
    integer,            intent(in)    :: ix(:)  ! Coodinates of requested point.

    integer :: ix2(MAX_DIM), idir

    forall (idir = 1:dim) ix2(idir) = ix(idir)
    forall (idir = dim + 1:MAX_DIM) ix2(idir) = 0

    ! FIXME 4D
    index = idx%Lxyz_inv(ix2(1), ix2(2), ix2(3))
    
  end function index_from_coords
  
  subroutine index_to_coords(idx, dim, ip, ix)
    type(index_t),      intent(in)    :: idx
    integer,            intent(in)    :: dim
    integer,            intent(in)    :: ip
    integer,            intent(out)   :: ix(:)  ! Coodinates of requested point.

    integer :: idir 

    forall (idir = 1:dim) ix(idir) = idx%Lxyz(ip, idir)

  end subroutine index_to_coords

end module mesh_lib_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
