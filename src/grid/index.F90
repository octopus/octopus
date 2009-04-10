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
#define MIN_DIM 3

module index_m
  use global_m
  use hypercube_m
  use simul_box_m

  implicit none

  private
  public ::                  &
       index_t,              &
       index_from_coords,    &
       index_to_coords
  
  type index_t
    type(hypercube_t)          :: hypercube
    type(simul_box_t), pointer :: sb
    integer                    :: nr(2, MAX_DIM)   ! dimensions of the box where the points are contained
    integer                    :: ll(MAX_DIM)      ! literally nr(2,:) - nr(1,:) + 1 - 2*enlarge(:)
    integer, pointer           :: Lxyz(:,:)        ! return x, y and z for each point
    integer, pointer           :: Lxyz_inv(:,:,:)  ! return points # for each xyz
    integer, pointer           :: Lxyz_tmp(:,:,:)  ! init_1 and init_2
    integer                    :: enlarge(MAX_DIM) ! number of points to add for boundary conditions
  end type index_t

contains

  !
  ! This function takes care of the boundary conditions for a given
  ! vector of integer coordinates it returns the true _global_ index
  ! of the point.
  !
  integer function index_from_coords(idx, dim, ix) result(index)
    type(index_t),      intent(in)    :: idx
    integer,            intent(in)    :: dim
    integer,            intent(in)    :: ix(:)

    integer :: ix2(MAX_DIM), idir

    forall (idir = 1:dim) ix2(idir) = ix(idir)
    forall (idir = dim + 1:MAX_DIM) ix2(idir) = 0

    if(idx%sb%box_shape /= HYPERCUBE) then
      index = idx%Lxyz_inv(ix2(1), ix2(2), ix2(3))
    else
      call hypercube_x_to_i(idx%hypercube, dim, idx%nr, idx%enlarge(1), ix, index)
    end if
    
  end function index_from_coords

  !
  ! Given a _global_ point index, this function returns the set of
  ! integer coordinates of the point.
  !
  pure subroutine index_to_coords(idx, dim, ip, ix)
    type(index_t),      intent(in)    :: idx
    integer,            intent(in)    :: dim
    integer,            intent(in)    :: ip
    integer,            intent(out)   :: ix(:)

    integer :: idir 

    ! We set all ix to zero first (otherwise the non-existent dimesions would be 
    ! undefined on exit).
    ix = 0
    if(idx%sb%box_shape /= HYPERCUBE) then
      forall (idir = 1:dim) ix(idir) = idx%Lxyz(ip, idir)
    else
      call hypercube_i_to_x(idx%hypercube, dim, idx%nr, idx%enlarge(1), ip, ix)
    end if
  end subroutine index_to_coords

  !
  ! This function returns .true. if the coordinates are inside the
  ! inner cube (without the enlargement).
  !
  logical pure function coords_in_inner_cube(idx, dim, ix) result(inside)
    type(index_t),      intent(in)    :: idx
    integer,            intent(in)    :: dim
    integer,            intent(in)    :: ix(:)
    
    integer :: idir

    inside = .true.
    do idir = 1, dim
      inside = inside &
           .and. (ix(idir) > idx%nr(1, idir) + idx%enlarge(idir)) &
           .and. (ix(idir) < idx%nr(2, idir) - idx%enlarge(idir))
    end do

  end function coords_in_inner_cube

end module index_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
