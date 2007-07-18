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
!! $Id: td_transport.F90 3030 2007-06-25 16:45:05Z marques $

! Management of the mesh points in the interface region.

#include "global.h"

module td_trans_intf_m
  use derivatives_m
  use global_m
  use grid_m
  use math_m
  use mesh_m
  use messages_m

  implicit none

  private
  public ::        &
    intface_t,     &
    intface_init,  &
    intface_end,   &
    intface_index, &
    member_of_intface
  
  integer, parameter, public :: &
    LEFT      = 1,              & ! Lead indices,
    RIGHT     = 2,              & ! L=1, R=2.
    NLEADS    = 2,              & ! Number of leads.
    TRANS_DIR = 1                 ! Transport is in x-direction.

  character(len=5), dimension(NLEADS), parameter, public :: lead_name = &
    (/'left ', 'right'/)

  ! Important: the point numbers of the interface are dense and
  ! sorted in increasing order.
  type intface_t
    integer          :: extent
    integer          :: np
    integer, pointer :: index(:, :)            ! (np_interface, NLEADS)
    integer          :: index_range(2, NLEADS)
  end type intface_t

contains

  ! ---------------------------------------------------------
  ! Calculate the member points of the interface region.
  subroutine intface_init(gr, intface)
    type(grid_t),    intent(in)  :: gr
    type(intface_t), intent(out) :: intface

    logical :: ok

    integer :: il, i
    integer :: from(3, NLEADS)
    integer :: to(3, NLEADS)

    call push_sub('td_trans_intf.intface_init')

    intface%extent = stencil_extent(gr%f_der%der_discr, TRANS_DIR)

    intface%np = intface%extent*gr%m%l(2)*gr%m%l(3)

    ALLOCATE(intface%index(intface%np, NLEADS), intface%np*NLEADS)

    ! Extract submeshes for the interface regions.
    ! 
    ! In 2D:
    !    +-----to(1)-----------------+----from(2)
    !    |       |                   |       |
    !    |       |                   |       |
    ! from(1)----+-----------------to(2)-----+
    from(:, 1) = gr%m%nr(1, :) + gr%m%enlarge
    from(:, 2) = gr%m%nr(2, :) - gr%m%enlarge
    
    to(:, 1) = from(:, 1) + (/intface%extent, gr%m%l(2), gr%m%l(3)/) - 1
    to(:, 2) = from(:, 2) - (/intface%extent, gr%m%l(2), gr%m%l(3)/) + 1

    do il = 1, NLEADS
      call mesh_subset_indices(gr%m, from(:, il), to(:, il), intface%index(:, il))
    end do

    ! Extract begin and end point numbers of interface region.
    ! Sort the array first to make translations between interface point
    ! number and global point number faster.
    do il = 1, NLEADS
      call sort(intface%index(:, il))
      intface%index_range(1, il) = intface%index(1, il)
      intface%index_range(2, il) = intface%index(intface%np, il)
    end do

    ! In debug mode, check if x_intface is dense, i. e. no index
    ! between min(x_index) and max(x_index) is omitted.
    ! As transport is along x-axis, which is the index running slowest,
    ! this should be the case.
    if(in_debug_mode) then
      ok = .true.
      do il = 1, NLEADS
        do i = intface%index_range(1, il), intface%index_range(2, il)
          if(intface%index(i-intface%index_range(1, il)+1, il).ne.i) then
            ok = .false.
          end if
        end do
      end do
      if(.not.ok) then
        message(1) = 'Failed assertion:'
        message(2) = 'td_transport.intface_init: the interface region'
        message(3) = 'has holes.'
        call write_fatal(3)
      end if
    end if

    call pop_sub()
  end subroutine intface_init


  ! ---------------------------------------------------------
  ! Return the interface point number of the global point idx.
  ! For this to work, the intface%index(:, il) array has to
  ! be sorted in increasing order.
  integer function intface_index(idx, intface, il)
    integer,         intent(in) :: idx
    type(intface_t), intent(in) :: intface
    integer,         intent(in) :: il

    call push_sub('td_trans_intf.intface_index')

    intface_index = idx - intface%index_range(1, il) + 1

    call pop_sub()
  end function intface_index
  

  ! ---------------------------------------------------------
  ! Checks if point number idx is an interface point of interface
  ! number il.
  logical function member_of_intface(idx, intface, il)
    integer,         intent(in) :: idx
    type(intface_t), intent(in) :: intface
    integer,         intent(in) :: il

    call push_sub('td_trans_intf.member_of_intface')

    member_of_intface =                       &
      idx.ge.intface%index_range(1, il) .and. &
      idx.le.intface%index_range(2, il)

    call pop_sub()
  end function member_of_intface


  ! ---------------------------------------------------------
  ! Free memory.
  subroutine intface_end(intface)
    type(intface_t), intent(inout) :: intface

    call push_sub('td_trans_intf.intface_end')

    if(associated(intface%index)) then
      deallocate(intface%index)
      nullify(intface%index)
    end if
    
    call pop_sub()
  end subroutine intface_end
end module td_trans_intf_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
