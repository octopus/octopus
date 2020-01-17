!! Copyright (C) 2009 N. Helbig, X. Andrade
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module hypercube_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                           &
    hypercube_t,                   &
    hypercube_init,                &
    hypercube_nullify,             &
    hypercube_end,                 &
    hypercube_i_to_x,              &
    hypercube_x_to_i,              &
    hypercube_number_inner_points, &
    hypercube_number_total_points

  type hypercube_t
    private
    integer, pointer :: boxdim(:)
  end type hypercube_t

contains

  subroutine hypercube_init(this, ndim, nr, enlarge)
    type(hypercube_t), intent(out) :: this
    integer,           intent(in)  :: ndim
    integer,           intent(in)  :: nr(:, :)
    integer,           intent(in)  :: enlarge

    integer, allocatable :: npoints(:)
    integer :: ii, jj

    PUSH_SUB(hypercube_init)

    SAFE_ALLOCATE(this%boxdim(1:ndim + 1))
    SAFE_ALLOCATE(npoints(1:ndim))

    forall (ii = 1:ndim) npoints(ii) = nr(2,ii) - nr(1,ii) + 1

    !number of points in each box
    this%boxdim = 1

    !first box is inner points
    do jj = 1, ndim
      this%boxdim(1) = this%boxdim(1)*(npoints(jj)-2*enlarge)
    end do

    !all other boxes are boundary points

    do ii=2, ndim+1
      do jj=1, ii-2
        this%boxdim(ii)=this%boxdim(ii)*(npoints(jj)-2*enlarge)
      end do
      this%boxdim(ii)=this%boxdim(ii)*2*enlarge
      do jj=ii, ndim
        this%boxdim(ii)=this%boxdim(ii)*npoints(jj)
      end do
    end do

    SAFE_DEALLOCATE_A(npoints)

    POP_SUB(hypercube_init)
  end subroutine hypercube_init


  ! ---------------------------------------------------------
  subroutine hypercube_nullify(this)
    type(hypercube_t), intent(inout) :: this

    PUSH_SUB(hypercube_nullify)

    nullify(this%boxdim)

    POP_SUB(hypercube_nullify)
  end subroutine hypercube_nullify


  ! ---------------------------------------------------------
  subroutine hypercube_end(this)
    type(hypercube_t), intent(inout) :: this

    PUSH_SUB(hypercube_end)

    SAFE_DEALLOCATE_P(this%boxdim)

    POP_SUB(hypercube_end)
  end subroutine hypercube_end


  ! ---------------------------------------------------------
  subroutine hypercube_x_to_i(this, ndim, nr, enlarge, coord, icoord)
    type(hypercube_t), intent(in)  :: this
    integer,           intent(in)  :: ndim
    integer,           intent(in)  :: nr(:,:) !< (1:2,1:ndim)
    integer,           intent(in)  :: enlarge
    integer,           intent(in)  :: coord(:) !< (1:ndim)
    integer,           intent(out) :: icoord

    integer :: boxnumb
    integer :: border(1:MAX_DIM), npoints(1:MAX_DIM), lowerb(1:MAX_DIM), tempcoord(1:MAX_DIM)
    integer :: ii, jj

    ! no push sub, called too frequently

    do ii = 1, ndim
      npoints(ii) = nr(2, ii) - nr(1, ii) + 1
      border(ii) = nr(1, ii) + 2*enlarge
      lowerb(ii) = nr(1, ii)
    end do

    !move coordinates such that inner box is in the upper right corner
    do ii = 1, ndim
      tempcoord(ii) = coord(ii)
      tempcoord(ii) = tempcoord(ii) + enlarge - nr(1,ii)
      tempcoord(ii) = mod(tempcoord(ii), npoints(ii))
      tempcoord(ii) = tempcoord(ii) + nr(1,ii)
    end do

    !determine which box we are in
    boxnumb = 1
    do ii = 1, ndim
      if(tempcoord(ii) < border(ii)) then
        boxnumb = ii + 1
        exit
      end if
    end do

    !transform coordinates
    icoord = 0

    if(boxnumb == 1) then
      npoints = npoints - 2*enlarge
      do ii = ndim, 1, -1
        icoord = icoord*npoints(ii)
        icoord = icoord + tempcoord(ii) - border(ii)
      end do
      icoord = icoord+1
      if(icoord > this%boxdim(1) .or. icoord < 1) then
        message(1) = "hypercube box point outside box"
        call messages_fatal(1)
      end if
    else
      do jj = 1, boxnumb - 2
        npoints(jj) = npoints(jj) - 2*enlarge
        lowerb(jj) = nr(1, jj) + 2*enlarge
      end do
      npoints(boxnumb-1) = 2*enlarge
      do jj=ndim, 1, -1
        icoord = icoord*npoints(jj)
        icoord = icoord + (tempcoord(jj) - lowerb(jj))
      end do
      icoord = icoord + 1
      if(icoord > this%boxdim(boxnumb) .or. icoord < 1) then
        message(1) = "hypercube box point outside box"
        call messages_fatal(1)
      else
        do jj = 1, boxnumb - 1
          icoord = icoord + this%boxdim(jj)
        end do
      end if
    end if

  end subroutine hypercube_x_to_i


  ! ---------------------------------------------------------
  pure subroutine hypercube_i_to_x(this, ndim, nr, enlarge, icoord, coord)
    type(hypercube_t), intent(in)  :: this
    integer,           intent(in)  :: ndim
    integer,           intent(in)  :: nr(:,:) !< (1:2,1:ndim)
    integer,           intent(in)  :: enlarge
    integer,           intent(in)  :: icoord
    integer,           intent(out) :: coord(:) !< (1:ndim)

    integer :: boxnumb, tempcoord
    integer :: border(1:MAX_DIM), npoints(1:MAX_DIM), lowerb(1:MAX_DIM)
    integer :: ii, jj

    boxnumb = 1

    jj = 0

    do ii = 1, ndim
      jj = jj + this%boxdim(ii)
      if(icoord > jj) then
        boxnumb = ii + 1
      end if
    end do

    do ii = 1, ndim
      npoints(ii) = nr(2, ii) - nr(1, ii)+1
      border(ii) = nr(1, ii) + 2*enlarge
      lowerb(ii) = nr(1, ii)
    end do

    tempcoord=icoord

    if(boxnumb == 1) then
      npoints = npoints - 2*enlarge
      tempcoord = tempcoord - 1
      do ii = 1, ndim
        coord(ii) = mod(tempcoord, npoints(ii))
        tempcoord = tempcoord - coord(ii)
        tempcoord = tempcoord/npoints(ii)
        coord(ii) = coord(ii) + border(ii)
      end do
    else
      do ii = 1, boxnumb - 2
        npoints(ii) = npoints(ii) - 2*enlarge
        lowerb(ii) = nr(1,ii) + 2*enlarge
        tempcoord = tempcoord - this%boxdim(ii)
      end do
      npoints(boxnumb - 1) = 2*enlarge
      tempcoord = tempcoord - this%boxdim(boxnumb-1)
      tempcoord = tempcoord - 1
      do ii = 1, ndim
        coord(ii) = mod(tempcoord, npoints(ii))
        tempcoord = tempcoord - coord(ii)
        tempcoord = tempcoord/npoints(ii)
        coord(ii) = coord(ii) + lowerb(ii)
      end do
    end if

    do ii = 1, ndim
      npoints(ii) = nr(2,ii) - nr(1,ii) + 1
    end do

    !move inner box back to the middle

    do ii = 1, ndim
      coord(ii) = coord(ii) - nr(1,ii) - enlarge
      if(coord(ii) < 0) then
        coord(ii) = coord(ii) + npoints(ii)
      end if
      coord(ii) = coord(ii) + nr(1,ii)
    end do

  end subroutine hypercube_i_to_x


  ! ---------------------------------------------------------
  pure integer function hypercube_number_inner_points(this) result(number)
    type(hypercube_t), intent(in)  :: this

    number = this%boxdim(1)
  end function hypercube_number_inner_points


  ! ---------------------------------------------------------
  pure integer function hypercube_number_total_points(this) result(number)
    type(hypercube_t), intent(in)  :: this

    number = sum(this%boxdim)
  end function hypercube_number_total_points

end module hypercube_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
