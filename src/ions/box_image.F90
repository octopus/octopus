!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 M. Oliveira
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

module box_image_oct_m
  use box_shape_oct_m
  use iso_c_binding
  use gdlib_oct_m
  use global_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  public :: box_image_t

  !> Class implementing a box generated from a 2D image.
  type, extends(box_shape_t) :: box_image_t
    private
    integer                       :: image_size(2) !< size of the image in each direction in pixels
    FLOAT, public                 :: pixel_size(2) !< size of a pixel in atomic units
    type(c_ptr)                   :: image         !< libgd handler
    character(len=:), allocatable :: filename      !< name of image file
  contains
    procedure :: contains_points => box_image_contains_points
    procedure :: write_info => box_image_write_info
    procedure :: write_short_info => box_image_write_short_info
    final     :: box_image_finalize
  end type box_image_t

  interface box_image_t
    procedure box_image_constructor
  end interface box_image_t

contains

  !--------------------------------------------------------------
  function box_image_constructor(center, lsize, filename, periodic_dim, namespace) result(box)
    FLOAT,             intent(in)    :: center(2)
    FLOAT,             intent(inout) :: lsize(2) !< Length of the image along each Cartesian
                                                 !! direction. Since it might be modified here, we currently
                                                 !! need to return the new value, as this is still needed in
                                                 !! simul_box_t. Later this should be changed.
    character(len=*),  intent(in)    :: filename
    integer,           intent(in)    :: periodic_dim
    type(namespace_t), intent(in)    :: namespace
    class(box_image_t), pointer :: box

    logical :: found
    integer :: idir, box_npts

    PUSH_SUB(box_image_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, 2, center)

    box%filename = trim(filename)
    inquire(file=trim(box%filename), exist=found)
    if (.not. found) then
      message(1) = "Could not find file '" // trim(box%filename) // "' for BoxShape = box_image."

      deallocate(box%filename)
      box%filename = trim(conf%share) // '/' // trim(filename)
      inquire(file=trim(box%filename), exist=found)

      if (.not. found) call messages_fatal(1, namespace=namespace)
    end if

#ifdef HAVE_GDLIB
    box%image = gdlib_image_create_from(box%filename)
#endif
    if (.not. c_associated(box%image)) then
      message(1) = "Could not open file '" // trim(box%filename) // "' for BoxShape = box_image."
      call messages_fatal(1, namespace=namespace)
    end if
#ifdef HAVE_GDLIB
    box%image_size(1) = gdlib_image_sx(box%image)
    box%image_size(2) = gdlib_image_sy(box%image)
#endif

    ! If necessary, adjust lsize to ensure that we always have a pixel at the
    ! origin and that the edges always fall between two pixels.
    do idir = 1, 2
      box_npts = box%image_size(idir)
      if((idir >  periodic_dim .and. even(box%image_size(idir))) .or. &
        (idir <= periodic_dim .and.  odd(box%image_size(idir)))) then
        box_npts = box_npts + 1
        lsize(idir) = lsize(idir) * box_npts / box%image_size(idir)
      end if
    end do

    ! Calculate the size of a pixel. To have one grid point = one pixel the
    ! spacing must be the same as the pixel size.
    box%pixel_size = (M_TWO*lsize)/box%image_size

    POP_SUB(box_image_constructor)
  end function box_image_constructor

  !--------------------------------------------------------------
  subroutine box_image_finalize(this)
    type(box_image_t), intent(inout) :: this

    PUSH_SUB(box_image_finalize)

    call box_shape_end(this)
    if (allocated(this%filename)) then
      deallocate(this%filename)
    end if
#ifdef HAVE_GDLIB
    call gdlib_imagedestroy(this%image)
#endif

    POP_SUB(box_image_finalize)
  end subroutine box_image_finalize

  !--------------------------------------------------------------
  function box_image_contains_points(this, nn, xx) result(contained)
    class(box_image_t), intent(in)  :: this
    integer,            intent(in)  :: nn
    FLOAT,              intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip
    integer :: red, green, blue, ix, iy

    do ip = 1, nn
      ! Transform our cartesian coordinates into pixel coordinates.
      ! Why the minus sign for y? Explanation: http://biolinx.bios.niu.edu/bios546/gd_mod.htm
      ! For reasons that probably made sense to someone at some time, computer graphic coordinates are not the same
      ! as in standard graphing. ... The top left corner of the screen is (0,0).
      ix =   nint((xx(ip, 1) - this%center(1))/this%pixel_size(1)) + (this%image_size(1) - 1)/2
      iy = - nint((xx(ip, 2) - this%center(2))/this%pixel_size(2)) + (this%image_size(2) - 1)/2

#if defined(HAVE_GDLIB)
      call gdlib_image_get_pixel_rgb(this%image, ix, iy, red, green, blue)
#endif
      contained(ip) = (red == 255 .and. green == 255 .and. blue == 255) .neqv. this%is_inside_out()
    end do

  end function box_image_contains_points

  !--------------------------------------------------------------
  subroutine box_image_write_info(this, iunit)
    class(box_image_t), intent(in) :: this
    integer,            intent(in) :: iunit

    PUSH_SUB(box_image_write_info)

    write(iunit,'(2x,3a,i6,a,i6)') 'Type = defined by image "', trim(this%filename), '"', this%image_size(1), ' x ', &
      this%image_size(2)

    POP_SUB(box_image_write_info)
  end subroutine box_image_write_info

  !--------------------------------------------------------------
  subroutine box_image_write_short_info(this, iunit)
    class(box_image_t), intent(in) :: this
    integer,            intent(in) :: iunit

    PUSH_SUB(box_image_write_short_info)

    write(iunit, '(2a)') 'BoxShape = box_image; BoxShapeImage = ', trim(this%filename)

    POP_SUB(box_image_write_short_info)
  end subroutine box_image_write_short_info

end module box_image_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
