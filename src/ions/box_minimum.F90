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

module box_minimum_oct_m
  use box_oct_m
  use box_shape_oct_m
  use global_oct_m
  use lookup_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  public :: box_minimum_t

  !> Class implementing a box that is a union of spheres. We do this in a specific class
  !! instead of using the box_union class for performance reasons (although this
  !! should be benchmarked at some point).
  type, extends(box_shape_t) :: box_minimum_t
    private
    FLOAT, public :: radius = M_ZERO !< Radius of the boxes when they all have the same radius.

    integer :: n_site_types !< How many sites types there are. Each site type has a label and a radius.
    character(len=:), allocatable :: site_type_label(:)  !< Array storing the labels of the site types
    FLOAT,            allocatable :: site_type_radius(:) !< Array storing the radii of the site types
    integer :: n_sites !< How many sites there are.
    integer, allocatable :: site_type(:)       !< Site type for each site.
    FLOAT,   allocatable :: site_position(:,:) !< Site coordinates.

    type(lookup_t) :: site_lookup
  contains
    procedure :: contains_points => box_minimum_contains_points
    procedure :: write_info => box_minimum_write_info
    procedure :: short_info => box_minimum_short_info
    final     :: box_minimum_finalize
  end type box_minimum_t

  interface box_minimum_t
    procedure box_minimum_constructor
  end interface box_minimum_t

contains

  !--------------------------------------------------------------
  function box_minimum_constructor(dim, n_site_types, site_type_label, site_type_radius, n_sites, site_type, &
    site_position) result(box)
    integer,                  intent(in) :: dim
    integer,                  intent(in) :: n_site_types
    character(len=*),         intent(in) :: site_type_label(1:n_site_types)
    FLOAT,                    intent(in) :: site_type_radius(1:n_site_types)
    integer,                  intent(in) :: n_sites
    integer,                  intent(in) :: site_type(1:n_sites)
    FLOAT,                    intent(in) :: site_position(1:dim,1:n_sites)
    class(box_minimum_t), pointer :: box

    FLOAT :: center(dim)

    PUSH_SUB(box_minimum_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    center = M_ZERO
    call box_shape_init(box, dim, center)
    box%n_site_types = n_site_types
    SAFE_ALLOCATE_SOURCE(box%site_type_label, site_type_label)
    SAFE_ALLOCATE_SOURCE(box%site_type_radius, site_type_radius)
    box%n_sites = n_sites
    SAFE_ALLOCATE_SOURCE(box%site_type, site_type)
    SAFE_ALLOCATE_SOURCE(box%site_position, site_position)
    if (all(site_type_radius(:) == site_type_radius(1))) then
      box%radius = site_type_radius(1)
    end if

    call lookup_init(box%site_lookup, box%dim, box%n_sites, box%site_position)

    POP_SUB(box_minimum_constructor)
  end function box_minimum_constructor

  !--------------------------------------------------------------
  subroutine box_minimum_finalize(this)
    type(box_minimum_t), intent(inout) :: this

    PUSH_SUB(box_minimum_finalize)

    call lookup_end(this%site_lookup)

    SAFE_DEALLOCATE_A(this%site_type_label)
    SAFE_DEALLOCATE_A(this%site_type_radius)
    SAFE_DEALLOCATE_A(this%site_type)
    SAFE_DEALLOCATE_A(this%site_position)

    call box_shape_end(this)

    POP_SUB(box_minimum_finalize)
  end subroutine box_minimum_finalize

  !--------------------------------------------------------------
  function box_minimum_contains_points(this, nn, xx) result(contained)
    class(box_minimum_t), intent(in)  :: this
    integer,              intent(in)  :: nn
    FLOAT,                intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip, site, ilist
    integer, allocatable :: nlist(:)
    integer, allocatable :: list(:, :)
    FLOAT :: max_radius, dist2

    max_radius = maxval(this%site_type_radius) + BOX_BOUNDARY_DELTA

    SAFE_ALLOCATE(nlist(1:nn))

    if (this%radius > M_ZERO) then
      call lookup_get_list(this%site_lookup, nn, xx, max_radius, nlist)
    else
      call lookup_get_list(this%site_lookup, nn, xx, max_radius, nlist, list = list)
    end if

    if (this%radius > M_ZERO) then
      do ip = 1, nn
        contained(ip) = nlist(ip) /= 0 .neqv. this%is_inside_out()
      end do
    else
      do ip = 1, nn
        contained(ip) = .false.
        do ilist = 1, nlist(ip)
          site = list(ilist, ip)
          dist2 = sum((xx(ip, 1:this%dim) - this%site_position(1:this%dim, site))**2)
          if (dist2 <= (this%site_type_radius(this%site_type(site)) + BOX_BOUNDARY_DELTA)**2) then
            contained(ip) = .true.
            exit
          end if
        end do
        contained(ip) = contained(ip) .neqv. this%is_inside_out()
      end do
    end if

    SAFE_DEALLOCATE_A(nlist)
    SAFE_DEALLOCATE_A(list)

  end function box_minimum_contains_points

  !--------------------------------------------------------------
  subroutine box_minimum_write_info(this, iunit)
    class(box_minimum_t), intent(in) :: this
    integer,             intent(in) :: iunit

    integer :: itype

    PUSH_SUB(box_minimum_write_info)

    write(message(1), '(2x,a)') 'Type = minimum'
    call messages_info(1, iunit)
    if (this%radius > M_ZERO) then
      write(message(1),'(2x,3a,f7.3)') 'Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, this%radius)
      call messages_info(1, iunit)
    else
      do itype = 1, this%n_site_types
        write(message(1),'(2x,a,a5,5x,a,f7.3,2a)') 'Species = ', trim(this%site_type_label(itype)), 'Radius = ', &
          units_from_atomic(units_out%length, this%site_type_radius(itype)), ' ', trim(units_abbrev(units_out%length))
        call messages_info(1, iunit)
      end do
    end if

    POP_SUB(box_minimum_write_info)
  end subroutine box_minimum_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_minimum_short_info(this, unit_length) result(info)
    class(box_minimum_t), intent(in) :: this
    type(unit_t),         intent(in) :: unit_length

    PUSH_SUB(box_minimum_short_info)

    write(info,'(a,f11.6,a,a)') 'BoxShape = minimum; Radius =', units_from_atomic(unit_length, this%radius), ' ', &
      trim(units_abbrev(unit_length))

    POP_SUB(box_minimum_short_info)
  end function box_minimum_short_info

end module box_minimum_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
