!! Copyright (C) 2015 X. Andrade
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
!! $Id:$

#include "global.h"

module xml_oct_m

  ! This module implements a very limited and non-general interface to
  ! read xml files. Its main purpose it is to read pseudo-potentials
  ! in xml format.

  implicit none 

  private

  public ::                        &
    xml_file_t,                    &
    xml_tag_t,                     &
    xml_file_init,                 &
    xml_file_end,                  &
    xml_get_tag_value,             &
    xml_file_tag,                  &
    xml_tag_get_attribute_value,   &
    xml_tag_get_attribute_float,   &
    xml_tag_get_attribute_string,  &
    xml_tag_end

  type xml_file_t
    private
    integer, pointer :: dummy
  end type xml_file_t

  type xml_tag_t
    private
    integer, pointer :: dummy
  end type xml_tag_t

  interface
    integer function xml_file_init(this, filename)
      import :: xml_file_t
      type(xml_file_t), intent(out)   :: this
      character(len=*), intent(in)    :: filename
    end function xml_file_init

    subroutine xml_file_end(this)
      import :: xml_file_t
      type(xml_file_t), intent(inout) :: this
    end subroutine xml_file_end

    integer function xml_file_tag(this, tag_name, index, tag)
      import :: xml_file_t
      import :: xml_tag_t
      type(xml_file_t), intent(inout) :: this
      character(len=*), intent(in)    :: tag_name
      integer,          intent(in)    :: index
      type(xml_tag_t),  intent(out)   :: tag
    end function xml_file_tag

    subroutine xml_tag_end(this)
      import :: xml_tag_t
      type(xml_tag_t), intent(inout) :: this
    end subroutine xml_tag_end

    integer function xml_tag_get_attribute_value(this, att_name, val)
      import :: xml_tag_t
      type(xml_tag_t),  intent(in)    :: this
      character(len=*), intent(in)    :: att_name
      integer,          intent(out)   :: val
    end function xml_tag_get_attribute_value

    integer function xml_tag_get_attribute_float(this, att_name, val)
      import :: xml_tag_t
      type(xml_tag_t),  intent(in)    :: this
      character(len=*), intent(in)    :: att_name
      real(8),          intent(out)   :: val
    end function xml_tag_get_attribute_float
    
    integer function xml_tag_get_attribute_string(this, att_name, val)
      import :: xml_tag_t
      type(xml_tag_t),  intent(in)    :: this
      character(len=*), intent(in)    :: att_name
      character(len=*), intent(out)   :: val
    end function xml_tag_get_attribute_string
    
  end interface

  interface xml_get_tag_value
    module procedure xml_file_read_integer
    module procedure xml_file_read_float
    module procedure xml_file_read_string
    module procedure xml_tag_get_tag_value_array
    module procedure xml_tag_get_value_array
  end interface xml_get_tag_value

contains

  integer function xml_file_read_integer(this, tag, val) result(ierr)
    type(xml_file_t), intent(inout) :: this
    character(len=*), intent(in)    :: tag
    integer,          intent(out)   :: val

    interface
      integer function xml_file_read_integer_low(this, tag, val) result(ierr)
        import :: xml_file_t
        type(xml_file_t), intent(inout) :: this
        character(len=*), intent(in)    :: tag
        integer,          intent(out)   :: val
      end function xml_file_read_integer_low
    end interface

    ierr = xml_file_read_integer_low(this, '<'//trim(tag), val)

  end function xml_file_read_integer

  ! ----------------------------------------------------------------------

  integer function xml_file_read_float(this, tag, val) result(ierr)
    type(xml_file_t), intent(inout) :: this
    character(len=*), intent(in)    :: tag
    real(8),          intent(out)   :: val

    interface
      integer function xml_file_read_double_low(this, tag, val) result(ierr)
        import :: xml_file_t
        type(xml_file_t), intent(inout) :: this
        character(len=*), intent(in)    :: tag
        real(8),          intent(out)   :: val
      end function xml_file_read_double_low
    end interface

    ierr = xml_file_read_double_low(this, '<'//trim(tag), val)

  end function xml_file_read_float

  ! ----------------------------------------------------------------------

  integer function xml_file_read_string(this, tag, val) result(ierr)
    type(xml_file_t), intent(inout) :: this
    character(len=*), intent(in)    :: tag
    character(len=*), intent(out)   :: val

    interface
      integer function xml_file_read_string_low(this, tag, val) result(ierr)
        import :: xml_file_t
        type(xml_file_t), intent(inout) :: this
        character(len=*), intent(in)    :: tag
        character(len=*), intent(out)   :: val
      end function xml_file_read_string_low
    end interface

    ierr = xml_file_read_string_low(this, '<'//trim(tag), val)

  end function xml_file_read_string

  ! ----------------------------------------------------------------------

  integer function xml_tag_get_tag_value_array(this, tag_name, size, val) result(ierr)
    type(xml_tag_t),  intent(in)    :: this
    character(len=*), intent(in)    :: tag_name
    integer,          intent(in)    :: size
    real(8),          intent(out)   :: val(:)

    interface
      integer function xml_tag_get_tag_value_array_low(this, tag_name, size, val)
        import :: xml_tag_t
        type(xml_tag_t),  intent(in)    :: this
        character(len=*), intent(in)    :: tag_name
        integer,          intent(in)    :: size
        real(8),          intent(out)   :: val
      end function xml_tag_get_tag_value_array_low
    end interface

    if(size > 0) then
      ierr = xml_tag_get_tag_value_array_low(this, tag_name, size, val(1))
    else
      ierr = 0
    end if

  end function xml_tag_get_tag_value_array
  ! ----------------------------------------------------------------------

  integer function xml_tag_get_value_array(this, size, val) result(ierr)
    type(xml_tag_t),  intent(in)    :: this
    integer,          intent(in)    :: size
    real(8),          intent(out)   :: val(:)

    interface
      integer function xml_tag_get_value_array_low(this, size, val)
        import :: xml_tag_t
        type(xml_tag_t),  intent(in)    :: this
        integer,          intent(in)    :: size
        real(8),          intent(out)   :: val
      end function xml_tag_get_value_array_low
    end interface

    if(size > 0) then
      ierr = xml_tag_get_value_array_low(this, size, val(1))
    else
      ierr = 0
    end if

  end function xml_tag_get_value_array
  
end module xml_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
