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

module xml_file_m

  implicit none 

  private

  public ::                   &
    xml_file_t,               &
    xml_tag_t,                &
    xml_file_init,            &
    xml_file_end,             &
    xml_file_read,            &
    xml_file_tag,             &
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
      import xml_file_t
      type(xml_file_t), intent(out)   :: this
      character(len=*), intent(in)    :: filename
    end function xml_file_init

    subroutine xml_file_end(this)
      import xml_file_t
      type(xml_file_t), intent(inout) :: this
    end subroutine xml_file_end

    subroutine xml_file_tag(this, tag_name, tag)
      import xml_file_t
      import xml_tag_t
      type(xml_file_t), intent(inout) :: this
      character(len=*), intent(in)    :: tag_name
      type(xml_tag_t),  intent(out)   :: tag
    end subroutine xml_file_tag

    subroutine xml_tag_end(this)
      import xml_tag_t
      type(xml_tag_t), intent(inout) :: this
    end subroutine xml_tag_end

    integer function xml_tag_get_attribute(this, att_name, val)
      import xml_tag_t
      type(xml_tag_t),  intent(in)    :: this
      character(len=*), intent(in)    :: att_name
      integer,          intent(out)   :: val
    end function xml_tag_get_attribute

  end interface

 
  interface xml_file_read
    module procedure xml_file_read_integer
    module procedure xml_file_read_float
  end interface xml_file_read

contains

    integer function xml_file_read_integer(this, tag, val) result(ierr)
      type(xml_file_t), intent(inout) :: this
      character(len=*), intent(in)    :: tag
      integer,          intent(out)   :: val

      interface
        integer function xml_file_read_integer_low(this, tag, val) result(ierr)
          import xml_file_t
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
          import xml_file_t
          type(xml_file_t), intent(inout) :: this
          character(len=*), intent(in)    :: tag
          real(8),          intent(out)   :: val
        end function xml_file_read_double_low
      end interface

      ierr = xml_file_read_double_low(this, '<'//trim(tag), val)

    end function xml_file_read_float

end module xml_file_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
