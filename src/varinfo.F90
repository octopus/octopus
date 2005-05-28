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

module varinfo
  use string

  implicit none

  private
  public :: varinfo_init, varinfo_end, varinfo_print

  interface
    subroutine varinfo_init(filename)
      character(len=*), intent(in) :: filename
    end subroutine varinfo_init

    subroutine varinfo_end()
    end subroutine varinfo_end
  end interface
  
contains

  ! ---------------------------------------------------------
  subroutine varinfo_print(iunit, var)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var

    integer(POINTER_SIZE) :: handle, opt, name, type, section, desc
    logical :: first

    call varinfo_getvar(var, handle)
    if(handle.ne.0) then
      call varinfo_getinfo(handle, name, type, section, desc)
      call print_C_string(iunit, name, "Variable: ")
      call print_C_string(iunit, type, "Type:     ")
      call print_C_string(iunit, section, "Section:  ")
      write(iunit, '(a)') "Description:"
      call print_C_string(iunit, desc, "    ")

      opt = int(0, POINTER_SIZE)
      first = .true.
      do
        call varinfo_getopt(handle, opt)
        if(opt==0) then
          exit
        else
          if(first) then
            write(iunit, '(a)') "Available options:"
            first = .false.
          end if
          call varinfo_opt_getinfo(opt, name, desc)
          call print_C_string(iunit, name, "  ")
          call print_C_string(iunit, desc, "    ")
        end if
      end do

    end if
  end subroutine varinfo_print

end module varinfo
