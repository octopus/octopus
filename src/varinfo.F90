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
  public :: varinfo_init, &
            varinfo_end, &
            varinfo_print, &
            varinfo_print_option, &
            varinfo_valid_option

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
    integer :: value
    logical :: first

    call varinfo_getvar(var, handle)
    if(handle == 0) return

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
        call varinfo_opt_getinfo(opt, name, value, desc)
        call print_C_string(iunit, name, "  ")
        call print_C_string(iunit, desc, "    ")
      end if
    end do

  end subroutine varinfo_print

  logical function varinfo_valid_option(var, option) result(l)
    character(len=*), intent(in) :: var
    integer,          intent(in) :: option

    integer(POINTER_SIZE) :: handle, opt, name, desc
    integer :: value

    l = .false.

    call varinfo_getvar(var, handle)
    if(handle.eq.0) return

    opt = int(0, POINTER_SIZE)
    do
      call varinfo_getopt(handle, opt)
      if(opt == 0) exit
      call varinfo_opt_getinfo(opt, name, value, desc)
      if(value == option) then
        l = .true.; return
      end if
    end do    

  end function varinfo_valid_option

  ! ---------------------------------------------------------
  subroutine varinfo_print_option(iunit, var, option, pre)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,          intent(in) :: option
    character(len=*), intent(in) :: pre

    integer(POINTER_SIZE) :: handle, opt, name, desc
    integer :: value

    call varinfo_getvar(var, handle)
    if(handle.eq.0) return

    opt = int(0, POINTER_SIZE)
    do
      call varinfo_getopt(handle, opt)
      if(opt == 0) exit

      call varinfo_opt_getinfo(opt, name, value, desc)

      if(value == option) then
        write(iunit, '(a,a)', advance='no') pre, ' ['
        call print_C_string(iunit, name, advance='no')
        call print_C_string(iunit, desc, pre='] ')
        return
      end if
    end do

  end subroutine varinfo_print_option

end module varinfo
