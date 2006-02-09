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

module varinfo_m
  use string_m

  implicit none

  private
  public ::               &
    varinfo_init,         &
    varinfo_end,          &
    varinfo_print,        &
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


  ! ---------------------------------------------------------
  logical function varinfo_valid_option(var, option, is_flag) result(l)
    character(len=*),  intent(in) :: var
    integer,           intent(in) :: option
    logical, optional, intent(in) :: is_flag

    integer(POINTER_SIZE) :: handle, opt, name, desc
    integer :: value, option_
    logical :: is_flag_

    is_flag_ = .false.
    if(present(is_flag)) is_flag_ = is_flag
    option_ = option ! copy that we can change

    l = .false.

    call varinfo_getvar(var, handle)
    if(handle.eq.0) return

    opt = int(0, POINTER_SIZE)
    do
      call varinfo_getopt(handle, opt)
      if(opt == 0) exit
      call varinfo_opt_getinfo(opt, name, value, desc)

      if(is_flag_) then
        option_ = iand(option_, not(value))
      else
        if(value == option_) then
          l = .true.
          return
        end if
      end if

    end do

    if(is_flag_ .and. (option_ == 0)) l = .true.

  end function varinfo_valid_option


  ! ---------------------------------------------------------
  subroutine varinfo_print_option(iunit, var, option, pre)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,          intent(in) :: option
    character(len=*), intent(in), optional :: pre

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
        write(iunit, '(4a)', advance='no') "Input:", ' [', var, ' = '
        call print_C_string(iunit, name, advance='no')
        write(iunit, '(a)', advance='no') ']'
        if(present(pre)) then
          write(iunit, '(3a)') ' (', pre, ')'
        else
          write(iunit, '(1x)')
        end if
        ! uncoment to print the description of the options
        !call print_C_string(iunit, desc, pre='  > ')

        return
      end if
    end do

  end subroutine varinfo_print_option

end module varinfo_m
