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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module varinfo_m
  use string_m
  use c_pointer_m

  implicit none

  private
  public ::               &
    varinfo_init,         &
    varinfo_end,          &
    varinfo_print,        &
    varinfo_search,       &
    varinfo_print_option, &
    varinfo_valid_option, &
    varinfo_option

  interface
    subroutine varinfo_init(filename)
      character(len=*), intent(in) :: filename
    end subroutine varinfo_init

    subroutine varinfo_end()
    end subroutine varinfo_end
  end interface


contains

  ! ---------------------------------------------------------
  subroutine varinfo_print(iunit, var, ierr)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,optional, intent(out):: ierr

    type(c_ptr) :: handle, opt, name, type, section, desc
    integer :: val
    logical :: first

    call varinfo_getvar(var, handle)
    if(.not. c_associated(handle)) then
      if(present(ierr)) then
        ierr = -1
      else
        write(iunit, '(3a)') 'Could not find a variable named ', var, '.'
        write(iunit, '(a)') 'You should update your variable info.'
      end if
      return
    end if

    if(present(ierr)) ierr = 0
    call varinfo_getinfo(handle, name, type, section, desc)
    call print_C_string(iunit, name, "Variable: ")
    call print_C_string(iunit, type, "Type:     ")
    call print_C_string(iunit, section, "Section:  ")
    write(iunit, '(a)') "Description:"
    call print_C_string(iunit, desc, "    ")

    call set_null(opt)
    first = .true.
    do
      call varinfo_getopt(handle, opt)
      if(.not. c_associated(opt)) then
        exit
      else
        if(first) then
          write(iunit, '(a)') "Available options:"
          first = .false.
        end if
        call varinfo_opt_getinfo(opt, name, val, desc)
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

    type(c_ptr) :: handle, opt, name, desc
    integer :: val, option_
    logical :: is_flag_

    is_flag_ = .false.
    if(present(is_flag)) is_flag_ = is_flag
    option_ = option ! copy that we can change

    l = .false.

    call varinfo_getvar(var, handle)
    if(.not. c_associated(handle)) return

    call set_null(opt)
    do
      call varinfo_getopt(handle, opt)
      if(.not. c_associated(opt)) exit
      call varinfo_opt_getinfo(opt, name, val, desc)

      if(is_flag_) then
        option_ = iand(option_, not(val))
      else
        if(val == option_) then
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

    type(c_ptr) :: handle, opt, name, desc
    integer :: val

    call varinfo_getvar(var, handle)
    if(.not. c_associated(handle)) return

    call set_null(opt)
    do
      call varinfo_getopt(handle, opt)
      if(.not. c_associated(opt)) exit

      call varinfo_opt_getinfo(opt, name, val, desc)

      if(val == option) then
        write(iunit, '(4a)', advance='no') "Input:", ' [', var, ' = '
        call print_C_string(iunit, name, advance='no')
        write(iunit, '(a)', advance='no') ']'
        if(present(pre)) then
          write(iunit, '(3a)') ' (', trim(pre), ')'
        else
          write(iunit, '(1x)')
        end if
        ! uncomment to print the description of the options
        !call print_C_string(iunit, desc, pre='  > ')

        return
      end if
    end do

  end subroutine varinfo_print_option

  ! ---------------------------------------------------------
  subroutine varinfo_search(iunit, var, ierr)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,optional, intent(out):: ierr
    
    type(c_ptr) :: handle, name, type, section, desc
    
    call set_null(handle)
    if(present(ierr)) ierr = -1
    do 
      call varinfo_search_var(var, handle)

      if(c_associated(handle)) then 
        if(present(ierr)) ierr = 0
      else
        exit
      end if

      call varinfo_getinfo(handle, name, type, section, desc)
      call print_C_string(iunit, name)
      
    end do
      
  end subroutine varinfo_search

  ! ---------------------------------------------------------
  integer function varinfo_option(var, option) result(val)
    character(len=*),  intent(in) :: var
    character(len=*),  intent(in) :: option

    type(c_ptr) :: handle
    integer  :: ierr

    call varinfo_getvar(var, handle)
    call varinfo_search_option(handle, option, val, ierr)

    if(ierr /= 0) then
      ! we cannot use messages here :-(
      stop "Error: invalid option"
    end if

  end function varinfo_option


end module varinfo_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
