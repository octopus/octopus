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

#include "global.h"

module varinfo_oct_m
  use iso_c_binding
  use string_oct_m

  implicit none

  private
  public ::               &
    varinfo_init,         &
    varinfo_end,          &
    varinfo_print,        &
    varinfo_search,       &
    varinfo_print_option, &
    varinfo_valid_option, &
    varinfo_option,       &
    varinfo_exists

  interface varinfo_valid_option
    module procedure varinfo_valid_option_8
    module procedure varinfo_valid_option_4
  end interface varinfo_valid_option
  
  interface
    subroutine varinfo_init(filename)
      implicit none
      character(len=*), intent(in) :: filename
    end subroutine varinfo_init

    subroutine varinfo_getvar(name, var)
      use iso_c_binding
      implicit none
      character(len=*), intent(in)    :: name
      type(c_ptr),      intent(inout) :: var
    end subroutine varinfo_getvar

    subroutine varinfo_getinfo(var, name, type, default, section, desc)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in)  :: var
      type(c_ptr), intent(out) :: name
      type(c_ptr), intent(out) :: type
      type(c_ptr), intent(out) :: default
      type(c_ptr), intent(out) :: section
      type(c_ptr), intent(out) :: desc
    end subroutine varinfo_getinfo

    subroutine varinfo_getopt(var, opt)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in)  :: var
      type(c_ptr), intent(out) :: opt
    end subroutine varinfo_getopt

    subroutine varinfo_opt_getinfo(opt, name, val, desc)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in)  :: opt
      type(c_ptr), intent(out) :: name
      integer(8),  intent(out) :: val
      type(c_ptr), intent(out) :: desc
    end subroutine varinfo_opt_getinfo

    subroutine varinfo_search_var(name, var)
      use iso_c_binding
      implicit none
      character(len=*), intent(in)    :: name
      type(c_ptr),      intent(inout) :: var
    end subroutine varinfo_search_var

    subroutine varinfo_search_option(var, name, val, ierr)
      use iso_c_binding
      implicit none
      type(c_ptr),      intent(in)  :: var
      character(len=*), intent(in)  :: name
      integer,          intent(out) :: val
      integer,          intent(out) :: ierr
    end subroutine varinfo_search_option

    subroutine varinfo_end()
      implicit none
    end subroutine varinfo_end
  end interface


contains

  ! ---------------------------------------------------------
  subroutine varinfo_print(iunit, var, ierr)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,optional, intent(out):: ierr

    type(c_ptr) :: handle, opt, name, type, default, section, desc
    integer(8) :: val
    logical :: first

    call varinfo_getvar(var, handle)
    if(.not. c_associated(handle)) then
      if(present(ierr)) then
        ierr = -1
        return
      else
        write(iunit, '(3a)') 'ERROR: Could not find a variable named ', trim(var), '.'
        stop
      end if
    end if

    if(present(ierr)) ierr = 0
    call varinfo_getinfo(handle, name, type, default, section, desc)

    call print_C_string(iunit, name,    "Variable: ")
    call print_C_string(iunit, type,    "Type:     ")
    call print_C_string(iunit, default, "Default:  ")
    call print_C_string(iunit, section, "Section:  ")
    write(iunit, '(a)') "Description:"
    call print_C_string(iunit, desc, "    ")

    opt = c_null_ptr
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
  logical function varinfo_valid_option_8(var, option, is_flag) result(l)
    character(len=*),  intent(in) :: var
    integer(8),        intent(in) :: option
    logical, optional, intent(in) :: is_flag

    type(c_ptr) :: handle, opt, name, desc
    integer(8) :: val, option_
    logical :: is_flag_

    is_flag_ = .false.
    if(present(is_flag)) is_flag_ = is_flag
    option_ = option ! copy that we can change

    l = .false.

    call varinfo_getvar(var, handle)
    if(.not. c_associated(handle)) then
      write(0, '(3a)') 'ERROR: Could not find a variable named ', trim(var), '.'
      stop
    endif

    opt = c_null_ptr
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

  end function varinfo_valid_option_8

  ! ---------------------------------------------------------
  
  logical function varinfo_valid_option_4(var, option, is_flag) result(l)
    character(len=*),  intent(in) :: var
    integer,           intent(in) :: option
    logical, optional, intent(in) :: is_flag

    l = varinfo_valid_option_8(var, int(option, 8), is_flag)
    
  end function varinfo_valid_option_4

  ! ---------------------------------------------------------
  subroutine varinfo_print_option(iunit, var, option, pre)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,          intent(in) :: option
    character(len=*), intent(in), optional :: pre

    type(c_ptr) :: handle, opt, name, desc
    integer(8) :: val
    logical :: option_found
    
    call varinfo_getvar(var, handle)
    if(.not. c_associated(handle)) then
      write(iunit, '(3a)') 'ERROR: Could not find a variable named ', trim(var), '.'
      stop
    endif

    option_found = .false.
    opt = c_null_ptr
    do
      call varinfo_getopt(handle, opt)
      if(.not. c_associated(opt)) exit

      call varinfo_opt_getinfo(opt, name, val, desc)

      if(val == int(option, 8)) then
        option_found = .true.
        exit
      endif
    end do

    write(iunit, '(4a)', advance='no') "Input:", ' [', var, ' = '
    
    if(option_found) then
      call print_C_string(iunit, name, advance='no')
    else
      write(iunit,'(i6,a)', advance='no') option, " (INVALID)"
    endif
    write(iunit, '(a)', advance='no') ']'
    if(present(pre)) then
      write(iunit, '(3a)') ' (', trim(pre), ')'
    else
      write(iunit, '(1x)')
    end if
    ! uncomment to print the description of the options
    !call print_C_string(iunit, desc, pre='  > ')

    if(.not. option_found) then
      ! we cannot use messages here :-(
      write(iunit,'(a,i6,2a)') "ERROR: invalid option ", option, " for variable ", trim(var)
      stop
    end if
    
  end subroutine varinfo_print_option

  ! ---------------------------------------------------------
  subroutine varinfo_search(iunit, var, ierr)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,optional, intent(out):: ierr
    
    type(c_ptr) :: handle, name, type, default, section, desc

    handle = c_null_ptr
    if(present(ierr)) ierr = -1
    do 
      call varinfo_search_var(var, handle)

      if(c_associated(handle)) then 
        if(present(ierr)) ierr = 0
      else
        exit
      end if

      call varinfo_getinfo(handle, name, type, default, section, desc)
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
      write(0,'(4a)') "ERROR: invalid option ", trim(option), " for variable ", trim(var)
      stop
    end if

  end function varinfo_option

  ! ----------------------------------------------------------

  logical function varinfo_exists(var) result(exists)
    character(len=*),  intent(in) :: var

    type(c_ptr) :: handle

    handle = c_null_ptr
    call varinfo_search_var(var, handle)

    exists = c_associated(handle)

  end function varinfo_exists
  

end module varinfo_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
