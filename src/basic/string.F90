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

module string_oct_m
  
  use iso_c_binding
  use loct_oct_m
  
  implicit none

  private
  public ::          &
    compact,         &
    add_last_slash,  &
    str_trim,        &
    str_center,      &
    print_C_string,  &
    conv_to_C_string 

contains

  ! ---------------------------------------------------------
  !> Removes all spaces from a string
  !! \date 15-OCT-2000: First version
  !! \author Fernando Nogueira
  subroutine compact(str)
    character(len=*), intent(inout) :: str

    integer :: i, j

    j = 1
    do i = 1, len(str)
      if(str(i:i) /= ' ') then
        str(j:j) = str(i:i)
        j = j + 1
      end if
    end do
    do i = j, len(str)
      str(i:i) = ' '
    end do

  end subroutine compact

  ! ---------------------------------------------------------
  !> Adds a '/' in the end of the string, only if it missing.
  !! Useful for directories
  subroutine add_last_slash(str)
    character(len=*), intent(inout) :: str

    character(64) :: tmp_str
    
    if (index(str, '/', .true.) /= len_trim(str)) then
      tmp_str = str
      write(str,'(a,a1)') trim(tmp_str), '/'
    end if
  end subroutine add_last_slash


  ! ---------------------------------------------------------
  !> removes leading spaces from string
  subroutine str_trim(str)
    character (len=*), intent(inout) :: str
    integer :: i, j, l

    l = len(str)
    do i = 1, l
      if(str(i:i) /= ' ') exit
    end do

    do j = 1, l - i + 1
      str(j:j) = str(i:i)
      i = i + 1
    end do

    do i = j, l
      str(j:j) = ' '
    end do

  end subroutine str_trim

  ! ---------------------------------------------------------
  !> puts space around string, so that it is centered
  character(len=80) function str_center(s_in, l_in) result(s_out)
    character(len=*), intent(in) :: s_in
    integer,          intent(in) :: l_in

    integer :: pad, i, li, l

    l = min(80, l_in)
    li = len(s_in)
    if(l < li) then
      s_out(1:l) = s_in(1:l)
      return
    end if

    pad = (l - li)/2

    s_out = ""
    do i = 1, pad
      s_out(i:i) = " "
    end do
    s_out(pad + 1:pad + li) = s_in(1:li)
    do i = pad + li + 1, l
      s_out(i:i) = " "
    end do

  end function str_center

  ! ---------------------------------------------------------
  !> prints the C string given by the pointer str
  subroutine print_C_string(iunit, str, pre, advance)
    integer,                    intent(in) :: iunit
    type(c_ptr),                intent(in) :: str
    character(len=*), optional, intent(in) :: pre
    character(len=*), optional, intent(in) :: advance

    type(c_ptr)        :: s
    character(len=256) :: line
    character(len=5)   :: advance_

    advance_ = "yes"
    if(present(advance)) advance_ = advance

    s = c_null_ptr
    do
      call loct_break_C_string(str, s, line)
      if (.not. c_associated(s)) exit
      if(present(pre)) then
        write(iunit, '(a,a)', advance=advance_) pre, trim(line)
      else
        write(iunit, '(a)',   advance=advance_) trim(line)
      end if
    end do
  end subroutine print_C_string

  ! ---------------------------------------------------------
  !> converts to c string
  subroutine conv_to_C_string(str)
    character(len=*), intent(out) :: str
    
    integer :: j

    j = len(trim(str))
    str(j+1:j+1) = achar(0) 
  end subroutine conv_to_C_string


end module string_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
