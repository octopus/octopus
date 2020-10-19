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
  public ::           &
    compact,          &
    add_last_slash,   &
    str_trim,         &
    str_center,       &
    print_C_string,   &
    conv_to_C_string, &
    string_f_to_c,    &
    string_c_to_f,    &
    string_c_ptr_to_f

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

  ! Helper functions to convert between C and Fortran strings
  ! Based on the routines by Joseph M. Krahn

  ! ---------------------------------------------------------
  function string_f_to_c(f_string) result(c_string)
    character(len=*), intent(in) :: f_string
    character(kind=c_char,len=1) :: c_string(len_trim(f_string)+1)

    integer :: i, strlen

    strlen = len_trim(f_string)

    do i = 1, strlen
      c_string(i) = f_string(i:i)
    end do
    c_string(strlen+1) = C_NULL_CHAR

  end function string_f_to_c

  ! ---------------------------------------------------------
  subroutine string_c_to_f(c_string, f_string)
    character(kind=c_char,len=1), intent(in)  :: c_string(*)
    character(len=*),             intent(out) :: f_string

    integer :: i

    i = 1
    do while(c_string(i) /= C_NULL_CHAR .and. i <= len(f_string))
      f_string(i:i) = c_string(i)
      i = i + 1
    end do
    if (i < len(f_string)) f_string(i:) = ' '

  end subroutine string_c_to_f

  ! ---------------------------------------------------------
  subroutine string_c_ptr_to_f(c_string, f_string)
    type(c_ptr),      intent(in)  :: c_string
    character(len=*), intent(out) :: f_string

    character(len=1, kind=c_char), pointer :: p_chars(:)
    integer :: i

    if (.not. c_associated(c_string)) then
      f_string = ' '
    else
      call c_f_pointer(c_string, p_chars, [huge(0)])
      i = 1
      do while(p_chars(i) /= C_NULL_CHAR .and. i <= len(f_string))
        f_string(i:i) = p_chars(i)
        i = i + 1
      end do
      if (i < len(f_string)) f_string(i:) = ' '
    end if

  end subroutine string_c_ptr_to_f

end module string_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
