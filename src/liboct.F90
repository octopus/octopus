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

#include "global.h"

module liboct

  implicit none
  
  interface
    FLOAT function oct_gamma(x)
      FLOAT, intent(in) :: x
    end function oct_gamma

    FLOAT function oct_bessel(n, x)
      integer,  intent(in) :: n
      FLOAT, intent(in)  :: x
    end function oct_bessel
    
    FLOAT function oct_asinh(x)
      FLOAT, intent(in) :: x
    end function oct_asinh

    FLOAT function oct_erf(x)
      FLOAT, intent(in) :: x
    end function oct_erf

    FLOAT function oct_erfc(x)
      FLOAT, intent(in) :: x
    end function oct_erfc

    FLOAT function oct_ylm(x, y, z, l, m)
      FLOAT, intent(in) :: x, y, z
      integer, intent(in) :: l, m
    end function oct_ylm

    subroutine oct_fft_optimize(n, p, par)
      integer, intent(inout) :: n
      integer, intent(in) :: p, par
    end subroutine oct_fft_optimize

    subroutine oct_mkdir(name)
      character(len=*), intent(in) :: name
    end subroutine oct_mkdir

    subroutine oct_rm(name)
      character(len=*), intent(in) :: name
    end subroutine oct_rm
    
    subroutine oct_getcwd(name)
      character(len=*), intent(out) :: name
    end subroutine oct_getcwd

    subroutine oct_wfs_list(str, l)
      character(len=*), intent(in) :: str
      integer, intent(out) :: l(32)
    end subroutine oct_wfs_list

    FLOAT function oct_ran_gaussian(r, sigma)
      integer(POINTER_SIZE), intent(in) :: r
      FLOAT, intent(in) :: sigma
    end function oct_ran_gaussian

    subroutine oct_ran_init(r)
      integer(POINTER_SIZE), intent(out) :: r
    end subroutine oct_ran_init

    subroutine oct_ran_end(r)
      integer(POINTER_SIZE), intent(out) :: r
    end subroutine oct_ran_end

    FLOAT function oct_clock()
    end function oct_clock

    integer function oct_getmem()
    end function oct_getmem

    subroutine oct_sysname(name)
      character(len=*), intent(out) :: name
    end subroutine oct_sysname

    integer function print_file(filename)
      character(len=*), intent(in) :: filename
    end function print_file

    integer function number_of_lines(filename)
      character(len=*), intent(in) :: filename
    end function number_of_lines

  end interface

end module liboct
