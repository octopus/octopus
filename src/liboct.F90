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

#include "config_F90.h"

module liboct

  implicit none
  
  interface
    real(8) function oct_gamma(x)
      real(8), intent(in) :: x
    end function oct_gamma

    real(8) function oct_bessel(n, x)
      integer,  intent(in) :: n
      real(8), intent(in)  :: x
    end function oct_bessel
    
    real(8) function oct_asinh(x)
      real(8), intent(in) :: x
    end function oct_asinh

    real(8) function oct_erf(x)
      real(8), intent(in) :: x
    end function oct_erf

    real(8) function oct_erfc(x)
      real(8), intent(in) :: x
    end function oct_erfc

    real(8) function oct_ylm(x, y, z, l, m)
      real(8), intent(in) :: x, y, z
      integer, intent(in) :: l, m
    end function oct_ylm

    subroutine oct_fft_optimize(n, p, par)
      integer, intent(inout) :: n
      integer, intent(in) :: p, par
    end subroutine oct_fft_optimize

    integer function oct_parse_init(file_in, file_out)
      character(len=*), intent(in)  :: file_in, file_out
    end function oct_parse_init

    subroutine oct_parse_end()
    end subroutine oct_parse_end

    integer function oct_parse_isdef(name)
      character(len=*), intent(in) :: name
    end function oct_parse_isdef

    subroutine oct_parse_int(name, def, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: def
      integer, intent(out)         :: res
    end subroutine oct_parse_int

    subroutine oct_parse_double(name, def, res)
      character(len=*), intent(in) :: name
      real(8), intent(in)          :: def
      real(8), intent(out)         :: res
    end subroutine oct_parse_double
    
    subroutine oct_parse_complex(name, def, res)
      character(len=*), intent(in) :: name
      complex(8), intent(in)       :: def
      complex(8), intent(out)      :: res
    end subroutine oct_parse_complex

    subroutine oct_parse_string(name, def, res)
      character(len=*), intent(in) :: name, def
      character(len=*), intent(out):: res
    end subroutine oct_parse_string

    integer function oct_parse_block_n(name)
      character(len=*), intent(in) :: name
    end function oct_parse_block_n

    subroutine oct_parse_block_int(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      integer, intent(out)         :: res
    end subroutine oct_parse_block_int

    subroutine oct_parse_block_double(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      real(8), intent(out)         :: res
    end subroutine oct_parse_block_double

    subroutine oct_parse_block_complex(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      complex(8), intent(out)      :: res
    end subroutine oct_parse_block_complex

    subroutine oct_parse_block_string(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      character(len=*), intent(out):: res
    end subroutine oct_parse_block_string

    real(8) function oct_parse_potential(x, y, z, r, pot)
      real(8), intent(in) :: x, y, z, r
      character(len=*), intent(in) :: pot
    end function oct_parse_potential

    subroutine oct_mkdir(name)
      character(len=*), intent(in) :: name
    end subroutine oct_mkdir

    subroutine oct_getcwd(name)
      character(len=*), intent(out) :: name
    end subroutine oct_getcwd

    subroutine oct_wfs_list(str, l)
      character(len=*), intent(in) :: str
      integer, intent(out) :: l(32)
    end subroutine oct_wfs_list

    real(8) function oct_ran_gaussian(r, sigma)
      integer(POINTER_SIZE), intent(in) :: r
      real(8), intent(in) :: sigma
    end function oct_ran_gaussian

    subroutine oct_ran_init(r)
      integer(POINTER_SIZE), intent(out) :: r
    end subroutine oct_ran_init

    subroutine oct_ran_end(r)
      integer(POINTER_SIZE), intent(out) :: r
    end subroutine oct_ran_end

    real(8) function oct_clock()
    end function oct_clock

    integer function oct_getmem()
    end function oct_getmem

    integer function print_file(filename)
      character(len=*), intent(in) :: filename
    end function print_file

    integer function number_of_lines(filename)
      character(len=*), intent(in) :: filename
    end function number_of_lines

    subroutine oct_sysname(name)
      character(len=*), intent(out) :: name
    end subroutine oct_sysname
  end interface

contains
  
  ! logical is a FORTRAN type, so we emulate the routine with integers
  subroutine oct_parse_logical(name, def, res)
    character(len=*), intent(IN) :: name
    logical, intent(in) :: def
    logical, intent(out) :: res
    
    integer :: idef, ires
    
    idef = 0
    if(def) idef = 1
    
    call oct_parse_int(name, idef, ires)
    res = (ires .ne. 0)
    
  end subroutine oct_parse_logical

  subroutine oct_parse_block_logical(name, l, c, res)
    character(len=*), intent(IN) :: name
    integer, intent(in)          :: l, c
    logical, intent(out)         :: res
    
    integer :: ires
    
    call oct_parse_block_int(name, l, c, ires)
    res = (ires .ne. 0)
    
  end subroutine oct_parse_block_logical

end module liboct
