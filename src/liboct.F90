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

module liboct
  use global

  implicit none
  
  interface
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

  end interface

  private :: oct_parse_string, oct_parse_block_string
contains
  
  ! logical is a FORTRAN type, so we emulate the routine with integers
  subroutine oct_parse_logical(name, def, res)
    character(len=*), intent(IN) :: name
    logical, intent(in) :: def
    logical, intent(out) :: res
    
    integer :: idef, ires
    
    idef = 0
    if(def) idef = 1
    
    call oct_parse_int(C_String(name), idef, ires)
    res = (ires .ne. 0)
    
  end subroutine oct_parse_logical

  ! to avoid errors
  subroutine oct_parse_str(name, def, res)
    character(len=*), intent(in)  :: name, def
    character(len=*), intent(out) :: res

    call clear_str(res)
    call oct_parse_string(C_string(name), C_string(def), res)
  end subroutine oct_parse_str

  subroutine oct_parse_block_logical(name, l, c, res)
    character(len=*), intent(IN) :: name
    integer, intent(in)          :: l, c
    logical, intent(out)         :: res
    
    integer :: ires
    
    call oct_parse_block_int(C_String(name), l, c, ires)
    res = (ires .ne. 0)
    
  end subroutine oct_parse_block_logical

  subroutine oct_parse_block_str(name, l, c, res)
    character(len=*), intent(in)  :: name
    integer, intent(in) :: l, c
    character(len=*), intent(out) :: res

    character(len=80) :: name1

    name1 = C_string(name)
    call clear_str(res)
    call oct_parse_block_string(name1, l, c, res)
  end subroutine oct_parse_block_str

  ! we need this routine to communicate strings to C
  character(len=80) function C_string(si) result(so)
    character(len=*), intent(IN) :: si
    
    integer :: i
    
    i = len(trim(si))
    call clear_str(so)
    so(1:i) = si(1:i)
    so(i+1:i+1) = achar(0)
  end function C_string

  subroutine clear_str(str)
    character(len=*), intent(out) :: str

    integer :: i
    do i = 1, len(str)
      str(i:i) = " "
    end do
  end subroutine clear_str

end module liboct
