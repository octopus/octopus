!! Copyright (C) 2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module lib_oct_parser
  implicit none

  ! Define the which routines can be seen from the outside
  private
  public :: loct_parse_init, loct_parse_end
  public :: loct_parse_isdef, loct_parse_int, loct_parse_float, loct_parse_cmplx, &
      loct_parse_string, loct_parse_logical
  public :: loct_parse_block_n, loct_parse_block_cols
  public :: loct_parse_block_int, loct_parse_block_float, loct_parse_block_cmplx, &
      loct_parse_block_string, loct_parse_block_logical
  public :: loct_parse_potential
  
  interface loct_parse_init
    integer function oct_parse_init(file_in, file_out)
      character(len=*), intent(in)  :: file_in, file_out
    end function oct_parse_init
  end interface

  interface loct_parse_end
    subroutine oct_parse_end()
    end subroutine oct_parse_end
  end interface

  interface loct_parse_isdef
    integer function oct_parse_isdef(name)
      character(len=*), intent(in) :: name
    end function oct_parse_isdef
  end interface

  interface loct_parse_int
    subroutine oct_parse_int(name, def, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: def
      integer, intent(out)         :: res
    end subroutine oct_parse_int
  end interface

  interface loct_parse_float
    module procedure oct_parse_double4
    
    subroutine oct_parse_double(name, def, res)
      character(len=*), intent(in)  :: name
      real(8),          intent(in)  :: def
      real(8),          intent(out) :: res
    end subroutine oct_parse_double
  end interface

  interface loct_parse_cmplx
    module procedure oct_parse_complex4

    subroutine oct_parse_complex(name, def, res)
      character(len=*), intent(in) :: name
      complex(8), intent(in)       :: def
      complex(8), intent(out)      :: res
    end subroutine oct_parse_complex
  end interface

  interface loct_parse_string
    subroutine oct_parse_string(name, def, res)
      character(len=*), intent(in) :: name, def
      character(len=*), intent(out):: res
    end subroutine oct_parse_string
  end interface

  interface loct_parse_block_n
    integer function oct_parse_block_n(name)
      character(len=*), intent(in) :: name
    end function oct_parse_block_n
  end interface

  interface loct_parse_block_cols
    integer function oct_parse_block_cols(name, line)
      character(len=*), intent(in) :: name
      integer, intent(in) :: line
    end function oct_parse_block_cols
  end interface

  interface loct_parse_block_int
    subroutine oct_parse_block_int(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      integer, intent(out)         :: res
    end subroutine oct_parse_block_int
  end interface

  interface loct_parse_block_float
    module procedure oct_parse_block_double4

    subroutine oct_parse_block_double(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      real(8), intent(out)         :: res
    end subroutine oct_parse_block_double
  end interface

  interface loct_parse_block_cmplx
    module procedure oct_parse_block_complex4

    subroutine oct_parse_block_complex(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      complex(8), intent(out)      :: res
    end subroutine oct_parse_block_complex
  end interface

  interface loct_parse_block_string
    subroutine oct_parse_block_string(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      character(len=*), intent(out):: res
    end subroutine oct_parse_block_string
  end interface

  interface loct_parse_potential
    module procedure oct_parse_potential4
    
    real(8) function oct_parse_potential(x, y, z, r, pot)
      real(8), intent(in) :: x, y, z, r
      character(len=*), intent(in) :: pot
    end function oct_parse_potential
  end interface

contains

  ! logical is a FORTRAN type, so we emulate the routine with integers
  subroutine loct_parse_logical(name, def, res)
    character(len=*), intent(IN) :: name
    logical, intent(in) :: def
    logical, intent(out) :: res
    
    integer :: idef, ires
    
    idef = 0
    if(def) idef = 1
    
    call oct_parse_int(name, idef, ires)
    res = (ires .ne. 0)
    
  end subroutine loct_parse_logical
  
  subroutine loct_parse_block_logical(name, l, c, res)
    character(len=*), intent(IN) :: name
    integer, intent(in)          :: l, c
    logical, intent(out)         :: res
    
    integer :: ires, r
    
    call oct_parse_block_int(name, l, c, ires)
    res = (ires .ne. 0)
    
  end subroutine loct_parse_block_logical

  ! The code may want to compile in single precision mode
  ! As I did not want to change the parser library, these
  ! driver functions just convert their arguments.
  subroutine oct_parse_double4(name, def4, res4)
    character(len=*), intent(in) :: name
    real(4), intent(in)          :: def4
    real(4), intent(out)         :: res4

    real(8) :: res8
    call oct_parse_double(name, real(def4, 8), res8)
    res4 = real(res8, kind=4)
  end subroutine oct_parse_double4

  subroutine oct_parse_complex4(name, def4, res4)
    character(len=*), intent(in) :: name
    complex(4), intent(in)       :: def4
    complex(4), intent(out)      :: res4

    complex(8) :: res8
    call oct_parse_complex(name, cmplx(def4, kind=8), res8)
    res4 = real(res8, kind=4)
  end subroutine oct_parse_complex4
 
  subroutine oct_parse_block_double4(name, l, c, res4)
    character(len=*), intent(in) :: name
    integer, intent(in)          :: l, c
    real(4), intent(out)         :: res4

    real(8) :: res8
    call oct_parse_block_double(name, l, c, res8)
    res4 = real(res8, kind=4)
  end subroutine oct_parse_block_double4

  subroutine oct_parse_block_complex4(name, l, c, res4)
    character(len=*), intent(in) :: name
    integer, intent(in)          :: l, c
    complex(4), intent(out)      :: res4

    complex(8) :: res8
    call oct_parse_block_complex(name, l, c, res8)
    res4 = cmplx(res8, kind=4)
  end subroutine oct_parse_block_complex4

  real(4) function oct_parse_potential4(x4, y4, z4, r4, pot)
    real(4), intent(in) :: x4, y4, z4, r4
    character(len=*), intent(in) :: pot

    real(8) :: res
    res = oct_parse_potential(real(x4, 8), real(y4, 8), real(z4, 8), real(r4, 8), pot)
    oct_parse_potential4 = real(res, 4)
  end function oct_parse_potential4

end module lib_oct_parser
