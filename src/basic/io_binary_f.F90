!! Copyright (C) 2009 X. Andrade
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
!! $Id: global.F90 3587 2007-11-22 16:43:00Z xavier $

#include "global.h"
#include "io_binary.h"

module io_binary_m
  use messages_m

  implicit none 

  private

  public ::             &
    io_binary_write,    &
    io_binary_read,     &
    io_binary_get_info

  interface io_binary_write
    module procedure swrite_binary, dwrite_binary, cwrite_binary, zwrite_binary, iwrite_binary, lwrite_binary
    module procedure iwrite_binary2, lwrite_binary2
  end interface

  interface io_binary_read
    module procedure sread_binary, dread_binary, cread_binary, zread_binary, iread_binary, lread_binary
    module procedure iread_binary2, lread_binary2
  end interface

contains

  subroutine swrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(4),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT

    call push_sub('io_binary_f.swrite_binary')

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine swrite_binary

  subroutine dwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    call push_sub('io_binary_f.dwrite_binary')

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine dwrite_binary

  subroutine cwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    call push_sub('io_binary_f.cwrite_binary')

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine cwrite_binary

  subroutine zwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    call push_sub('io_binary_f.zwrite_binary')

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine zwrite_binary

  subroutine iwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32

    call push_sub('io_binary_f.iwrite_binary')

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine iwrite_binary

  subroutine lwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64

    call push_sub('io_binary_f.lwrite_binary')

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine lwrite_binary

  subroutine iwrite_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(4),          intent(in)  :: ff(:, :)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32

    call push_sub('io_binary_f.iwrite_binary2')

    ierr = 0
    call write_binary(np, ff(1, 1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine iwrite_binary2

  subroutine lwrite_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(8),          intent(in)  :: ff(:, :)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64

    call push_sub('io_binary_f.lwrite_binary2')

    ierr = 0
    call write_binary(np, ff(1, 1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine lwrite_binary2

  !------------------------------------------------------

  subroutine sread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)   :: fname
    integer,             intent(in)   :: np
    real(4),             intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr

    integer, parameter :: type = TYPE_FLOAT

    call push_sub('io_binary_f.sread_binary')

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine sread_binary

  subroutine dread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(out) :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    call push_sub('io_binary_f.dread_binary')

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine dread_binary

  subroutine cread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    call push_sub('io_binary_f.cread_binary')

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine cread_binary

  subroutine zread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    call push_sub('io_binary_f.zread_binary')

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine zread_binary

  subroutine iread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(4),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32

    call push_sub('io_binary_f.iread_binary')

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine iread_binary

  subroutine lread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(8),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64

    call push_sub('io_binary_f.lread_binary')

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine lread_binary

  subroutine iread_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(4),          intent(out) :: ff(:, :)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32

    call push_sub('io_binary_f.iread_binary2')

    ierr = 0
    call read_binary(np, ff(1, 1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine iread_binary2

  subroutine lread_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(8),          intent(out) :: ff(:, :)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64

    call push_sub('io_binary_f.lread_binary2')

    ierr = 0
    call read_binary(np, ff(1, 1), type, ierr, trim(fname))

    call pop_sub()
  end subroutine lread_binary2

  subroutine io_binary_get_info(fname, np, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(inout) :: np
    integer,             intent(out)   :: ierr

    integer :: type
    
    call push_sub('io_binary_f.io_binary_get_info')

    type = 0
    ierr = 0
    call get_info_binary(np, type, ierr, trim(fname))

    call pop_sub()
  end subroutine io_binary_get_info

end module io_binary_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
