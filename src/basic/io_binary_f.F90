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

  implicit none 

  private

  public ::             &
    io_binary_write,    &
    io_binary_read

  interface io_binary_write
    module procedure swrite_binary, dwrite_binary, cwrite_binary, zwrite_binary
  end interface

  interface io_binary_read
    module procedure sread_binary, dread_binary, cread_binary, zread_binary
  end interface

contains

  subroutine swrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(4),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))
  end subroutine swrite_binary

  subroutine dwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))
  end subroutine dwrite_binary

  subroutine cwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))
  end subroutine cwrite_binary

  subroutine zwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))
  end subroutine zwrite_binary

  !------------------------------------------------------

  subroutine sread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)   :: fname
    integer,             intent(in)   :: np
    real(4),             intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr

    integer, parameter :: type = TYPE_FLOAT

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))
  end subroutine sread_binary

  subroutine dread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))
  end subroutine dread_binary

  subroutine cread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))
  end subroutine cread_binary

  subroutine zread_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    ierr = 0
    call read_binary(np, ff(1), type, ierr, trim(fname))
  end subroutine zread_binary

end module io_binary_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
