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
  use global_m
  use messages_m

  implicit none 

  private

  public ::             &
    io_binary_write,    &
    io_binary_read,     &
    io_binary_get_info

  interface io_binary_write
    module procedure swrite_binary, dwrite_binary, cwrite_binary, zwrite_binary, iwrite_binary, lwrite_binary
    module procedure iwrite_binary2, lwrite_binary2, zwrite_binary3, cwrite_binary3, dwrite_binary3
  end interface

  interface io_binary_read
    module procedure sread_binary, dread_binary, cread_binary, zread_binary, iread_binary, lread_binary
    module procedure iread_binary2, lread_binary2, zread_binary3, cread_binary3, dread_binary3
  end interface

contains

  ! ------------------------------------------------------

  subroutine swrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(4),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT

    PUSH_SUB(swrite_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    POP_SUB(swrite_binary)
  end subroutine swrite_binary

  !------------------------------------------------------

  subroutine dwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    PUSH_SUB(dwrite_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    POP_SUB(dwrite_binary)
  end subroutine dwrite_binary

  !------------------------------------------------------

  subroutine cwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    PUSH_SUB(cwrite_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    POP_SUB(cwrite_binary)
  end subroutine cwrite_binary

  !------------------------------------------------------

  subroutine zwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(zwrite_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    POP_SUB(zwrite_binary)
  end subroutine zwrite_binary

  !------------------------------------------------------

  subroutine dwrite_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),          intent(in)  :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(dwrite_binary3)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(dwrite_binary3)
  end subroutine dwrite_binary3

  subroutine zwrite_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(zwrite_binary3)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(zwrite_binary3)
  end subroutine zwrite_binary3

  subroutine cwrite_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(in)  :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    PUSH_SUB(cwrite_binary3)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(cwrite_binary3)
  end subroutine cwrite_binary3

  !------------------------------------------------------
  
  subroutine iwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32

    PUSH_SUB(iwrite_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    POP_SUB(iwrite_binary)
  end subroutine iwrite_binary

  !------------------------------------------------------

  subroutine lwrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64

    PUSH_SUB(lwrite_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    POP_SUB(lwrite_binary)
  end subroutine lwrite_binary

  !------------------------------------------------------

  subroutine iwrite_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(4),          intent(in)  :: ff(:, :)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32

    PUSH_SUB(iwrite_binary2)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1, 1), type, ierr, trim(fname))

    POP_SUB(iwrite_binary2)
  end subroutine iwrite_binary2

  !------------------------------------------------------

  subroutine lwrite_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(8),          intent(in)  :: ff(:, :)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64

    PUSH_SUB(lwrite_binary2)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1, 1), type, ierr, trim(fname))

    POP_SUB(lwrite_binary2)
  end subroutine lwrite_binary2

  !------------------------------------------------------

  subroutine sread_binary(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)   :: fname
    integer,             intent(in)   :: np
    real(4),             intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr
    integer, optional,   intent(in)   :: offset

    integer, parameter :: type = TYPE_FLOAT

    PUSH_SUB(sread_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, optional_default(offset, 0), ff(1), type, ierr, trim(fname))

    POP_SUB(sread_binary)
  end subroutine sread_binary

  !------------------------------------------------------
 
  subroutine dread_binary(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(out) :: ff(:)
    integer,             intent(out) :: ierr
    integer, optional,   intent(in)  :: offset

    integer, parameter :: type = TYPE_DOUBLE

    PUSH_SUB(dread_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, optional_default(offset, 0), ff(1), type, ierr, trim(fname))

    POP_SUB(dread_binary)
  end subroutine dread_binary

  !------------------------------------------------------

  subroutine cread_binary(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr
    integer, optional,   intent(in)  :: offset

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    PUSH_SUB(cread_binary)

    ierr = 0
    call read_binary(np, optional_default(offset, 0), ff(1), type, ierr, trim(fname))

    POP_SUB(cread_binary)
  end subroutine cread_binary

  !------------------------------------------------------

  subroutine zread_binary(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr
    integer, optional,   intent(in)  :: offset

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(zread_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, optional_default(offset, 0), ff(1), type, ierr, trim(fname))
    
    POP_SUB(zread_binary)
  end subroutine zread_binary

  !------------------------------------------------------
  subroutine dread_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),          intent(out) :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(dread_binary3)
   
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(dread_binary3)
  end subroutine dread_binary3

  subroutine zread_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(out) :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(zread_binary3)
   
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(zread_binary3)
  end subroutine zread_binary3

  subroutine cread_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(out) :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    PUSH_SUB(cread_binary3)
   
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(cread_binary3)
  end subroutine cread_binary3

  !------------------------------------------------------
 
  subroutine iread_binary(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(4),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr
    integer, optional,   intent(in)  :: offset

    integer, parameter :: type = TYPE_INT_32

    PUSH_SUB(iread_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, optional_default(offset, 0), ff(1), type, ierr, trim(fname))

    POP_SUB(iread_binary)
  end subroutine iread_binary

  !------------------------------------------------------

  subroutine lread_binary(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(8),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr
    integer, optional,   intent(in)  :: offset

    integer, parameter :: type = TYPE_INT_64

    PUSH_SUB(lread_binary)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, optional_default(offset, 0), ff(1), type, ierr, trim(fname))

    POP_SUB(lread_binary)
  end subroutine lread_binary

  !------------------------------------------------------

  subroutine iread_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(4),          intent(out) :: ff(:, :)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32

    PUSH_SUB(iread_binary2)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1, 1), type, ierr, trim(fname))

    POP_SUB(iread_binary2)
  end subroutine iread_binary2

  !------------------------------------------------------

  subroutine lread_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    integer(8),          intent(out) :: ff(:, :)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64

    PUSH_SUB(lread_binary2)

    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1, 1), type, ierr, trim(fname))

    POP_SUB(lread_binary2)
  end subroutine lread_binary2

  !------------------------------------------------------

  subroutine io_binary_get_info(fname, np, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(inout) :: np
    integer,             intent(out)   :: ierr

    integer :: type
    
    PUSH_SUB(io_binary_get_info)

    type = 0
    ierr = 0
    call get_info_binary(np, type, ierr, trim(fname))

    POP_SUB(io_binary_get_info)
  end subroutine io_binary_get_info

end module io_binary_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
