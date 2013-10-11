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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: global.F90 3587 2007-11-22 16:43:00Z xavier $

#include "global.h"
#include "io_binary.h"

module io_binary_m
  use global_m
  use messages_m
  use mpi_m

  implicit none 

  private

  public ::                   &
    io_binary_write_header,   &
    io_binary_write,          &
    io_binary_write_parallel, &
    io_binary_read,           &
    io_binary_read_parallel,  &
    io_binary_get_info

  interface io_binary_write
    module procedure swrite_binary, dwrite_binary, cwrite_binary, zwrite_binary, iwrite_binary, lwrite_binary
    module procedure iwrite_binary2, lwrite_binary2, zwrite_binary3, cwrite_binary3, dwrite_binary3
  end interface io_binary_write
  
  interface io_binary_write_parallel
    module procedure swrite_parallel, dwrite_parallel, cwrite_parallel,  zwrite_parallel, iwrite_parallel, lwrite_parallel
  end interface io_binary_write_parallel

  interface io_binary_read
    module procedure sread_binary, dread_binary, cread_binary, zread_binary, iread_binary, lread_binary
    module procedure iread_binary2, lread_binary2, zread_binary3, cread_binary3, dread_binary3
  end interface io_binary_read

  interface io_binary_read_parallel
    module procedure sread_parallel, dread_parallel, cread_parallel,  zread_parallel, iread_parallel, lread_parallel
  end interface io_binary_read_parallel

  interface io_binary_write_header
    module procedure iwrite_header
  end interface io_binary_write_header

contains

  !> Interface to C to write the header of Integer type
  subroutine iwrite_header(fname, np_global, int_size, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np_global
    integer,             intent(in)  :: int_size
    integer,             intent(out) :: ierr
    
    PUSH_SUB(iwrite_header)

    ASSERT(np_global > 0)
    call write_header(np_global, int_size, ierr, trim(fname))

    POP_SUB(iwrite_header)
  end subroutine iwrite_header


  ! ------------------------------------------------------

  subroutine swrite_binary(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(4),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT

    PUSH_SUB(swrite_binary)

    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1), type, ierr, trim(fname))

    POP_SUB(zwrite_binary)
  end subroutine zwrite_binary

  !------------------------------------------------------

  subroutine dwrite_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(in)  :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    PUSH_SUB(dwrite_binary3)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(dwrite_binary3)
  end subroutine dwrite_binary3

  !------------------------------------------------------

  subroutine zwrite_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(zwrite_binary3)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(zwrite_binary3)
  end subroutine zwrite_binary3

  !------------------------------------------------------

  subroutine cwrite_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(in)  :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    PUSH_SUB(cwrite_binary3)

    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1, 1), type, ierr, trim(fname))

    POP_SUB(lwrite_binary2)
  end subroutine lwrite_binary2

  ! HASIERA **********************************************

  ! ------------------------------------------------------

  subroutine swrite_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    real(4),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status

    PUSH_SUB(swrite_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    !! I assume 4 byte real (ubuntu x86_64, gcc-4.4/4.6, openmpi) =>
    offset = (xlocal-1)*16+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_write_ordered(file_handle, ff(1), np, MPI_REAL4, status, mpi_err)
    ierr = mpi_err
#endif

    POP_SUB(swrite_parallel)
  end subroutine swrite_parallel

  !------------------------------------------------------

  subroutine dwrite_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    real(8),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status

    PUSH_SUB(dwrite_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0

    PUSH_SUB(cwrite_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    !! I assume 8 byte real (ubuntu x86_64, gcc-4.4/4.6, openmpi) => 
    offset = (xlocal-1)*8+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_write_ordered(file_handle, ff(1), np, MPI_REAL8, status, mpi_err)
    ierr = mpi_err
#endif

    POP_SUB(dwrite_parallel)
  end subroutine dwrite_parallel

  !------------------------------------------------------

  subroutine cwrite_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    complex(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status

    PUSH_SUB(cwrite_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    !! I assume 8 byte complex (ubuntu x86_64, gcc-4.4/4.6, openmpi) => MPI_DOUBLE_COMPLEX
    offset = (xlocal-1)*8+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_write_ordered(file_handle, ff(1), np, MPI_COMPLEX, status, mpi_err)
    ierr = mpi_err
#endif

    POP_SUB(cwrite_parallel)
  end subroutine cwrite_parallel

  !------------------------------------------------------

  subroutine zwrite_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status

    PUSH_SUB(zwrite_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    !! I assume 16 byte complex (ubuntu x86_64, gcc-4.4/4.6, openmpi) => MPI_DOUBLE_COMPLEX
    offset = (xlocal-1)*16+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_write_ordered(file_handle, ff(1), np, MPI_DOUBLE_COMPLEX, status, mpi_err)
    ierr = mpi_err
#endif

    POP_SUB(zwrite_parallel)
  end subroutine zwrite_parallel

  !------------------------------------------------------
  
  subroutine iwrite_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    integer(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status

    PUSH_SUB(iwrite_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    offset = (xlocal-1)*FC_INTEGER_SIZE+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_write_ordered(file_handle, ff(1), np, MPI_INTEGER, status, mpi_err)
    ierr = mpi_err
#endif

    POP_SUB(iwrite_parallel)
  end subroutine iwrite_parallel

  !------------------------------------------------------

  subroutine lwrite_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    integer(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status

    PUSH_SUB(lwrite_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    offset = (xlocal-1)*8+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_write_ordered(file_handle, ff(1), np, MPI_INTEGER4, status, mpi_err)
    ierr = mpi_err
#endif

    POP_SUB(lwrite_parallel)
  end subroutine lwrite_parallel

  !--- BUKAERA **************************************

  !------------------------------------------------------

  subroutine sread_binary(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)   :: fname
    integer,             intent(in)   :: np
    real(4),             intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr
    integer, optional,   intent(in)   :: offset

    integer, parameter :: type = TYPE_FLOAT

    PUSH_SUB(sread_binary)

    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

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

    ASSERT(np > 0)
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

    integer, parameter :: type = TYPE_DOUBLE

    PUSH_SUB(dread_binary3)
   
    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(dread_binary3)
  end subroutine dread_binary3

  !------------------------------------------------------

  subroutine zread_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(out) :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(zread_binary3)
   
    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1,1,1), type, ierr, trim(fname))

    POP_SUB(zread_binary3)
  end subroutine zread_binary3

  !------------------------------------------------------

  subroutine cread_binary3(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(4),          intent(out) :: ff(:,:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX

    PUSH_SUB(cread_binary3)
   
    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
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

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1, 1), type, ierr, trim(fname))

    POP_SUB(lread_binary2)
  end subroutine lread_binary2

  !------------------------------------------------------

  !--- HASIERA **************************************
  ! ------------------------------------------------------

  subroutine sread_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    real(4),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status, read_count

    PUSH_SUB(sread_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    !! I assume 4 byte real (ubuntu x86_64, gcc-4.4/4.6, openmpi) =>
    offset = (xlocal-1)*16+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_read(file_handle, ff(1), np, MPI_REAL4, status, mpi_err)
    call MPI_Get_count(status, MPI_REAL4, read_count, mpi_err)
    if (read_count /= np) then 
      write(message(1),'(a,i8,a,i8)') " read elements=", read_count, " instead of", np
      call messages_fatal(1)
    end if
    ierr = mpi_err
#endif

    POP_SUB(sread_parallel)
  end subroutine sread_parallel

  !------------------------------------------------------

  subroutine dread_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    real(8),             intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset 
#endif
    integer :: status, read_count

    PUSH_SUB(dread_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0

    PUSH_SUB(cread_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    !! I assume 8 byte real (ubuntu x86_64, gcc-4.4/4.6, openmpi) => 
    offset = (xlocal-1)*8+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_read(file_handle, ff(1), np, MPI_REAL8, status, mpi_err)
    call MPI_Get_count(status, MPI_REAL8, read_count, mpi_err)
    if (read_count /= np) then 
      write(message(1),'(a,i8,a,i8)') " read elements=", read_count, " instead of", np
      call messages_fatal(1)
    end if
    ierr = mpi_err
#endif

    POP_SUB(dread_parallel)
  end subroutine dread_parallel

  !------------------------------------------------------

  subroutine cread_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    complex(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_FLOAT_COMPLEX
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status, read_count

    PUSH_SUB(cread_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    !! I assume 8 byte complex (ubuntu x86_64, gcc-4.4/4.6, openmpi) => MPI_DOUBLE_COMPLEX
    offset = (xlocal-1)*8+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_read(file_handle, ff(1), np, MPI_COMPLEX, status, mpi_err)
    call MPI_Get_count(status, MPI_COMPLEX, read_count, mpi_err)
    if (read_count /= np) then 
      write(message(1),'(a,i8,a,i8)') " read elements=", read_count, " instead of", np
      call messages_fatal(1)
    end if
    ierr = mpi_err
#endif

    POP_SUB(cread_parallel)
  end subroutine cread_parallel

  !------------------------------------------------------

  subroutine zread_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status, read_count

    PUSH_SUB(zread_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    !! I assume 16 byte complex (ubuntu x86_64, gcc-4.4/4.6, openmpi) => MPI_DOUBLE_COMPLEX
    offset = (xlocal-1)*16+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_read(file_handle, ff(1), np, MPI_DOUBLE_COMPLEX, status, mpi_err)
    call MPI_Get_count(status, MPI_DOUBLE_COMPLEX, read_count, mpi_err)
    if (read_count /= np) then 
      write(message(1),'(a,i8,a,i8)') " read elements=", read_count, " instead of", np
      call messages_fatal(1)
    end if
    ierr = mpi_err
#endif

    POP_SUB(zread_parallel)
  end subroutine zread_parallel

  !------------------------------------------------------
  
  subroutine iread_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    integer(4),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_32
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status, read_count

    PUSH_SUB(iread_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    offset = (xlocal-1)*FC_INTEGER_SIZE+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_read(file_handle, ff(1), np, MPI_INTEGER, status, mpi_err)
    call MPI_Get_count(status, MPI_INTEGER, read_count, mpi_err)
    if (read_count /= np) then 
      write(message(1),'(a,i8,a,i8)') " read elements=", read_count, " instead of", np
      call messages_fatal(1)
    end if
    ierr = mpi_err
#endif

    POP_SUB(iread_parallel)
  end subroutine iread_parallel

  !------------------------------------------------------

  subroutine lread_parallel(file_handle, xlocal, np, ff, ierr)
    integer,             intent(in)  :: file_handle
    integer,             intent(in)  :: xlocal
    integer,             intent(in)  :: np
    integer(8),          intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_INT_64
#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset    
#endif
    integer :: status, read_count

    PUSH_SUB(lread_parallel)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
#ifdef HAVE_MPI2
    offset = (xlocal-1)*8+64
    call MPI_File_set_atomicity(file_handle, .true., mpi_err)
    call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    call MPI_File_read(file_handle, ff(1), np, MPI_INTEGER4, status, mpi_err)
    call MPI_Get_count(status, MPI_INTEGER4, read_count, mpi_err)
    if (read_count /= np) then 
      write(message(1),'(a,i8,a,i8)') " read elements=", read_count, " instead of", np
      call messages_fatal(1)
    end if
    ierr = mpi_err
#endif

    POP_SUB(lread_parallel)
  end subroutine lread_parallel

  !--- BUKAERA **************************************

  subroutine io_binary_get_info(fname, np, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(out)   :: np
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
