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
!! $Id$

#include "global.h"
#include "io_binary.h"

module io_binary_m
  use global_m
  use messages_m
  use mpi_m
  use profiling_m

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
    module procedure iwrite_binary2, lwrite_binary2, dwrite_binary2, zwrite_binary2
    module procedure zwrite_binary3, cwrite_binary3, dwrite_binary3
  end interface io_binary_write
  
  interface io_binary_write_parallel
    module procedure swrite_parallel, dwrite_parallel, cwrite_parallel,  zwrite_parallel, iwrite_parallel, lwrite_parallel
  end interface io_binary_write_parallel

  interface io_binary_read
    module procedure sread_binary, dread_binary, cread_binary, zread_binary, iread_binary, lread_binary
    module procedure iread_binary2, lread_binary2, zread_binary2, dread_binary2
    module procedure zread_binary3, cread_binary3, dread_binary3
  end interface io_binary_read

  interface io_binary_read_parallel
    module procedure sread_parallel, dread_parallel, cread_parallel,  zread_parallel, iread_parallel, lread_parallel
  end interface io_binary_read_parallel

  !> Interfaces to C to write the header
  interface io_binary_write_header
    module procedure swrite_header, dwrite_header, cwrite_header,  zwrite_header, iwrite_header, lwrite_header
  end interface io_binary_write_header

  interface
    subroutine get_info_binary(np, type, ierr, fname)
      integer,             intent(out)   :: np
      integer,             intent(out)   :: type
      integer,             intent(out)   :: ierr      
      character(len=*),    intent(in)    :: fname
    end subroutine get_info_binary
  end interface
    

contains

  ! ------------------------------------------------------

  subroutine swrite_header(fname, np_global, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np_global
    real(4),             intent(in)  :: ff
    integer,             intent(out) :: ierr
    
    integer, parameter :: type = TYPE_FLOAT
    
    PUSH_SUB(swrite_header)

    ASSERT(np_global > 0)
    call write_header(np_global, type, ierr, trim(fname))
    
    POP_SUB(swrite_header)
  end subroutine swrite_header

  ! ------------------------------------------------------

  subroutine dwrite_header(fname, np_global, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np_global
    real(8),             intent(in)  :: ff
    integer,             intent(out) :: ierr
    
    integer, parameter :: type = TYPE_DOUBLE
    
    PUSH_SUB(dwrite_header)

    ASSERT(np_global > 0)
    call write_header(np_global, type, ierr, trim(fname))
    
    POP_SUB(dwrite_header)
  end subroutine dwrite_header

  ! ------------------------------------------------------

  subroutine cwrite_header(fname, np_global, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np_global
    complex(4),          intent(in)  :: ff
    integer,             intent(out) :: ierr
    
    integer, parameter :: type = TYPE_FLOAT_COMPLEX
    
    PUSH_SUB(cwrite_header)

    ASSERT(np_global > 0)
    call write_header(np_global, type, ierr, trim(fname))
    
    POP_SUB(cwrite_header)
  end subroutine cwrite_header

  ! ------------------------------------------------------

  subroutine zwrite_header(fname, np_global, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np_global
    complex(8),          intent(in)  :: ff
    integer,             intent(out) :: ierr
    
    integer, parameter :: type = TYPE_DOUBLE_COMPLEX
    
    PUSH_SUB(zwrite_header)

    ASSERT(np_global > 0)
    call write_header(np_global, type, ierr, trim(fname))
    
    POP_SUB(zwrite_header)
  end subroutine zwrite_header

  ! ------------------------------------------------------

  subroutine iwrite_header(fname, np_global, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np_global
    integer(4),          intent(in)  :: ff
    integer,             intent(out) :: ierr
    
    integer, parameter :: type = TYPE_INT_32
    
    PUSH_SUB(iwrite_header)

    ASSERT(np_global > 0)
    call write_header(np_global, type, ierr, trim(fname))
    
    POP_SUB(iwrite_header)
  end subroutine iwrite_header

  ! ------------------------------------------------------

  subroutine lwrite_header(fname, np_global, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np_global
    integer(8),          intent(in)  :: ff
    integer,             intent(out) :: ierr
    
    integer, parameter :: type = TYPE_INT_64
    
    PUSH_SUB(lwrite_header)

    ASSERT(np_global > 0)
    call write_header(np_global, type, ierr, trim(fname))
    
    POP_SUB(lwrite_header)
  end subroutine lwrite_header

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

  !------------------------------------------------------

  subroutine dwrite_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(in)  :: ff(:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    PUSH_SUB(dwrite_binary2)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1,1), type, ierr, trim(fname))

    POP_SUB(dwrite_binary2)
  end subroutine dwrite_binary2

  !------------------------------------------------------

  subroutine zwrite_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(in)  :: ff(:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(zwrite_binary2)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call write_binary(np, ff(1,1), type, ierr, trim(fname))

    POP_SUB(zwrite_binary2)
  end subroutine zwrite_binary2

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

  ! ------------------------------------------------------

  subroutine io_binary_parallel_start(fname, file_handle, comm, xlocal, np, sizeof_ff, is_write, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(out)   :: file_handle
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    integer*8,           intent(in)    :: sizeof_ff !< should be same type as offset?
    logical,             intent(in)    :: is_write !< if false, is read.
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer(MPI_OFFSET_KIND) :: offset
#endif
    integer :: amode, mpi_info

    PUSH_SUB(io_binary_parallel_start)

    ASSERT(np > 0)

    ! FIXME: everything will fail if sizeof is not available!

    ierr = 0
#ifdef HAVE_MPI2
    offset = (xlocal-1)*sizeof_ff+64
    
    if(is_write) then
      amode = IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND)
    else
      amode = MPI_MODE_RDONLY
    endif
    call MPI_File_open(comm, fname, amode, MPI_INFO_NULL, file_handle, mpi_err)

    if(mpi_err == 0) then
      call MPI_File_set_atomicity(file_handle, .true., mpi_err)
      call MPI_File_seek(file_handle, offset, MPI_SEEK_SET, mpi_err)
    endif
    ierr = mpi_err
#else
    message(1) = "Internal error: cannot call io_binary parallel routines without MPI2."
    call messages_fatal(1)
#endif

    POP_SUB(io_binary_parallel_start)
  end subroutine io_binary_parallel_start

  ! ------------------------------------------------------

  subroutine io_binary_parallel_end(file_handle)
    integer, intent(out)   :: file_handle

    logical :: finalized

    PUSH_SUB(io_binary_parallel_end)

#ifdef HAVE_MPI2
    call MPI_Finalized(finalized, mpi_err)
    if (.not. finalized) then
      call MPI_File_close(file_handle, mpi_err)
    end if
#else
    message(1) = "Internal error: cannot call io_binary parallel routines without MPI2."
    call messages_fatal(1)
#endif

    ! how will mpi do error handling with mpi_err?
    ! do read_count here

    POP_SUB(io_binary_parallel_end)
  end subroutine io_binary_parallel_end


  ! ------------------------------------------------------

  subroutine swrite_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    real(4),             intent(in)    :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: file_handle

    PUSH_SUB(swrite_parallel)

    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .true., ierr)
    ASSERT(product(ubound(ff)) >= np)

#ifdef HAVE_MPI2
    if(ierr == 0) call MPI_File_write_ordered(file_handle, ff(1), np, MPI_REAL4, status, mpi_err)
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(swrite_parallel)
  end subroutine swrite_parallel

  !------------------------------------------------------

  subroutine dwrite_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    real(8),             intent(in)    :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: file_handle

    PUSH_SUB(dwrite_parallel)

    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .true., ierr)
    ASSERT(product(ubound(ff)) >= np)

#ifdef HAVE_MPI2
    if(ierr == 0) call MPI_File_write_ordered(file_handle, ff(1), np, MPI_REAL8, status, mpi_err)
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(dwrite_parallel)
  end subroutine dwrite_parallel

  !------------------------------------------------------

  subroutine cwrite_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    complex(4),          intent(in)    :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: file_handle

    PUSH_SUB(cwrite_parallel)

    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .true., ierr)
    ASSERT(product(ubound(ff)) >= np)

#ifdef HAVE_MPI2
    if(ierr == 0) call MPI_File_write_ordered(file_handle, ff(1), np, MPI_COMPLEX, status, mpi_err)
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(cwrite_parallel)
  end subroutine cwrite_parallel

  !------------------------------------------------------

  subroutine zwrite_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    complex(8),          intent(in)    :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: file_handle

    PUSH_SUB(zwrite_parallel)

    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .true., ierr)
    ASSERT(product(ubound(ff)) >= np)

#ifdef HAVE_MPI2
    if(ierr == 0) call MPI_File_write_ordered(file_handle, ff(1), np, MPI_DOUBLE_COMPLEX, status, mpi_err)
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(zwrite_parallel)
  end subroutine zwrite_parallel

  !------------------------------------------------------
  
  subroutine iwrite_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm   
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    integer(4),          intent(in)    :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: file_handle

    PUSH_SUB(iwrite_parallel)

    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .true., ierr)
    ASSERT(product(ubound(ff)) >= np)

#ifdef HAVE_MPI2
    if(ierr == 0) call MPI_File_write_ordered(file_handle, ff(1), np, MPI_INTEGER, status, mpi_err)
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(iwrite_parallel)
  end subroutine iwrite_parallel

  !------------------------------------------------------

  subroutine lwrite_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    integer(8),          intent(in)    :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif   
    integer :: file_handle

    PUSH_SUB(lwrite_parallel)

    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .true., ierr)
    ASSERT(product(ubound(ff)) >= np)

#ifdef HAVE_MPI2
    if(ierr == 0) call MPI_File_write_ordered(file_handle, ff(1), np, MPI_INTEGER4, status, mpi_err)
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(lwrite_parallel)
  end subroutine lwrite_parallel

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
    integer :: read_np, number_type
    real(8), allocatable :: read_ff(:)

    PUSH_SUB(zread_binary)

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0

    call get_info_binary(read_np, number_type, ierr, fname)
    ! if the type of the file is real, then read real numbers and convert to complex
    ! @TODO other casts are missing
    if (number_type /= type) then
      SAFE_ALLOCATE(read_ff(1:np))
      call dread_binary(fname, np, read_ff, ierr, offset)
      ff = read_ff
      SAFE_DEALLOCATE_A(read_ff)
    else
      call read_binary(np, optional_default(offset, 0), ff(1), type, ierr, trim(fname))
    endif
    
    POP_SUB(zread_binary)
  end subroutine zread_binary

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

  subroutine dread_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    real(8),             intent(out) :: ff(:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    PUSH_SUB(dread_binary2)
   
    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1,1), type, ierr, trim(fname))

    POP_SUB(dread_binary2)
  end subroutine dread_binary2

  !------------------------------------------------------

  subroutine zread_binary2(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(out) :: ff(:,:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX

    PUSH_SUB(zread_binary2)
   
    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = 0
    call read_binary(np, 0, ff(1,1), type, ierr, trim(fname))

    POP_SUB(zread_binary2)
  end subroutine zread_binary2

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

  subroutine sread_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    real(4),             intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: read_count, file_handle

    PUSH_SUB(sread_parallel)

    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .false., ierr)
    ASSERT(product(ubound(ff)) >= np)

#ifdef HAVE_MPI2
    call MPI_File_read(file_handle, ff(1), np, MPI_REAL4, status, mpi_err)
    call MPI_Get_count(status, MPI_REAL4, read_count, mpi_err)
    if (read_count /= np) then 
      write(message(1),'(a,i8,a,i8)') " real(4) read elements=", read_count, " instead of", np
      write(message(2), '(a,a)') " of file= ", fname
      call messages_fatal(2)
    end if
#endif
    
    call io_binary_parallel_end(file_handle)

    POP_SUB(sread_parallel)
  end subroutine sread_parallel

  !------------------------------------------------------

  subroutine dread_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    real(8),             intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: read_count, file_handle

    PUSH_SUB(dread_parallel)

    ASSERT(product(ubound(ff)) >= np)
    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .false., ierr)

#ifdef HAVE_MPI2
    if(ierr == 0) then
      call MPI_File_read(file_handle, ff(1), np, MPI_REAL8, status, mpi_err)
      call MPI_Get_count(status, MPI_REAL8, read_count, mpi_err)
      if (read_count /= np) then 
        write(message(1),'(a,i8,a,i8)') " real(8) read elements=", read_count, " instead of", np
        write(message(2), '(a,a)') " of file= ", fname
        call messages_fatal(2)
      end if
    endif
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(dread_parallel)
  end subroutine dread_parallel

  !------------------------------------------------------

  subroutine cread_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    complex(4),          intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: read_count, file_handle

    PUSH_SUB(cread_parallel)

    ASSERT(product(ubound(ff)) >= np)
    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .false., ierr)

#ifdef HAVE_MPI2
    if(ierr == 0) then
      call MPI_File_read(file_handle, ff(1), np, MPI_COMPLEX, status, mpi_err)
      call MPI_Get_count(status, MPI_COMPLEX, read_count, mpi_err)
      if (read_count /= np) then 
        write(message(1),'(a,i8,a,i8)') " complex(4) read elements=", read_count, " instead of", np
        write(message(2), '(a,a)') " of file= ", fname
        call messages_fatal(2)
      end if
    endif
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(cread_parallel)
  end subroutine cread_parallel

  !------------------------------------------------------

  subroutine zread_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    complex(8),          intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

    integer, parameter :: type = TYPE_DOUBLE_COMPLEX
#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: read_np, number_type
    integer :: read_count, file_handle
    real(8), allocatable :: read_ff(:)

    PUSH_SUB(zread_parallel)

    ASSERT(product(ubound(ff)) >= np)

    read_count = 0
#ifdef HAVE_MPI2
    call get_info_binary(read_np, number_type, ierr, fname)
    ! if the type of the file is real, then read real numbers and convert to complex
    ! @TODO other casts are missing
    if (number_type /= type) then
      if (in_debug_mode) then
        write(message(1),'(a,i2,a,i2)') "Debug: Found type = ", number_type, " instead of ", type
        call messages_info(1)
      end if
      SAFE_ALLOCATE(read_ff(1:np))
      call dread_parallel(fname, comm, xlocal, np, read_ff, ierr)
      ff = read_ff
      SAFE_DEALLOCATE_A(read_ff)
    else
      call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .false., ierr)
      if(ierr == 0) then
        call MPI_File_read(file_handle, ff(1), np, MPI_DOUBLE_COMPLEX, status, mpi_err)
        call MPI_Get_count(status, MPI_DOUBLE_COMPLEX, read_count, mpi_err)
        if (read_count /= np) then 
          write(message(1),'(a,i8,a,i8)') " complex(8) read elements=", read_count, " instead of", np
          write(message(2), '(a,a)') " of file= ", fname
          call messages_fatal(2)
        end if
      endif
    endif
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(zread_parallel)
  end subroutine zread_parallel

  !------------------------------------------------------
  
  subroutine iread_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    integer(4),          intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: read_count, file_handle

    PUSH_SUB(iread_parallel)

    ASSERT(product(ubound(ff)) >= np)
    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .false., ierr)

#ifdef HAVE_MPI2
    if(ierr == 0) then
      call MPI_File_read(file_handle, ff(1), np, MPI_INTEGER, status, mpi_err)
      call MPI_Get_count(status, MPI_INTEGER, read_count, mpi_err)
      if (read_count /= np) then 
        write(message(1),'(a,i8,a,i8)') " integer(4) read elements=", read_count, " instead of", np
        write(message(2), '(a,a)') " of file= ", fname
        call messages_fatal(2)
      end if
    endif
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(iread_parallel)
  end subroutine iread_parallel

  !------------------------------------------------------

  subroutine lread_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    integer(8),          intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: read_count, file_handle

    PUSH_SUB(lread_parallel)

    ASSERT(product(ubound(ff)) >= np)
    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, int(sizeof(ff(1)), kind=8), .false., ierr)

#ifdef HAVE_MPI2
    if(ierr == 0) then
      call MPI_File_read(file_handle, ff(1), np, MPI_INTEGER4, status, mpi_err)
      call MPI_Get_count(status, MPI_INTEGER4, read_count, mpi_err)
      if (read_count /= np) then 
        write(message(1),'(a,i8,a,i8)') " integer(8) read elements=", read_count, " instead of", np
        write(message(2), '(a,a)') " of file= ", fname
        call messages_fatal(2)
      end if
    endif
#endif

    call io_binary_parallel_end(file_handle)

    POP_SUB(lread_parallel)
  end subroutine lread_parallel

  !------------------------------------------------------

  subroutine io_binary_get_info(fname, np, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(out)   :: np
    integer,             intent(out)   :: ierr

    integer :: type
    
    PUSH_SUB(io_binary_get_info)

    ierr = 0
    call get_info_binary(np, type, ierr, trim(fname))

    POP_SUB(io_binary_get_info)
  end subroutine io_binary_get_info

end module io_binary_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
