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
    module procedure dwrite_binary, zwrite_binary, iwrite_binary
    module procedure iwrite_binary2, dwrite_binary2, zwrite_binary2
    module procedure zwrite_binary3, dwrite_binary3, iwrite_binary3
  end interface io_binary_write
  
  interface io_binary_write_parallel
    module procedure dwrite_parallel, zwrite_parallel, iwrite_parallel
  end interface io_binary_write_parallel

  interface io_binary_read
    module procedure dread_binary, zread_binary, iread_binary
    module procedure iread_binary2, zread_binary2, dread_binary2
    module procedure zread_binary3, iread_binary3, dread_binary3
  end interface io_binary_read

  interface io_binary_read_parallel
    module procedure dread_parallel, zread_parallel, iread_parallel
  end interface io_binary_read_parallel

  !> Interfaces to C to write the header
  interface io_binary_write_header
    module procedure dwrite_header, zwrite_header, iwrite_header
  end interface io_binary_write_header

  interface
    subroutine get_info_binary(np, type, ierr, fname)
      integer,             intent(out)   :: np
      integer,             intent(out)   :: type
      integer,             intent(out)   :: ierr      
      character(len=*),    intent(in)    :: fname
    end subroutine get_info_binary
  end interface

  interface
    subroutine write_header(np_global, type, ierr, fname)
      integer,             intent(in)  :: np_global
      integer,             intent(in)  :: type
      integer,             intent(out) :: ierr      
      character(len=*),    intent(in)  :: fname
    end subroutine write_header
  end interface

  ! no interfaces for read_binary, write_binary since we call them with different types

contains

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

  subroutine try_dread_binary(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    complex(8),          intent(out) :: ff(:)
    integer,             intent(out) :: ierr
    integer, optional,   intent(in)  :: offset

    integer :: read_np, number_type
    real(8), allocatable :: read_ff(:)

    PUSH_SUB(try_dread_binary)

    call get_info_binary(read_np, number_type, ierr, fname)
    ! if the type of the file is real, then read real numbers and convert to complex
    ! @TODO other casts are missing
    if (number_type /= TYPE_DOUBLE_COMPLEX) then
      if (in_debug_mode) then
        write(message(1),'(a,i2,a,i2)') "Debug: Found type = ", number_type, " instead of ", TYPE_DOUBLE_COMPLEX
        call messages_info(1)
      end if

      SAFE_ALLOCATE(read_ff(1:np))
      call dread_binary(fname, np, read_ff, ierr, offset)
      ff = read_ff
      SAFE_DEALLOCATE_A(read_ff)
    else
      ierr = -1
    endif
    ! ierr will be 0 if dread_binary succeeded

    POP_SUB(try_dread_binary)
  end subroutine try_dread_binary

  !------------------------------------------------------

  subroutine try_dread_parallel(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    complex(8),          intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

    integer :: read_np, number_type
    real(8), allocatable :: read_ff(:)

    PUSH_SUB(try_dread_parallel)

    call get_info_binary(read_np, number_type, ierr, fname)
    ! if the type of the file is real, then read real numbers and convert to complex
    ! @TODO other casts are missing
    if (number_type /= TYPE_DOUBLE_COMPLEX) then
      if (in_debug_mode) then
        write(message(1),'(a,i2,a,i2)') "Debug: Found type = ", number_type, " instead of ", TYPE_DOUBLE_COMPLEX
        call messages_info(1)
      end if
      SAFE_ALLOCATE(read_ff(1:np))
      call dread_parallel(fname, comm, xlocal, np, read_ff, ierr)
      ff = read_ff
      SAFE_DEALLOCATE_A(read_ff)
    else
      ierr = -1
    endif
    ! ierr will be 0 if dread_parallel succeeded

    POP_SUB(try_dread_parallel)
  end subroutine try_dread_parallel

  !------------------------------------------------------

  subroutine io_binary_get_info(fname, np, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(out)   :: np
    integer,             intent(out)   :: ierr

    integer :: type
    
    PUSH_SUB(io_binary_get_info)

    call get_info_binary(np, type, ierr, trim(fname))

    POP_SUB(io_binary_get_info)
  end subroutine io_binary_get_info

#include "complex.F90"
#include "io_binary_f_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "io_binary_f_inc.F90"

#include "undef.F90"

#include "integer.F90"
#include "io_binary_f_inc.F90"

end module io_binary_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
