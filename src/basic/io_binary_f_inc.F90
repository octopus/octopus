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

  subroutine X(write_header)(fname, np_global, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np_global
    R_TYPE,              intent(in)  :: ff
    integer,             intent(out) :: ierr
    
    PUSH_SUB(X(write_header))

    call write_header(np_global, R_TYPE_IOBINARY, ierr, trim(fname))
    
    POP_SUB(X(write_header))
  end subroutine X(write_header)

  ! ------------------------------------------------------

  subroutine X(write_binary)(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    R_TYPE,              intent(in)  :: ff(:)
    integer,             intent(out) :: ierr

    PUSH_SUB(X(write_binary))

    ASSERT(product(ubound(ff)) >= np)

    call write_binary(np, ff(1), R_TYPE_IOBINARY, ierr, trim(fname))

    POP_SUB(X(write_binary))
  end subroutine X(write_binary)

  !------------------------------------------------------

  subroutine X(write_binary2)(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    R_TYPE,              intent(in)  :: ff(:, :)
    integer,             intent(out) :: ierr

    PUSH_SUB(X(write_binary2))

    ASSERT(product(ubound(ff)) >= np)

    call write_binary(np, ff(1, 1), R_TYPE_IOBINARY, ierr, trim(fname))

    POP_SUB(X(write_binary2))
  end subroutine X(write_binary2)

  !------------------------------------------------------

  subroutine X(write_binary3)(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    R_TYPE,              intent(in)  :: ff(:,:,:)
    integer,             intent(out) :: ierr

    PUSH_SUB(X(write_binary3))

    ASSERT(product(ubound(ff)) >= np)

    call write_binary(np, ff(1,1,1), R_TYPE_IOBINARY, ierr, trim(fname))

    POP_SUB(X(write_binary3))
  end subroutine X(write_binary3)

  ! ------------------------------------------------------

  subroutine X(write_parallel)(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    R_TYPE,              intent(in)    :: ff(:)
    integer,             intent(out)   :: ierr
 
 #ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
 #endif
    integer :: file_handle
 
    PUSH_SUB(X(write_parallel))
 
    call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, R_SIZEOF, .true., ierr)
    ASSERT(product(ubound(ff)) >= np)
 
 #ifdef HAVE_MPI2
    if(ierr == 0) call MPI_File_write_ordered(file_handle, ff(1), np, R_MPITYPE, status, mpi_err)
 #endif
 
    call io_binary_parallel_end(file_handle)
 
    POP_SUB(X(write_parallel))
  end subroutine X(write_parallel)

  !------------------------------------------------------

  subroutine X(read_binary)(fname, np, ff, ierr, offset)
    character(len=*),    intent(in)   :: fname
    integer,             intent(in)   :: np
    R_TYPE,              intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr
    integer, optional,   intent(in)   :: offset

    PUSH_SUB(X(read_binary))

    ASSERT(np > 0)
    ASSERT(product(ubound(ff)) >= np)

    ierr = -1
#ifdef R_TCOMPLEX
    call try_dread_binary(fname, np, ff, ierr, offset)
#endif
    if(ierr == -1) &
      call read_binary(np, optional_default(offset, 0), ff(1), R_TYPE_IOBINARY, ierr, trim(fname))

    POP_SUB(X(read_binary))
  end subroutine X(read_binary)

  !------------------------------------------------------ 

  subroutine X(read_parallel)(fname, comm, xlocal, np, ff, ierr)
    character(len=*),    intent(in)    :: fname
    integer,             intent(in)    :: comm
    integer,             intent(in)    :: xlocal
    integer,             intent(in)    :: np
    R_TYPE,              intent(inout) :: ff(:)
    integer,             intent(out)   :: ierr

#ifdef HAVE_MPI2
    integer :: status(MPI_STATUS_SIZE)
#endif
    integer :: read_count, file_handle

    PUSH_SUB(X(read_parallel))

    ASSERT(np > 0)
    ierr = -1
#ifdef R_TCOMPLEX
    call try_dread_parallel(fname, comm, xlocal, np, ff, ierr)
#endif
    if(ierr == -1) then
      call io_binary_parallel_start(fname, file_handle, comm, xlocal, np, R_SIZEOF, .false., ierr)
      ASSERT(product(ubound(ff)) >= np)

#ifdef HAVE_MPI2
      call MPI_File_read(file_handle, ff(1), np, R_MPITYPE, status, mpi_err)
      call MPI_Get_count(status, R_MPITYPE, read_count, mpi_err)
      if (read_count /= np) then
        write(message(1),'(1x,2a,i8,a,i8)') TOSTRING(R_TYPE), " read elements=", read_count, " instead of", np
        write(message(2), '(a,a)') " of file= ", fname
        call messages_fatal(2)
      end if
#endif
    endif
    
    call io_binary_parallel_end(file_handle)

    POP_SUB(X(read_parallel))
  end subroutine X(read_parallel)

  !------------------------------------------------------

  subroutine X(read_binary2)(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    R_TYPE,              intent(out) :: ff(:,:)
    integer,             intent(out) :: ierr

    PUSH_SUB(X(read_binary2))

    ASSERT(product(ubound(ff)) >= np)

    call read_binary(np, 0, ff(1,1), R_TYPE_IOBINARY, ierr, trim(fname))

    POP_SUB(X(read_binary2))
  end subroutine X(read_binary2)

  !------------------------------------------------------

  subroutine X(read_binary3)(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer,             intent(in)  :: np
    R_TYPE,              intent(out) :: ff(:,:,:)
    integer,             intent(out) :: ierr

    PUSH_SUB(X(read_binary3))

    ASSERT(product(ubound(ff)) >= np)

    call read_binary(np, 0, ff(1,1,1), R_TYPE_IOBINARY, ierr, trim(fname))

    POP_SUB(X(read_binary3))
  end subroutine X(read_binary3)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
