!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2014 M. Oliveira
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


! ---------------------------------------------------------
subroutine X(restart_write_mesh_function)(restart, filename, mesh, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  type(mesh_t),     intent(in)  :: mesh
  R_TYPE,           intent(in)  :: ff(:)
  integer,          intent(out) :: ierr

  PUSH_SUB(X(restart_write_mesh_function))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_DUMP)

  call X(io_function_output)(restart%format, trim(restart%pwd), trim(filename), mesh, ff(:), unit_one, ierr, is_tmp=.true.)
  ! all restart files are in atomic units

  if (ierr /= 0) then
    message(1) = "Unable to write restart function to '"//trim(restart%pwd)//"/"//trim(filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_write_mesh_function))
end subroutine X(restart_write_mesh_function)


! ---------------------------------------------------------
!> In domain parallel case each process reads a part of the file.
!! At the end all the processes have the corresponding mesh part
subroutine X(restart_read_mesh_function)(restart, filename, mesh, ff, ierr)
  type(restart_t),  intent(in)    :: restart
  character(len=*), intent(in)    :: filename
  type(mesh_t),     intent(in)    :: mesh
  R_TYPE, target,   intent(inout) :: ff(:)
  integer,          intent(out)   :: ierr

  integer :: ip, np, offset
  R_TYPE, pointer :: read_ff(:)
  type(profile_t), save :: prof_io
  type(batch_t) :: ffb
  type(profile_t), save :: prof_comm

  PUSH_SUB(X(restart_read_mesh_function))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_LOAD)  

  nullify(read_ff)

  if (associated(restart%map) .and. mesh%parallel_in_domains) then 
    ! for the moment we do not do this directly
    call X(io_function_input) (trim(restart%pwd)//'/'//trim(filename)//'.obf', mesh, ff(1:mesh%np), ierr, &
                               is_tmp=.true., map = restart%map)

    POP_SUB(X(restart_read_mesh_function))
    return
  end if

  if (associated(restart%map)) then
    call io_binary_get_info(trim(restart%pwd)//'/'//trim(filename)//'.obf', np, ierr)

    if (ierr /= 0) then
      POP_SUB(X(restart_read_mesh_function))
      return
    end if

    ASSERT(np > 0)
    SAFE_ALLOCATE(read_ff(1:np))
  else
    np = mesh%np
    read_ff => ff
  end if

  offset = 0
  !in the parallel case, each node reads a part of the file
  if(mesh%parallel_in_domains) then
    offset = mesh%vp%xlocal - 1
  end if

  ASSERT(associated(read_ff))

  call profiling_in(prof_io, "RESTART_READ_IO")

#ifdef HAVE_MPI2
  ! Ensure that xlocal has a proper value
  ASSERT(mesh%vp%xlocal >= 0 .and. mesh%vp%xlocal <= mesh%np_part_global)
  call io_binary_read_parallel(trim(restart%pwd)//'/'//trim(filename)//'.obf', mesh%mpi_grp%comm, mesh%vp%xlocal, &
                               np, read_ff, ierr)
#else
  call io_binary_read(trim(restart%pwd)//'/'//trim(filename)//'.obf', np, read_ff, ierr, offset = offset)
#endif
  call profiling_count_transfers(np, read_ff(1))
  call profiling_out(prof_io)

  if(mesh%parallel_in_domains) then
    call profiling_in(prof_comm, "RESTART_READ_COMM")
    ! this is the global index of the points we read

    ff(1:mesh%np) = read_ff(1:mesh%np)

    call batch_init(ffb, 1)
    call batch_add_state(ffb, ff)
    call X(mesh_batch_exchange_points)(mesh, ffb, backward_map = .true.)
    call batch_end(ffb)
    
    call profiling_out(prof_comm)
  end if

  if (associated(restart%map)) then
    ff(1:mesh%np_global) = M_ZERO
    do ip = 1, min(np, ubound(restart%map, dim = 1))
      if (restart%map(ip) > 0) ff(restart%map(ip)) = read_ff(ip)
    end do
    
    SAFE_DEALLOCATE_P(read_ff)
  end if

  if (ierr /= 0) then
    message(1) = "Unable to read mesh function from '"//trim(restart%pwd)//"/"//trim(filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_read_mesh_function))
end subroutine X(restart_read_mesh_function)


! ---------------------------------------------------------
subroutine X(restart_write_binary1)(restart, filename, np, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: np
  R_TYPE,           intent(in)  :: ff(:)
  integer,          intent(out) :: ierr

  PUSH_SUB(X(restart_write_binary1))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_DUMP)

  !Only the root node writes
  if (mpi_grp_is_root(restart%mpi_grp)) then
    call io_binary_write(trim(restart%pwd)//"/"//trim(filename)//".obf", np, ff, ierr)
  end if

#if defined(HAVE_MPI)
  call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, restart%mpi_grp%comm, mpi_err)
  call MPI_Barrier(restart%mpi_grp%comm, mpi_err)
#endif

  if (ierr /= 0) then
    message(1) = "Unable to write restart information to '"//trim(restart%pwd)//"/"//trim(filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_write_binary1))
end subroutine X(restart_write_binary1)

! ---------------------------------------------------------
subroutine X(restart_write_binary3)(restart, filename, np, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: np
  R_TYPE,           intent(in)  :: ff(:,:,:)
  integer,          intent(out) :: ierr

  PUSH_SUB(X(restart_write_binary3))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_DUMP)

  !Only the root node writes
  if (mpi_grp_is_root(restart%mpi_grp)) then
    call io_binary_write(trim(restart%pwd)//"/"//trim(filename)//".obf", np, ff, ierr)
  end if

#if defined(HAVE_MPI)
  call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, restart%mpi_grp%comm, mpi_err)
  call MPI_Barrier(restart%mpi_grp%comm, mpi_err)
#endif

  if (ierr /= 0) then
    message(1) = "Unable to write restart information to '"//trim(restart%pwd)//"/"//trim(filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_write_binary3))
end subroutine X(restart_write_binary3)


! ---------------------------------------------------------
subroutine X(restart_read_binary1)(restart, filename, np, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: np
  R_TYPE,           intent(out) :: ff(:)
  integer,          intent(out) :: ierr

  PUSH_SUB(X(restart_read_binary1))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_LOAD)

  call io_binary_read(trim(restart%pwd)//"/"//trim(filename)//".obf", np, ff, ierr)

  if (ierr /= 0) then
    message(1) = "Unable to read restart information from '"//trim(restart%pwd)//"/"//trim(filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_read_binary1))
end subroutine X(restart_read_binary1)


! ---------------------------------------------------------
subroutine X(restart_read_binary3)(restart, filename, np, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: np
  R_TYPE,           intent(out) :: ff(:,:,:)
  integer,          intent(out) :: ierr

  PUSH_SUB(X(restart_read_binary3))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_LOAD)

  call io_binary_read(trim(restart%pwd)//"/"//trim(filename)//".obf", np, ff, ierr)

  if (ierr /= 0) then
    message(1) = "Unable to read restart information from '"//trim(restart%pwd)//"/"//trim(filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_read_binary3))
end subroutine X(restart_read_binary3)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
