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

! ---------------------------------------------------------
subroutine X(restart_write_mesh_function)(restart, space, filename, mesh, ff, ierr, root)
  type(restart_t),   intent(in)  :: restart
  type(space_t),     intent(in)  :: space
  character(len=*),  intent(in)  :: filename
  type(mesh_t),      intent(in)  :: mesh
  R_TYPE,  target,   intent(in)  :: ff(:)
  integer,           intent(out) :: ierr
  integer, optional, intent(in)  :: root(:) !< which process is going to write the data

  logical :: i_am_root, in_line
  integer :: root_(1:P_STRATEGY_MAX)
  R_TYPE, pointer :: ff_global(:)
  
  PUSH_SUB(X(restart_write_mesh_function))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_DUMP)
  ASSERT(restart%has_mesh)  

  root_(1:P_STRATEGY_MAX) = 0
  if (present(root)) then
    ASSERT(root(P_STRATEGY_DOMAINS) >= 0)
    where(root >= 0)
      root_ = root
    elsewhere
      root_ = restart%mc%who_am_i
    end where
  end if
    
  ierr = 0
  i_am_root = all(root_ == restart%mc%who_am_i)
  in_line = (root_(P_STRATEGY_STATES)  == restart%mc%who_am_i(P_STRATEGY_STATES)  .and. &
             root_(P_STRATEGY_KPOINTS) == restart%mc%who_am_i(P_STRATEGY_KPOINTS) .and. &
             root_(P_STRATEGY_OTHER)   == restart%mc%who_am_i(P_STRATEGY_OTHER)  )

  if (i_am_root) then
    if (mesh%parallel_in_domains) then
      SAFE_ALLOCATE(ff_global(1:mesh%np_global))
    else
      ff_global => ff
    end if
  end if

  if (in_line .and. mesh%parallel_in_domains) then
    if (i_am_root) then
      call vec_gather(mesh%vp, root_(P_STRATEGY_DOMAINS), ff, ff_global)
    else
      call vec_gather(mesh%vp, root_(P_STRATEGY_DOMAINS), ff)
    end if
  end if
  
  if (i_am_root) then
    call X(io_function_output_global)(restart%format, restart%pwd, filename, restart%namespace, &
      space, mesh, ff_global, unit_one, ierr)
    ! all restart files are in atomic units

    if (mesh%parallel_in_domains) then
      SAFE_DEALLOCATE_P(ff_global)
    else
      nullify(ff_global)
    end if

    if (ierr /= 0) then
      message(1) = "Unable to write restart function to '"//trim(restart%pwd)//"/"//trim(filename)//"'."
      call messages_warning(1)
    end if
  end if
  
  if (mesh%parallel_in_domains .and. in_line) then
#ifdef HAVE_MPI
    ! I have to broadcast the error code
    call MPI_Bcast(ierr, 1, MPI_INTEGER, root_(P_STRATEGY_DOMAINS), mesh%mpi_grp%comm, mpi_err)
#endif
  end if
    
  POP_SUB(X(restart_write_mesh_function))
end subroutine X(restart_write_mesh_function)

! ---------------------------------------------------------
!> In domain parallel case each process reads a part of the file.
!! At the end all the processes have the corresponding mesh part
subroutine X(restart_read_mesh_function)(restart, space, filename, mesh, ff, ierr)
  type(restart_t),            intent(in)    :: restart
  type(space_t),              intent(in)    :: space
  character(len=*),           intent(in)    :: filename
  type(mesh_t),               intent(in)    :: mesh
  R_TYPE, target, contiguous, intent(inout) :: ff(:)
  integer,                    intent(out)   :: ierr

  integer :: ip, np, offset, file_size
  R_TYPE, pointer :: read_ff(:)
  type(profile_t), save :: prof_io
  type(batch_t) :: ffb
  type(profile_t), save :: prof_comm
  character(len=MAX_PATH_LEN) :: full_filename
  
  PUSH_SUB(X(restart_read_mesh_function))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_LOAD)
  ASSERT(restart%has_mesh)

  nullify(read_ff)
  full_filename = trim(restart%pwd)//'/'//trim(filename)//'.obf'
  
  if (restart_has_map(restart) .and. mesh%parallel_in_domains) then 
    ! for the moment we do not do this directly
    call X(io_function_input) (full_filename, restart%namespace, space, mesh, ff(1:mesh%np), ierr, &
                               map = restart%map)

    POP_SUB(X(restart_read_mesh_function))
    return
  end if

  if (restart_has_map(restart)) then
    call io_binary_get_info(io_workpath(full_filename, restart%namespace), np, file_size, ierr)

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

  call profiling_in(prof_io, TOSTRING(X(RESTART_READ_IO)))

#ifdef HAVE_MPI2
  if(mesh%parallel_in_domains) then
    ! Ensure that xlocal has a proper value
    ASSERT(mesh%vp%xlocal >= 0 .and. mesh%vp%xlocal <= mesh%np_part_global)
    call io_binary_read_parallel(io_workpath(full_filename, restart%namespace), &
      mesh%mpi_grp%comm, mesh%vp%xlocal, np, read_ff, ierr)
  else
#endif
    call io_binary_read(io_workpath(full_filename, restart%namespace), np, &
      read_ff, ierr, offset = offset)
#ifdef HAVE_MPI2
  end if
#endif
  call profiling_count_transfers(np, read_ff(1))
  call profiling_out(prof_io)

  if(mesh%parallel_in_domains) then
    call profiling_in(prof_comm, TOSTRING(X(RESTART_READ_COMM)))
    ! this is the global index of the points we read

    ff(1:mesh%np) = read_ff(1:mesh%np)

    call batch_init(ffb, ff)
    call X(mesh_batch_exchange_points)(mesh, ffb, backward_map = .true.)
    call ffb%end()
    
    call profiling_out(prof_comm)
  end if

  if (restart_has_map(restart)) then
    ff(1:mesh%np_global) = M_ZERO
    do ip = 1, min(np, ubound(restart%map, dim = 1))
      if (restart%map(ip) > 0) ff(restart%map(ip)) = read_ff(ip)
    end do
    
    SAFE_DEALLOCATE_P(read_ff)
  end if

  if (ierr /= 0) then
    message(1) = "Unable to read mesh function from '"//&
      trim(io_workpath(full_filename, restart%namespace))//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_read_mesh_function))
end subroutine X(restart_read_mesh_function)


! ---------------------------------------------------------
subroutine X(restart_write_binary1)(restart, filename, np, ff, ierr, root)
  type(restart_t),   intent(in)  :: restart
  character(len=*),  intent(in)  :: filename
  integer,           intent(in)  :: np
  R_TYPE,            intent(in)  :: ff(:)
  integer,           intent(out) :: ierr
  integer, optional, intent(in)  :: root(:) !< which process is going to write the data

  character(len=MAX_PATH_LEN) :: full_filename
  integer :: root_(1:P_STRATEGY_MAX)
    
  PUSH_SUB(X(restart_write_binary1))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_DUMP)

  full_filename = trim(io_workpath(restart%pwd, restart%namespace))//"/"//trim(filename)//".obf"

  root_(1:P_STRATEGY_MAX) = 0
  if (present(root)) then
    ASSERT(root(P_STRATEGY_DOMAINS) >= 0)
    where(root >= 0)
      root_ = root
    elsewhere
      root_ = restart%mc%who_am_i
    end where
  end if
  
  ierr = 0
  !Only the root node writes
  if (all(root_ == restart%mc%who_am_i)) then
    call io_binary_write(full_filename, np, ff, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write restart information to '"//trim(full_filename)//"'."
      call messages_warning(1, all_nodes=.true.)
    end if
  end if

  POP_SUB(X(restart_write_binary1))
end subroutine X(restart_write_binary1)

! ---------------------------------------------------------
subroutine X(restart_write_binary2)(restart, filename, np, ff, ierr, root)
  type(restart_t),   intent(in)  :: restart
  character(len=*),  intent(in)  :: filename
  integer,           intent(in)  :: np
  R_TYPE,            intent(in)  :: ff(:,:)
  integer,           intent(out) :: ierr
  integer, optional, intent(in)  :: root(:) !< which process is going to write the data

  character(len=MAX_PATH_LEN) :: full_filename
  integer :: root_(1:P_STRATEGY_MAX)

  PUSH_SUB(X(restart_write_binary2))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_DUMP)

  full_filename = trim(io_workpath(restart%pwd, restart%namespace))//"/"//trim(filename)//".obf"

  root_(1:P_STRATEGY_MAX) = 0
  if (present(root)) then
    ASSERT(root(P_STRATEGY_DOMAINS) >= 0)
    where(root >= 0)
      root_ = root
    elsewhere
      root_ = restart%mc%who_am_i
    end where
  end if
  
  ierr = 0
  !Only the root node writes
  if (all(root_ == restart%mc%who_am_i)) then
    call io_binary_write(full_filename, np, ff, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write restart information to '"//trim(full_filename)//"'."
      call messages_warning(1, all_nodes=.true.)
    end if
  end if

  POP_SUB(X(restart_write_binary2))
end subroutine X(restart_write_binary2)

! ---------------------------------------------------------
subroutine X(restart_write_binary3)(restart, filename, np, ff, ierr, root)
  type(restart_t),   intent(in)  :: restart
  character(len=*),  intent(in)  :: filename
  integer,           intent(in)  :: np
  R_TYPE,            intent(in)  :: ff(:,:,:)
  integer,           intent(out) :: ierr
  integer, optional, intent(in)  :: root(:) !< which process is going to write the data
  
  character(len=MAX_PATH_LEN) :: full_filename
  integer :: root_(1:P_STRATEGY_MAX)

  PUSH_SUB(X(restart_write_binary3))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_DUMP)

  full_filename = trim(io_workpath(restart%pwd, restart%namespace))//"/"//trim(filename)//".obf"

  root_(1:P_STRATEGY_MAX) = 0
  if (present(root)) then
    ASSERT(root(P_STRATEGY_DOMAINS) >= 0)
    where(root >= 0)
      root_ = root
    elsewhere
      root_ = restart%mc%who_am_i
    end where
  end if
  
  ierr = 0
  !Only the root node writes
  if (all(root_ == restart%mc%who_am_i)) then
    call io_binary_write(full_filename, np, ff, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write restart information to '"//trim(full_filename)//"'."
      call messages_warning(1, all_nodes=.true.)
    end if
  end if

  POP_SUB(X(restart_write_binary3))
end subroutine X(restart_write_binary3)

! ---------------------------------------------------------
subroutine X(restart_write_binary5)(restart, filename, np, ff, ierr, root)
  type(restart_t),   intent(in)  :: restart
  character(len=*),  intent(in)  :: filename
  integer,           intent(in)  :: np
  R_TYPE,            intent(in)  :: ff(:,:,:,:,:)
  integer,           intent(out) :: ierr
  integer, optional, intent(in)  :: root(:) !< which process is going to write the data

  character(len=MAX_PATH_LEN) :: full_filename
  integer :: root_(1:P_STRATEGY_MAX)

  PUSH_SUB(X(restart_write_binary5))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_DUMP)

  full_filename = trim(io_workpath(restart%pwd, restart%namespace))//"/"//trim(filename)//".obf"

  root_(1:P_STRATEGY_MAX) = 0
  if (present(root)) then
    ASSERT(root(P_STRATEGY_DOMAINS) >= 0)
    where(root >= 0)
      root_ = root
    elsewhere
      root_ = restart%mc%who_am_i
    end where
  end if
  
  ierr = 0
  !Only the root node writes
  if (all(root_ == restart%mc%who_am_i)) then
    call io_binary_write(full_filename, np, ff, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write restart information to '"//trim(full_filename)//"'."
      call messages_warning(1, all_nodes=.true.)
    end if
  end if

  POP_SUB(X(restart_write_binary5))
end subroutine X(restart_write_binary5)

! ---------------------------------------------------------
subroutine X(restart_read_binary1)(restart, filename, np, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: np
  R_TYPE,           intent(out) :: ff(:)
  integer,          intent(out) :: ierr

  character(len=MAX_PATH_LEN) :: full_filename

  PUSH_SUB(X(restart_read_binary1))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_LOAD)

  full_filename = trim(io_workpath(restart%pwd, restart%namespace))//"/"//trim(filename)//".obf"

  call io_binary_read(full_filename, np, ff, ierr)

  if (ierr /= 0) then
    message(1) = "Unable to read restart information from '"//trim(full_filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_read_binary1))
end subroutine X(restart_read_binary1)


! ---------------------------------------------------------
subroutine X(restart_read_binary2)(restart, filename, np, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: np
  R_TYPE,           intent(out) :: ff(:,:)
  integer,          intent(out) :: ierr

  character(len=MAX_PATH_LEN) :: full_filename

  PUSH_SUB(X(restart_read_binary2))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_LOAD)

  full_filename = trim(io_workpath(restart%pwd, restart%namespace))//"/"//trim(filename)//".obf"

  call io_binary_read(full_filename, np, ff, ierr)

  if (ierr /= 0) then
    message(1) = "Unable to read restart information from '"//trim(full_filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_read_binary2))
end subroutine X(restart_read_binary2)


! ---------------------------------------------------------
subroutine X(restart_read_binary3)(restart, filename, np, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: np
  R_TYPE,           intent(out) :: ff(:,:,:)
  integer,          intent(out) :: ierr

  character(len=MAX_PATH_LEN) :: full_filename

  PUSH_SUB(X(restart_read_binary3))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_LOAD)

  full_filename = trim(io_workpath(restart%pwd, restart%namespace))//"/"//trim(filename)//".obf"

  call io_binary_read(full_filename, np, ff, ierr)

  if (ierr /= 0) then
    message(1) = "Unable to read restart information from '"//trim(full_filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_read_binary3))
end subroutine X(restart_read_binary3)

! ---------------------------------------------------------
subroutine X(restart_read_binary5)(restart, filename, np, ff, ierr)
  type(restart_t),  intent(in)  :: restart
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: np
  R_TYPE,           intent(out) :: ff(:,:,:,:,:)
  integer,          intent(out) :: ierr

  character(len=MAX_PATH_LEN) :: full_filename

  PUSH_SUB(X(restart_read_binary5))

  ASSERT(.not. restart%skip)
  ASSERT(restart%type == RESTART_TYPE_LOAD)

  full_filename = trim(io_workpath(restart%pwd, restart%namespace))//"/"//trim(filename)//".obf"

  call io_binary_read(full_filename, np, ff, ierr)

  if (ierr /= 0) then
    message(1) = "Unable to read restart information from '"//trim(full_filename)//"'."
    call messages_warning(1)
  end if

  POP_SUB(X(restart_read_binary5))
end subroutine X(restart_read_binary5)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
