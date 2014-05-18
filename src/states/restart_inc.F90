!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
subroutine X(restart_write_function)(dir, filename, mesh, ff, ierr)
  character(len=*), intent(in)  :: dir
  character(len=*), intent(in)  :: filename
  type(mesh_t),     intent(in)  :: mesh
  R_TYPE,           intent(in)  :: ff(:)
  integer,          intent(out) :: ierr

  PUSH_SUB(X(restart_write_function))

  call X(io_function_output)(restart_format, trim(dir), trim(filename), mesh, ff(:), unit_one, ierr, is_tmp=.true.)
  ! all restart files are in atomic units

  POP_SUB(X(restart_write_function))
end subroutine X(restart_write_function)


! ---------------------------------------------------------
!> In domain parallel case each process reads a part of the file.
!! At the end all the processes have the corresponding mesh part
subroutine X(restart_read_function)(dir, filename, mesh, ff, ierr, map)
  character(len=*), intent(in)    :: dir
  character(len=*), intent(in)    :: filename
  type(mesh_t),     intent(in)    :: mesh
  R_TYPE, target,   intent(inout) :: ff(:)
  integer,          intent(out)   :: ierr
  integer, optional, intent(in)   :: map(:)

  integer :: ip, np, offset
  R_TYPE, pointer :: read_ff(:)
  type(profile_t), save :: prof_io
  type(batch_t) :: ffb
  type(profile_t), save :: prof_comm

  PUSH_SUB(X(restart_read_function))
  
  nullify(read_ff)

  if(present(map) .and. mesh%parallel_in_domains) then 
    ! for the moment we do not do this directly
    call X(io_function_input) (trim(dir)//'/'//trim(filename)//'.obf', mesh, ff(1:mesh%np), ierr, is_tmp=.true., map = map)

    POP_SUB(X(restart_read_function))
    return
  end if

  if(present(map)) then
    call io_binary_get_info(trim(dir)//'/'//trim(filename)//'.obf', np, ierr)

    if(ierr /= 0) then
      POP_SUB(X(restart_read_function))
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

  call io_binary_read(trim(dir)//'/'//trim(filename)//'.obf', np, read_ff, ierr, offset = offset)
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

  if(present(map)) then
    ff(1:mesh%np_global) = M_ZERO
    do ip = 1, min(np, ubound(map, dim = 1))
      if(map(ip) > 0) ff(map(ip)) = read_ff(ip)
    end do
    
    SAFE_DEALLOCATE_P(read_ff)
  end if

  POP_SUB(X(restart_read_function))
end subroutine X(restart_read_function)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
