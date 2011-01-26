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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$


! ---------------------------------------------------------
subroutine X(restart_write_function)(dir, filename, gr, ff, ierr, size)
  character(len=*), intent(in)  :: dir
  character(len=*), intent(in)  :: filename
  type(grid_t),     intent(in)  :: gr
  integer,          intent(out) :: ierr
  integer,          intent(in)  :: size
  R_TYPE,           intent(in)  :: ff(size)

  PUSH_SUB(X(restart_write_function))

  call X(output_function)(restart_format, trim(dir), trim(filename), gr%mesh, ff(:), unit_one, ierr, is_tmp=.true.)
  ! all restart files are in atomic units

  POP_SUB(X(restart_write_function))
end subroutine X(restart_write_function)


! ---------------------------------------------------------
subroutine X(restart_read_function)(dir, filename, mesh, ff, ierr, map)
  character(len=*), intent(in)    :: dir
  character(len=*), intent(in)    :: filename
  type(mesh_t),     intent(in)    :: mesh
  R_TYPE, target,   intent(inout) :: ff(1:mesh%np)
  integer,          intent(out)   :: ierr
  integer, optional, intent(in)   :: map(:)

  integer :: ip, np, il, ipart, offset
  integer, allocatable :: list(:)
  R_TYPE, pointer :: read_ff(:)

  PUSH_SUB(X(restart_read_function))

#ifdef HAVE_MPI
  if(present(map) .and. mesh%parallel_in_domains) then 
    ! for the moment we do not do this directly
    call X(input_function) (trim(dir)//'/'//trim(filename)//'.obf', mesh, ff(1:mesh%np), ierr, is_tmp=.true., map = map)

    POP_SUB(X(restart_read_function))
    return
  end if
#endif

  if(present(map)) then
    call io_binary_get_info(trim(dir)//'/'//trim(filename)//'.obf', np, ierr)
    SAFE_ALLOCATE(read_ff(1:np))
    offset = 0
  else
    np = mesh%np
    read_ff => ff
    offset = 0
  end if

#ifdef HAVE_MPI
  !in the parallel case, each node reads a part of the file
  if(mesh%parallel_in_domains) then
    SAFE_ALLOCATE(read_ff(1:np))
    offset = mesh%vp%xlocal(mesh%vp%partno) - 1
  end if
#endif

  call io_binary_read(trim(dir)//'/'//trim(filename)//'.obf', np, read_ff, ierr, offset = offset)
  call profiling_count_transfers(np, read_ff(1))

#ifdef HAVE_MPI
  if(mesh%parallel_in_domains) then
    ! no we need to move around the points that we read
    ! this is probably not the most efficient way of doing it
    do ipart = 1, mesh%vp%npart
      SAFE_ALLOCATE(list(1:mesh%vp%np_local(ipart)))

      do il = 1, mesh%vp%np_local(ipart)
        ! this is the global index of the point we read
        list(il) = il + mesh%vp%xlocal(ipart) - 1
      end do

      call X(vec_selective_scatter)(mesh%vp, mesh%vp%np_local(ipart), list, ipart - 1, read_ff, ff)      
      
      SAFE_DEALLOCATE_A(list)
    end do

    SAFE_DEALLOCATE_P(read_ff)
  end if
#endif

  if(present(map)) then
    ff(1:mesh%np_global) = M_ZERO
    do ip = 1, min(np, ubound(map, dim = 1))
      if(map(ip) > 0) ff(map(ip)) = read_ff(ip)
    end do
    
    SAFE_DEALLOCATE_P(read_ff)
  end if

  POP_SUB(X(restart_read_function))
end subroutine X(restart_read_function)


! ---------------------------------------------------------
subroutine X(restart_write_lr_rho)(lr, gr, nspin, restart_dir, rho_tag)
  type(lr_t),        intent(in)    :: lr
  type(grid_t),      intent(in)    :: gr
  integer,           intent(in)    :: nspin
  character(len=*),  intent(in)    :: restart_dir
  character(len=*),  intent(in)    :: rho_tag

  character(len=100) :: fname
  integer :: is, ierr

  PUSH_SUB(X(restart_write_lr_rho))

  call block_signals()
  do is = 1, nspin
    write(fname, '(a,i1,a)') trim(rho_tag)//'_', is
    call X(restart_write_function)(trim(tmpdir)//trim(RESTART_DIR), fname, gr,&
         lr%X(dl_rho)(:, is), ierr, size(lr%X(dl_rho), 1))
  end do
  call unblock_signals()

  POP_SUB(X(restart_write_lr_rho))
end subroutine X(restart_write_lr_rho)


! ---------------------------------------------------------
subroutine X(restart_read_lr_rho)(lr, gr, nspin, restart_subdir, rho_tag, ierr)
  type(lr_t),        intent(inout) :: lr
  type(grid_t),      intent(in)    :: gr
  integer,           intent(in)    :: nspin
  character(len=*),  intent(in)    :: restart_subdir
  character(len=*),  intent(in)    :: rho_tag
  integer,           intent(out)   :: ierr

  character(len=80) :: fname
  integer :: is, s_ierr

  PUSH_SUB(X(restart_read_lr_rho))

  ierr = 0
  do is = 1, nspin
    write(fname, '(a, i1,a)') trim(rho_tag)//'_', is
    call X(restart_read_function)(trim(restart_dir)//trim(restart_subdir), fname, gr%mesh,&
         lr%X(dl_rho)(:, is), s_ierr)
    if( s_ierr /=0 ) ierr = s_ierr
  end do


  if( ierr == 0 ) then 
    write(message(1),'(a)') 'Loaded restart density '//rho_tag
    call write_info(1)

  else

    write(message(1),'(a)') 'Could not load restart '//rho_tag
    call write_info(1)

  end if

  POP_SUB(X(restart_read_lr_rho))
end subroutine X(restart_read_lr_rho)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
