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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$


! ---------------------------------------------------------
subroutine X(restart_write_function)(dir, filename, gr, f, ierr, size)
  character(len=*), intent(in)  :: dir
  character(len=*), intent(in)  :: filename
  type(grid_t),     intent(in)  :: gr
  integer,          intent(out) :: ierr
  integer,          intent(in)  :: size
  R_TYPE,           intent(in)  :: f(size)

  call push_sub('restart_inc.Xrestart_write_function')

  call X(output_function) (restart_format, trim(dir), trim(filename), &
    gr%m, gr%sb, f(:), M_ONE, ierr, is_tmp=.true.)

  call pop_sub()
end subroutine X(restart_write_function)


! ---------------------------------------------------------
subroutine X(restart_read_function)(dir, filename, m, f, ierr)
  character(len=*), intent(in)  :: dir
  character(len=*), intent(in)  :: filename
  type(mesh_t),     intent(in)  :: m
  R_TYPE,           intent(out) :: f(1:m%np)
  integer,          intent(out) :: ierr

  call push_sub('restart_inc.Xrestart_read_function')

  ! try first to load plain binary files
  call X(input_function) (trim(dir)//'/'//trim(filename), m, f(:), ierr, is_tmp=.true.)

  ! if we do not succeed try NetCDF
  if(ierr>0) call X(input_function) (trim(dir)//'/'//trim(filename)//'.ncdf', m, f(:), ierr, is_tmp=.true.)

  call pop_sub()
end subroutine X(restart_read_function)
