!! Copyright (C) 2010 X. Andrade
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
!! $Id: opencl_inc.F90 3587 2007-11-22 16:43:00Z xavier $


subroutine X(opencl_write_buffer)(this, size, data, offset)
  type(opencl_mem_t),               intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(in)    :: data(:)
  integer(SIZEOF_SIZE_T), optional, intent(in)    :: offset

  integer(SIZEOF_SIZE_T) :: fsize, offset_
  
  fsize = size*R_SIZEOF
  offset_ = 0
  if(present(offset)) offset_ = offset*R_SIZEOF

  call f90_opencl_write_buffer(this%mem, opencl%env, fsize, offset_, data(1))

end subroutine X(opencl_write_buffer)

! -----------------------------------------------------------------------------

subroutine X(opencl_read_buffer)(this, size, data, offset)
  type(opencl_mem_t),               intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(out)   :: data(:)
  integer(SIZEOF_SIZE_T), optional, intent(in)    :: offset

  integer(SIZEOF_SIZE_T) :: fsize, offset_
  
  fsize = size*R_SIZEOF
  offset_ = 0
  if(present(offset)) offset_ = offset*R_SIZEOF

  call f90_opencl_read_buffer(this%mem, opencl%env, fsize, offset_, data(1))

end subroutine X(opencl_read_buffer)

! ---------------------------------------------------------------------------

subroutine X(opencl_set_kernel_arg_data)(kernel, type, narg, data)
  type(c_ptr),        intent(inout) :: kernel
  integer,            intent(in)    :: type
  integer,            intent(in)    :: narg
  R_TYPE,             intent(in)    :: data
  
  call f90_opencl_set_kernel_arg_data(kernel, type, narg, R_SIZEOF, data)
  
end subroutine X(opencl_set_kernel_arg_data)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
