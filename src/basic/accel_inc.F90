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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

subroutine X(accel_write_buffer_0)(this, data, async)
  type(accel_mem_t),               intent(inout) :: this
  R_TYPE,                          intent(in)    :: data
  logical,               optional, intent(in)    :: async

  R_TYPE, allocatable :: data_vec(:)
  
  PUSH_SUB(X(accel_write_buffer_0))

  SAFE_ALLOCATE(data_vec(1:1))

  data_vec(1:1) = data

  call X(accel_write_buffer_1)(this, 1, data_vec, async=async)
  
  SAFE_DEALLOCATE_A(data_vec)
  
  POP_SUB(X(accel_write_buffer_0))
end subroutine X(accel_write_buffer_0)

! -----------------------------------------------------------------------------

subroutine X(accel_write_buffer_1)(this, size, data, offset, async)
  type(accel_mem_t),               intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(in)    :: data(:)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  integer(8) :: fsize, offset_
  logical :: async_
#ifdef HAVE_OPENCL
  integer :: ierr
#endif

  PUSH_SUB(X(accel_write_buffer_1))
  call profiling_in(prof_write, TOSTRING(X(CL_WRITE_BUFFER)))

  ! it does not make sense to write a buffer that the kernels cannot read
  ASSERT(this%flags /= ACCEL_MEM_WRITE_ONLY)

  fsize = int(size,8)*R_SIZEOF
  offset_ = 0
  if(present(offset)) offset_ = int(offset, 8)*R_SIZEOF

  async_ = optional_default(async, .false.)

  if(fsize > 0) then
    
#ifdef HAVE_OPENCL
    call clEnqueueWriteBuffer(accel%command_queue, this%mem, cl_bool(.true.), offset_, fsize, data(1), ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueWriteBuffer")
#endif
#ifdef HAVE_CUDA
    call cuda_memcpy_htod(this%mem, data(1), fsize, offset_, async_)
#endif
    
    call profiling_count_transfers(size, data(1))
    if (.not. async_) call accel_finish()

  end if
  
  call profiling_out(prof_write)
  POP_SUB(X(accel_write_buffer_1))

end subroutine X(accel_write_buffer_1)

! -----------------------------------------------------------------------------

subroutine X(accel_write_buffer_2)(this, size, data, offset, async)
  type(accel_mem_t),               intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(in)    :: data(:, :)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  integer(8) :: fsize, offset_
  logical :: async_
#ifdef HAVE_OPENCL
  integer :: ierr
#endif

  PUSH_SUB(X(accel_write_buffer_2))
  call profiling_in(prof_write, TOSTRING(X(CL_WRITE_BUFFER)))

  ! it does not make sense to write a buffer that the kernels cannot read
  ASSERT(this%flags /= ACCEL_MEM_WRITE_ONLY)

  fsize = int(size, 8)*R_SIZEOF
  offset_ = 0
  if(present(offset)) offset_ = int(offset, 8)*R_SIZEOF

  async_ = optional_default(async, .false.)

  if(fsize > 0) then
    
#ifdef HAVE_OPENCL
    call clEnqueueWriteBuffer(accel%command_queue, this%mem, cl_bool(.true.), offset_, fsize, data(1, 1), ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueWriteBuffer")
#endif
#ifdef HAVE_CUDA
    call cuda_memcpy_htod(this%mem, data(1, 1), fsize, offset_, async_)
#endif
    
    call profiling_count_transfers(size, data(1, 1))
    if (.not. async_) call accel_finish()
    
  end if
    
  call profiling_out(prof_write)
  POP_SUB(X(accel_write_buffer_2))

end subroutine X(accel_write_buffer_2)

! -----------------------------------------------------------------------------

subroutine X(accel_write_buffer_3)(this, size, data, offset, async)
  type(accel_mem_t),               intent(inout) :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(in)    :: data(:, :, :)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  integer(8) :: fsize, offset_
  logical :: async_
#ifdef HAVE_OPENCL
  integer :: ierr
#endif

  PUSH_SUB(X(accel_write_buffer_3))
  call profiling_in(prof_write, TOSTRING(X(CL_WRITE_BUFFER)))

  ! it does not make sense to write a buffer that the kernels cannot read
  ASSERT(this%flags /= ACCEL_MEM_WRITE_ONLY)

  fsize = int(size, 8)*R_SIZEOF
  offset_ = 0
  if(present(offset)) offset_ = int(offset, 8)*R_SIZEOF

  async_ = optional_default(async, .false.)

  if(fsize > 0) then
  
#ifdef HAVE_OPENCL
    call clEnqueueWriteBuffer(accel%command_queue, this%mem, cl_bool(.true.), offset_, fsize, data(1, 1, 1), ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueWriteBuffer")
#endif
#ifdef HAVE_CUDA
    call cuda_memcpy_htod(this%mem, data(1, 1, 1), fsize, offset_, async_)
#endif
    
    call profiling_count_transfers(size, data(1, 1, 1))
    if (.not. async_) call accel_finish()

  end if
  
  call profiling_out(prof_write)
  POP_SUB(X(accel_write_buffer_3))

end subroutine X(accel_write_buffer_3)

! -----------------------------------------------------------------------------

subroutine X(accel_read_buffer_1)(this, size, data, offset, async)
  type(accel_mem_t),               intent(in)    :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(out)   :: data(:)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  integer(8) :: fsize, offset_
  logical :: async_
#ifdef HAVE_OPENCL
  integer :: ierr
#endif

  PUSH_SUB(X(accel_read_buffer_1))
  call profiling_in(prof_read, TOSTRING(X(CL_READ_BUFFER)))

  ! it does not make sense to read a buffer that the kernels cannot write
  ASSERT(this%flags /= ACCEL_MEM_READ_ONLY)

  fsize = size*R_SIZEOF
  offset_ = 0
  if(present(offset)) offset_ = offset*R_SIZEOF

  async_ = optional_default(async, .false.)

  if(fsize > 0) then
    
#ifdef HAVE_OPENCL
    call clEnqueueReadBuffer(accel%command_queue, this%mem, cl_bool(.true.), offset_, fsize, data(1), ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueReadBuffer")
#endif
#ifdef HAVE_CUDA
    call cuda_memcpy_dtoh(this%mem, data(1), fsize, offset_, async_)
#endif
    
    call profiling_count_transfers(size, data(1))
    if (.not. async_) call accel_finish()

  end if
    
  call profiling_out(prof_read)
  POP_SUB(X(accel_read_buffer_1))

end subroutine X(accel_read_buffer_1)

! ---------------------------------------------------------------------------

subroutine X(accel_read_buffer_2)(this, size, data, offset, async)
  type(accel_mem_t),               intent(in)    :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(out)   :: data(:, :)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  integer(8) :: fsize, offset_
  logical :: async_
#ifdef HAVE_OPENCL
  integer :: ierr
#endif
  
  PUSH_SUB(X(accel_read_buffer_2))
  call profiling_in(prof_read, TOSTRING(X(CL_READ_BUFFER)))

  ! it does not make sense to read a buffer that the kernels cannot write
  ASSERT(this%flags /= ACCEL_MEM_READ_ONLY)

  fsize = size*R_SIZEOF
  offset_ = 0
  if(present(offset)) offset_ = offset*R_SIZEOF

  async_ = optional_default(async, .false.)

  if(fsize > 0) then
    
#ifdef HAVE_OPENCL
    call clEnqueueReadBuffer(accel%command_queue, this%mem, cl_bool(.true.), offset_, fsize, data(1, 1), ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueReadBuffer")
#endif
#ifdef HAVE_CUDA
    call cuda_memcpy_dtoh(this%mem, data(1, 1), fsize, offset_, async_)
#endif
    
    call profiling_count_transfers(size, data(1, 1))
    if (.not. async_) call accel_finish()

  end if
  
  call profiling_out(prof_read)
  POP_SUB(X(accel_read_buffer_2))
  
end subroutine X(accel_read_buffer_2)

! ---------------------------------------------------------------------------

subroutine X(accel_read_buffer_3)(this, size, data, offset, async)
  type(accel_mem_t),               intent(in)    :: this
  integer,                          intent(in)    :: size
  R_TYPE,                           intent(out)   :: data(:, :, :)
  integer,                optional, intent(in)    :: offset
  logical,                optional, intent(in)    :: async

  integer(8) :: fsize, offset_
  logical :: async_
#ifdef HAVE_OPENCL
  integer :: ierr
#endif
  
  PUSH_SUB(X(accel_read_buffer_3))
  call profiling_in(prof_read, TOSTRING(X(CL_READ_BUFFER)))

  ! it does not make sense to read a buffer that the kernels cannot write
  ASSERT(this%flags /= ACCEL_MEM_READ_ONLY)

  fsize = size*R_SIZEOF
  offset_ = 0
  if(present(offset)) offset_ = offset*R_SIZEOF

  async_ = optional_default(async, .false.)

  if(fsize > 0) then
    
#ifdef HAVE_OPENCL
    call clEnqueueReadBuffer(accel%command_queue, this%mem, cl_bool(.true.), offset_, fsize, data(1, 1, 1), ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueReadBuffer")
#endif
#ifdef HAVE_CUDA
    call cuda_memcpy_dtoh(this%mem, data(1, 1, 1), fsize, offset_, async_)
#endif
    
    call profiling_count_transfers(size, data(1, 1, 1))
    if (.not. async_) call accel_finish()

  end if
  
  call profiling_out(prof_read)
  POP_SUB(X(accel_read_buffer_3))
  
end subroutine X(accel_read_buffer_3)

! ---------------------------------------------------------------------------

subroutine X(accel_set_kernel_arg_data)(kernel, narg, data)
  type(accel_kernel_t), intent(inout) :: kernel
  integer,              intent(in)    :: narg
  R_TYPE,               intent(in)    :: data
 
#ifdef HAVE_OPENCL 
  integer :: ierr
#endif

  ! no push_sub, called too frequently
#ifdef HAVE_CUDA
  call cuda_kernel_set_arg_value(kernel%arguments, data, narg, types_get_size(R_TYPE_VAL))
#endif
  
#ifdef HAVE_OPENCL
  call clSetKernelArg(kernel%kernel, narg, data, ierr)
  if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "set_kernel_arg_data")
#endif
  
end subroutine X(accel_set_kernel_arg_data)

! ---------------------------------------------------------------------------

subroutine X(accel_get_device_pointer_1)(host_pointer, device_pointer, dimensions)
  R_TYPE, pointer, dimension(:), intent(inout) :: host_pointer
  type(accel_mem_t),             intent(in)    :: device_pointer
  integer,         dimension(:), intent(in)    :: dimensions

  type(c_ptr) :: tmp_pointer

  PUSH_SUB(X(accel_get_device_pointer_1))

  ! move device pointer to fortran pointer for usage with CUDA-aware MPI
  call cuda_deref(device_pointer%mem, tmp_pointer)
  call c_f_pointer(tmp_pointer, host_pointer, dimensions)

  POP_SUB(X(accel_get_device_pointer_1))
end subroutine X(accel_get_device_pointer_1)

subroutine X(accel_get_device_pointer_2)(host_pointer, device_pointer, dimensions)
  R_TYPE, pointer, dimension(:,:), intent(inout) :: host_pointer
  type(accel_mem_t),               intent(in)    :: device_pointer
  integer,         dimension(:),   intent(in)    :: dimensions

  type(c_ptr) :: tmp_pointer

  PUSH_SUB(X(accel_get_device_pointer_2))

  ! move device pointer to fortran pointer for usage with CUDA-aware MPI
  call cuda_deref(device_pointer%mem, tmp_pointer)
  call c_f_pointer(tmp_pointer, host_pointer, dimensions)

  POP_SUB(X(accel_get_device_pointer_2))
end subroutine X(accel_get_device_pointer_2)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
