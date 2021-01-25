!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M.Oliveira
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
!> The following routines convert the function between real space and Fourier space
!! Note that the dimensions of the function in fs are different depending on whether
!! f is real or complex, because the FFT representation is different (FFTW scheme).
subroutine X(cube_function_rs2fs)(cube, cf)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_rs2fs))

  ASSERT(associated(cube%fft))
  ASSERT(cube%fft%library /= FFTLIB_NONE)

  if(cf%in_device_memory) then
    call X(fft_forward)(cube%fft, cf%real_space_buffer, cf%fourier_space_buffer)
  else
    ASSERT(associated(cf%X(rs)))
    ASSERT(associated(cf%fs))
    
    call X(fft_forward)(cube%fft, cf%X(rs), cf%fs)
  end if

  POP_SUB(X(cube_function_rs2fs))
end subroutine X(cube_function_rs2fs)

! ---------------------------------------------------------
subroutine X(cube_function_fs2rs)(cube, cf)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_fs2rs))

  ASSERT(associated(cube%fft))
  ASSERT(cube%fft%library /= FFTLIB_NONE)

  if(cf%in_device_memory) then
    call X(fft_backward)(cube%fft, cf%fourier_space_buffer, cf%real_space_buffer)
  else
    ASSERT(associated(cf%X(rs)))
    ASSERT(associated(cf%fs))

    call X(fft_backward)(cube%fft, cf%fs, cf%X(rs))
  end if

  POP_SUB(X(cube_function_fs2rs))
end subroutine X(cube_function_fs2rs)

! ---------------------------------------------------------
subroutine X(fourier_space_op_init)(this, cube, op, in_device)
  type(fourier_space_op_t), intent(inout) :: this
  type(cube_t),             intent(in)    :: cube
  R_TYPE,                   intent(in)    :: op(:, :, :)
  logical, optional,        intent(in)    :: in_device

  integer :: ii, jj, kk, ii_linear, size
  R_TYPE, allocatable :: op_linear(:)

  PUSH_SUB(X(fourier_space_op_init))

  ASSERT(associated(cube%fft))
  ASSERT(cube%fft%library /= FFTLIB_NONE)

#ifdef R_TREAL
  this%real_op = .true.
#else
  this%real_op = .false.
#endif

  if(cube%fft%library /= FFTLIB_ACCEL .or. .not. optional_default(in_device, .true.)) then
    this%in_device_memory = .false.
    SAFE_ALLOCATE(this%X(op)(1:cube%fs_n(1), 1:cube%fs_n(2), 1:cube%fs_n(3)))
    do kk = 1,cube%fs_n(3)
      do jj = 1,cube%fs_n(2)
        do ii = 1,cube%fs_n(1)
          this%X(op)(ii, jj, kk) = op(ii, jj, kk)
        end do
      end do
    end do
  else
    this%in_device_memory = .true.

    ASSERT(all(cube%fs_n(1:3) == ubound(op)))
    
    size = product(cube%fs_n(1:3))

    SAFE_ALLOCATE(op_linear(1:size))

    do kk = 1, cube%fs_n(3)
      do jj = 1, cube%fs_n(2)
        do ii = 1, cube%fs_n(1)
          ii_linear = 1 + (ii - 1)*cube%fft%stride_fs(1) + (jj - 1)*cube%fft%stride_fs(2) + (kk - 1)*cube%fft%stride_fs(3)
          op_linear(ii_linear) = fft_scaling_factor(cube%fft)*op(ii, jj, kk)
        end do
      end do
    end do

    call accel_create_buffer(this%op_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, size)
    call accel_write_buffer(this%op_buffer, size, op_linear)
    
    SAFE_DEALLOCATE_A(op_linear)
  end if

  POP_SUB(X(fourier_space_op_init))
end subroutine X(fourier_space_op_init)

! ---------------------------------------------------------
!> Applies a multiplication factor to the Fourier space grid.
!! This is a local function.
subroutine X(fourier_space_op_apply)(this, cube, cf)
  type(fourier_space_op_t), intent(in)     :: this
  type(cube_t),             intent(in)     :: cube
  type(cube_function_t),    intent(inout)  :: cf
  
  integer :: ii, jj, kk
  integer :: bsize

  type(profile_t), save :: prof_g, prof

  PUSH_SUB(X(fourier_space_op_apply))

  ASSERT(associated(cube%fft))
  ASSERT(cube%fft%library /= FFTLIB_NONE)
  ASSERT(cf%in_device_memory .eqv. this%in_device_memory)

  call cube_function_alloc_fs(cube, cf)

  call profiling_in(prof, TOSTRING(X(OP_APPLY)))

  call X(cube_function_rs2fs)(cube, cf)
   
  call profiling_in(prof_g,TOSTRING(X(G_APPLY)))

  if(cube%fft%library == FFTLIB_PFFT) then
    !Note that the function in fourier space returned by PFFT is transposed
    ! Probably in this case this%X(op) should be also transposed
    if(this%real_op) then
      !$omp parallel do private(ii, jj, kk)
      do ii = 1, cube%fs_n(1)
        do jj = 1, cube%fs_n(2)
          do kk = 1, cube%fs_n(3)
            cf%fs(kk, ii, jj) = cf%fs(kk, ii, jj)*this%dop(ii, jj, kk)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(ii, jj, kk)
      do ii = 1, cube%fs_n(1)
        do jj = 1, cube%fs_n(2)
          do kk = 1, cube%fs_n(3)
            cf%fs(kk, ii, jj) = cf%fs(kk, ii, jj)*this%X(op) (ii, jj, kk)
          end do
        end do
      end do
      !$omp end parallel do
    end if
  else if(cube%fft%library == FFTLIB_FFTW .or. .not. this%in_device_memory) then
    if(this%real_op) then
      !$omp parallel do private(ii, jj, kk)
      do kk = 1, cube%fs_n(3)
        do jj = 1, cube%fs_n(2)
          do ii = 1, cube%fs_n(1)
            cf%fs(ii, jj, kk) = cf%fs(ii, jj, kk)*this%dop(ii, jj, kk)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(ii, jj, kk)
      do kk = 1, cube%fs_n(3)
        do jj = 1, cube%fs_n(2)
          do ii = 1, cube%fs_n(1)
            cf%fs(ii, jj, kk) = cf%fs(ii, jj, kk)*this%X(op)(ii, jj, kk)
          end do
        end do
      end do
      !$omp end parallel do
    end if
  else if(cube%fft%library == FFTLIB_ACCEL) then
    call accel_set_kernel_arg(X(zmul), 0, product(cube%fs_n(1:3)))
    call accel_set_kernel_arg(X(zmul), 1, this%op_buffer)
    call accel_set_kernel_arg(X(zmul), 2, cf%fourier_space_buffer)
    bsize = accel_kernel_workgroup_size(X(zmul))
    call accel_kernel_run(X(zmul), (/pad(product(cube%fs_n(1:3)), bsize)/), (/bsize/))
    call accel_finish()
  end if

#ifdef R_TREAL
  call profiling_count_operations(2*R_MUL*product(cube%fs_n(1:3)))
#else
  call profiling_count_operations(R_MUL*product(cube%fs_n(1:3)))
#endif

  call profiling_out(prof_g)
  
  call X(cube_function_fs2rs)(cube, cf)
  call cube_function_free_fs(cube, cf)
  call profiling_out(prof)

  POP_SUB(X(fourier_space_op_apply))
end subroutine X(fourier_space_op_apply)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
