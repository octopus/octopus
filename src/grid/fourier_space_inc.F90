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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

! ---------------------------------------------------------
!> The following routines convert the function between real space and Fourier space
!! Note that the dimensions of the function in fs are different depending on whether
!! f is real or complex, because the FFT representation is different (FFTW scheme).
subroutine X(cube_function_rs2fs)(cube, cf)
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_rs2fs))

  ASSERT(associated(cube%fft))
  ASSERT(cube%fft%library /= FFTLIB_NONE)

  if(cf%in_device_memory) then
    call X(fft_forward_cl)(cube%fft, cf%real_space_buffer, cf%fourier_space_buffer)
  else
    ASSERT(associated(cf%X(rs)))
    ASSERT(associated(cf%fs))
    
    call X(fft_forward)(cube%fft, cf%X(rs), cf%fs)
  end if

  POP_SUB(X(cube_function_rs2fs))
end subroutine X(cube_function_rs2fs)

! ---------------------------------------------------------
subroutine X(cube_function_fs2rs)(cube, cf)
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_fs2rs))

  ASSERT(associated(cube%fft))
  ASSERT(cube%fft%library /= FFTLIB_NONE)

  if(cf%in_device_memory) then
    call X(fft_backward_cl)(cube%fft, cf%fourier_space_buffer, cf%real_space_buffer)
  else
    ASSERT(associated(cf%X(rs)))
    ASSERT(associated(cf%fs))

    call X(fft_backward)(cube%fft, cf%fs, cf%X(rs))
  end if

  POP_SUB(X(cube_function_fs2rs))
end subroutine X(cube_function_fs2rs)

! ---------------------------------------------------------
subroutine X(fourier_space_op_init)(this, cube, op, in_device)
  type(fourier_space_op_t), intent(out) :: this
  type(cube_t),             intent(in)  :: cube
  R_TYPE,                   intent(in)  :: op(:, :, :)
  logical, optional,        intent(in)  :: in_device

  integer :: ii, jj, kk, ii_linear, size
  R_TYPE, allocatable :: op_linear(:)

  PUSH_SUB(X(fourier_space_op_init))

  ASSERT(associated(cube%fft))
  ASSERT(cube%fft%library /= FFTLIB_NONE)

  nullify(this%dop)
  nullify(this%zop)

  if(cube%fft%library /= FFTLIB_CLAMD .or. .not. optional_default(in_device, .true.)) then
    this%in_device_memory = .false.
    SAFE_ALLOCATE(this%X(op)(1:cube%fs_n(1), 1:cube%fs_n(2), 1:cube%fs_n(3)))
    forall (kk = 1:cube%fs_n(3), jj = 1:cube%fs_n(2), ii = 1:cube%fs_n(1)) 
      this%X(op)(ii, jj, kk) = op(ii, jj, kk)
    end forall
  else
    this%in_device_memory = .true.

    ASSERT(all(cube%fs_n(1:3) == ubound(op)))
    
    size = product(cube%fs_n(1:3))

    SAFE_ALLOCATE(op_linear(1:size))

    do kk = 1, cube%fs_n(3)
      do jj = 1, cube%fs_n(2)
        do ii = 1, cube%fs_n(1)
          ii_linear = 1 + (ii - 1)*cube%fft%stride_fs(1) + (jj - 1)*cube%fft%stride_fs(2) + (kk - 1)*cube%fft%stride_fs(3)
          op_linear(ii_linear) = op(ii, jj, kk)
        end do
      end do
    end do

#ifdef HAVE_OPENCL
    call opencl_create_buffer(this%op_buffer, CL_MEM_READ_ONLY, R_TYPE_VAL, size)
    call opencl_write_buffer(this%op_buffer, size, op_linear)
#endif
    
    SAFE_DEALLOCATE_A(op_linear)
  end if

  POP_SUB(X(fourier_space_op_init))
end subroutine X(fourier_space_op_init)

! ---------------------------------------------------------
!> Applies a multiplication factor to the Fourier space grid.
!! This is a local function.
subroutine X(fourier_space_op_apply)(this, cube, cf)
  type(fourier_space_op_t), intent(in)     :: this
  type(cube_t),             intent(inout)  :: cube
  type(cube_function_t),    intent(inout)  :: cf
  
  integer :: ii, jj, kk
#ifdef HAVE_OPENCL
  integer :: bsize
#endif

  type(profile_t), save :: prof_g, prof

  PUSH_SUB(X(fourier_space_op_apply))

  ASSERT(associated(cube%fft))
  ASSERT(cube%fft%library /= FFTLIB_NONE)
  ASSERT(cf%in_device_memory .eqv. this%in_device_memory)

  call cube_function_alloc_fs(cube, cf)

  call profiling_in(prof, "OP_APPLY")

  call X(cube_function_rs2fs)(cube, cf)
   
  call profiling_in(prof_g,"G_APPLY")

  if(cube%fft%library == FFTLIB_PFFT) then
    !Note that the function in fourier space returned by PFFT is transposed
    ! Probably in this case this%X(op) should be also transposed
    !$omp parallel do private(ii, jj, kk)
    do ii = 1, cube%fs_n(1)
      do jj = 1, cube%fs_n(2)
        do kk = 1, cube%fs_n(3)
          cf%fs(kk, ii, jj) = cf%fs(kk, ii, jj)*this%X(op) (ii, jj, kk)
        end do
      end do
    end do
    !$omp end parallel do
  else if(cube%fft%library == FFTLIB_FFTW .or. .not. this%in_device_memory) then
    !$omp parallel do private(ii, jj, kk)
    do kk = 1, cube%fs_n(3)
      do jj = 1, cube%fs_n(2)
        do ii = 1, cube%fs_n(1)
          cf%fs(ii, jj, kk) = cf%fs(ii, jj, kk)*this%X(op)(ii, jj, kk)
        end do
      end do
    end do
    !$omp end parallel do
  else if(cube%fft%library == FFTLIB_CLAMD) then
#ifdef HAVE_OPENCL
    call opencl_set_kernel_arg(X(zmul), 0, product(cube%fs_n(1:3)))
    call opencl_set_kernel_arg(X(zmul), 1, this%op_buffer)
    call opencl_set_kernel_arg(X(zmul), 2, cf%fourier_space_buffer)
    bsize = opencl_kernel_workgroup_size(X(zmul))
    call opencl_kernel_run(X(zmul), (/pad(product(cube%fs_n(1:3)), bsize)/), (/bsize/))
    call opencl_finish()
#endif
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

! ---------------------------------------------------------
!> The next two subroutines convert a function in Fourier space
!! between the normal mesh and the cube.
!!
!! Note that the function in the mesh should be defined
!! globally, not just in a partition (when running in
!! parallel in real-space domains).
! ---------------------------------------------------------
subroutine X(mesh_to_fourier) (mesh, mf, cube, cf)
  type(mesh_t),  intent(in)    :: mesh
  CMPLX,         intent(in)    :: mf(:)   !< mf(mesh%np_global)
  type(cube_t),  intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  integer :: ip, ix, iy, iz

  PUSH_SUB(X(mesh_to_fourier))

  ASSERT(associated(cube%fft))
  ASSERT(associated(cf%fs))

  cf%fs = M_z0

  do ip = 1, mesh%np_global
    ix = pad_feq(mesh%idx%lxyz(ip, 1), cube%rs_n_global(1), .false.)
    if(ix > cube%fs_n_global(1)) cycle ! negative frequencies are redundant
    iy = pad_feq(mesh%idx%lxyz(ip, 2), cube%rs_n_global(2), .false.)
    iz = pad_feq(mesh%idx%lxyz(ip, 3), cube%rs_n_global(3), .false.)

    cf%fs(ix, iy, iz) = mf(ip)
  end do

  POP_SUB(X(mesh_to_fourier))
end subroutine X(mesh_to_fourier)


! ---------------------------------------------------------
subroutine X(fourier_to_mesh) (cube, cf, mesh, mf)
  type(cube_t),  intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  type(mesh_t),  intent(in)  :: mesh
  CMPLX,         intent(out) :: mf(:) !< mf(mesh%np_global)

  integer :: ip, ix, iy, iz

  PUSH_SUB(X(fourier_to_mesh))

  ASSERT(associated(cube%fft))
  ASSERT(associated(cf%fs))

  do ip = 1, mesh%np_global
    ix = pad_feq(mesh%idx%lxyz(ip, 1), cube%rs_n_global(1), .false.)
    iy = pad_feq(mesh%idx%lxyz(ip, 2), cube%rs_n_global(2), .false.)
    iz = pad_feq(mesh%idx%lxyz(ip, 3), cube%rs_n_global(3), .false.)

#ifdef R_TREAL
    if(ix > cube%fs_n_global(1)) then
      ix = pad_feq(-mesh%idx%lxyz(ip, 1), cube%rs_n_global(1), .false.)
      mf(ip) = conjg(cf%fs(ix, iy, iz))
    else
      mf(ip) = cf%fs(ix, iy, iz)
    end if
#else
    mf(ip) = cf%fs(ix, iy, iz)
#endif
  end do

  POP_SUB(X(fourier_to_mesh))
end subroutine X(fourier_to_mesh)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
