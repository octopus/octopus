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

subroutine X(cube_function_alloc_FS)(cube, cf)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_alloc_FS))

  if (cube%fft_library /= FFTLIB_PFFT) then
    ASSERT(.not.associated(cf%FS))
    ASSERT(associated(cube%X(fftw)))

    SAFE_ALLOCATE(cf%FS(1:cube%fs_n(1), 1:cube%fs_n(2), 1:cube%fs_n(3)))
#ifdef HAVE_PFFT
  else
    ASSERT(.not.associated(cf%FS))
    cf%FS => cube%pfft%fs_data(1:cube%fs_n(3), 1:cube%fs_n(1), 1:cube%fs_n(2))
#endif
  end if

  POP_SUB(X(cube_function_alloc_FS))
end subroutine X(cube_function_alloc_FS)


! ---------------------------------------------------------
subroutine X(cube_function_free_FS)(cube, cf)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_free_FS))

  if (cube%fft_library /= FFTLIB_PFFT) then
    SAFE_DEALLOCATE_P(cf%FS)
  else
    nullify(cf%FS)
  end if

  POP_SUB(X(cube_function_free_FS))
end subroutine X(cube_function_free_FS)

! ---------------------------------------------------------
!> The following routines convert the function between real space and Fourier space
!! Note that the dimensions of the function in FS are different depending on whether
!! f is real or complex, because the FFT representation is different (FFTW scheme).
subroutine X(cube_function_RS2FS)(cube, cf)
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_RS2FS))

  ASSERT(cube%fft_library /= FFTLIB_NONE)

  if (cube%fft_library == FFTLIB_PFFT) then
#ifdef HAVE_PFFT
    ASSERT(associated(cf%X(RS)))
    ASSERT(associated(cf%FS))

    call pfft_forward_3d(cube%pfft)
#endif
  else
    ASSERT(associated(cf%X(RS)))
    ASSERT(associated(cf%FS))

    call X(fft_forward)(cube%X(fftw), cf%X(RS), cf%FS)
  end if

  POP_SUB(X(cube_function_RS2FS))
end subroutine X(cube_function_RS2FS)


! ---------------------------------------------------------
subroutine X(cube_function_FS2RS)(cube, cf)
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_FS2RS))

  ASSERT(cube%fft_library /= FFTLIB_NONE)

  if (cube%fft_library == FFTLIB_PFFT) then
#ifdef HAVE_PFFT
    ASSERT(associated(cf%X(RS)))
    ASSERT(associated(cf%FS))

    call pfft_backward_3d(cube%pfft)
#endif
  else
    ASSERT(associated(cf%X(RS)))
    ASSERT(associated(cf%FS))

    call X(fft_backward)(cube%X(fftw), cf%FS, cf%X(RS))
  end if

  POP_SUB(X(cube_function_FS2RS))

end subroutine X(cube_function_FS2RS)


! ---------------------------------------------------------
subroutine X(fourier_space_op_init)(this, cube, op)
  type(fourier_space_op_t), intent(out) :: this
  type(cube_t),             intent(in)  :: cube
  R_TYPE,                   intent(in)  :: op(:, :, :)

  integer :: ii, jj, kk

  PUSH_SUB(X(fourier_space_op_init))

  ASSERT(cube%fft_library /= FFTLIB_NONE)

  nullify(this%dop)
  nullify(this%zop)
  SAFE_ALLOCATE(this%X(op)(1:cube%fs_n(1), 1:cube%fs_n(2), 1:cube%fs_n(3)))
  forall (kk = 1:cube%fs_n(3), jj = 1:cube%fs_n(2), ii = 1:cube%fs_n(1)) 
    this%X(op)(ii, jj, kk) = op(ii, jj, kk)
  end forall

  POP_SUB(X(fourier_space_op_init))
end subroutine X(fourier_space_op_init)

! ---------------------------------------------------------
subroutine X(fourier_space_op_apply)(this, cube, cf)
  type(fourier_space_op_t), intent(in)     :: this
  type(cube_t),             intent(inout)  :: cube
  type(cube_function_t),    intent(inout)  :: cf
  
  integer :: ii, jj, kk, index

  type(profile_t), save :: prof_g, rs2fs_prof, fs2rs_prof, prof

  PUSH_SUB(X(fourier_space_op_apply))

  ASSERT(cube%fft_library /= FFTLIB_NONE)

  call X(cube_function_alloc_FS)(cube, cf)

  call profiling_in(prof, "OP_APPLY")
  call profiling_in(rs2fs_prof, "RS2FS")
  call X(cube_function_RS2FS)(cube, cf)
  call profiling_out(rs2fs_prof)
  
  call profiling_in(prof_g,"G_APPLY")
  if (cube%fft_library == FFTLIB_PFFT) then
    !Note that the function in fourier space returned by PFFT is transposed
    do kk = 1, cube%fs_n(3)
      do jj = 1, cube%fs_n(2)
        do ii = 1, cube%fs_n(1)
          cf%FS(kk, ii, jj) = cf%FS(kk, ii, jj)*this%X(op) (ii, jj, kk)
        end do
      end do
    end do
  else
    do kk = 1, cube%n(3)
      do jj = 1, cube%n(2)
        do ii = 1, cube%nx
          cf%FS(ii, jj, kk) = cf%FS(ii, jj, kk)*this%X(op)(ii, jj, kk)
        end do
      end do
    end do
  end if
  call profiling_out(prof_g)
  
  call profiling_in(fs2rs_prof, "FS2RS")
  call X(cube_function_FS2RS)(cube, cf)
  call profiling_out(fs2rs_prof)
  call X(cube_function_free_FS)(cube, cf)
  call profiling_out(prof)

  POP_SUB(X(fourier_space_op_apply))
end subroutine X(fourier_space_op_apply)

! ---------------------------------------------------------
!> The next two subroutines convert a function in Fourier space
!! between the normal mesh and the cube
!! Note that the function in the mesh should be defined
!! globally, not just in a partition (when running in
!! parallel in real-space domains).
! ---------------------------------------------------------
subroutine X(mesh_to_fourier) (mesh, mf, cube, cf)
  type(mesh_t),  intent(in)    :: mesh
  CMPLX,         intent(in)    :: mf(:)   ! mf(mesh%np_global)
  type(cube_t),  intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  integer :: ip, ix, iy, iz

  PUSH_SUB(X(mesh_to_fourier))

  ASSERT(associated(cf%FS))

  cf%FS = M_z0

  do ip = 1, mesh%np_global
    ix = pad_feq(mesh%idx%lxyz(ip, 1), cube%n(1), .false.)
    if(ix > cube%nx) cycle ! negative frequencies are redundant
    iy = pad_feq(mesh%idx%lxyz(ip, 2), cube%n(2), .false.)
    iz = pad_feq(mesh%idx%lxyz(ip, 3), cube%n(3), .false.)

    cf%FS(ix, iy, iz) = mf(ip)
  end do

  POP_SUB(X(mesh_to_fourier))
end subroutine X(mesh_to_fourier)


! ---------------------------------------------------------
subroutine X(fourier_to_mesh) (cube, cf, mesh, mf)
  type(cube_t),  intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  type(mesh_t),  intent(in)  :: mesh
  CMPLX,         intent(out) :: mf(:) ! mf(mesh%np_global)

  integer :: ip, ix, iy, iz

  PUSH_SUB(X(fourier_to_mesh))

  ASSERT(associated(cf%FS))

  do ip = 1, mesh%np_global
    ix = pad_feq(mesh%idx%lxyz(ip, 1), cube%n(1), .false.)
    iy = pad_feq(mesh%idx%lxyz(ip, 2), cube%n(2), .false.)
    iz = pad_feq(mesh%idx%lxyz(ip, 3), cube%n(3), .false.)

#ifdef R_TREAL
    if(ix > cube%nx) then
      ix = pad_feq(-mesh%idx%lxyz(ip, 1), cube%n(1), .false.)
      mf(ip) = conjg(cf%FS(ix, iy, iz))
    else
      mf(ip) = cf%FS(ix, iy, iz)
    end if
#else
    mf(ip) = cf%FS(ix, iy, iz)
#endif
  end do

  POP_SUB(X(fourier_to_mesh))
end subroutine X(fourier_to_mesh)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
