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

  !Save memory if PFFT is used
  if (cube%fft_library /= PFFT_LIB) then
    ASSERT(.not.associated(cf%FS))
    ASSERT(associated(cube%X(fftw)))

    SAFE_ALLOCATE(cf%FS(1:cube%nx, 1:cube%n(2), 1:cube%n(3)))
  end if

  POP_SUB(X(cube_function_alloc_FS))
end subroutine X(cube_function_alloc_FS)


! ---------------------------------------------------------
subroutine X(cube_function_free_FS)(cf)
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_free_FS))

  !Save memory if PFFT is used
  if (associated(cf%FS)) then
    SAFE_DEALLOCATE_P(cf%FS)
  end if

  POP_SUB(X(cube_function_free_FS))
end subroutine X(cube_function_free_FS)

! ---------------------------------------------------------
!> The following routines convert the function between real space and Fourier space
!! Note that the dimensions of the function in FS are different depending on whether
!! f is real or complex, because the FFT representation is different (FFTW scheme).
subroutine X(cube_function_RS2FS)(cube, cf)
  type(cube_t),          intent(inout)  :: cube
  type(cube_function_t), intent(inout)  :: cf

  !Save memory if PFFT is used
  if (cube%fft_library /= PFFT_LIB) then
    ASSERT(associated(cf%X(RS)))
    if(.not.associated(cf%FS)) call X(cube_function_alloc_FS)(cube, cf)
  end if
  
  if (cube%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
    call pfft_forward_3d(cube%pfft)
#endif
  else
    call X(fft_forward)(cube%X(fftw), cf%X(RS), cf%FS)
  end if
   
end subroutine X(cube_function_RS2FS)


! ---------------------------------------------------------
subroutine X(cube_function_FS2RS)(cube, cf)
  type(cube_t),          intent(inout)  :: cube
  type(cube_function_t), intent(inout)  :: cf

  !Save memory if PFFT is used
  if (cube%fft_library /= PFFT_LIB) then
    ASSERT(associated(cf%FS))
    if(.not.associated(cf%X(RS))) call X(cube_function_alloc_RS)(cube, cf)
  end if

  if (cube%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT 
    call pfft_backward_3d(cube%pfft)
#endif
  else
    call X(fft_backward)(cube%X(fftw), cf%FS, cf%X(RS))
  end if

end subroutine X(cube_function_FS2RS)


! ---------------------------------------------------------
subroutine X(fourier_space_op_init)(this, cube, op)
  type(fourier_space_op_t), intent(out) :: this
  type(cube_t),             intent(in)  :: cube
  R_TYPE,                   intent(in)  :: op(:, :, :)

  integer :: ii, jj, kk
  integer :: start(3),last(3)

  nullify(this%dop)
  nullify(this%zop)
  if (cube%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
    start = cube%pfft%local_o_start
    last = cube%pfft%local_o_start + cube%pfft%local_no - 1
    SAFE_ALLOCATE(this%X(op)(start(2):last(2),start(1):last(1),start(3):last(3)))
     do kk =  cube%pfft%local_o_start(2),cube%pfft%local_o_start(2)+cube%pfft%local_no(2)-1
      do jj = cube%pfft%local_o_start(1), cube%pfft%local_o_start(1)+cube%pfft%local_no(1)-1
        do ii = cube%pfft%local_o_start(3), cube%pfft%local_o_start(3)+cube%pfft%local_no(3)-1
          this%X(op)(kk, jj, ii) = op(kk-start(2)+1, jj-start(1)+1, ii-start(3)+1)
        end do
      end do
    end do
#endif
  else
    SAFE_ALLOCATE(this%X(op)(1:cube%nx, 1:cube%n(2), 1:cube%n(3)))
    do kk = 1, cube%n(3)
      do jj = 1, cube%n(2)
        do ii = 1, cube%nx
          this%X(op)(ii, jj, kk) = op(ii, jj, kk)
        end do
      end do
    end do
  end if
end subroutine X(fourier_space_op_init)

! ---------------------------------------------------------
subroutine X(fourier_space_op_apply)(this, cube, cf)
  type(fourier_space_op_t), intent(in)     :: this
  type(cube_t),             intent(inout)  :: cube
  type(cube_function_t),    intent(inout)  :: cf
  
  integer :: ii, jj, kk, index
  integer :: start(3), last(3)

  type(profile_t), save :: prof_g, rs2fs_prof, fs2rs_prof, prof
  
  !Save memory if PFFT is used
  if (cube%fft_library /= PFFT_LIB) then
    call X(cube_function_alloc_FS)(cube, cf)
  end if

  call profiling_in(prof, "OP_APPLY")
  call profiling_in(rs2fs_prof, "RS2FS")
  call X(cube_function_RS2FS)(cube, cf)
  call profiling_out(rs2fs_prof)
  
  call profiling_in(prof_g,"G_APPLY")
  if (cube%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
    index = 1
    start = cube%pfft%local_o_start
    last = cube%pfft%local_o_start + cube%pfft%local_no - 1
    do kk =  cube%pfft%local_o_start(2),cube%pfft%local_o_start(2)+cube%pfft%local_no(2)-1
      do jj = cube%pfft%local_o_start(1), cube%pfft%local_o_start(1)+cube%pfft%local_no(1)-1
        do ii = cube%pfft%local_o_start(3), cube%pfft%local_o_start(3)+cube%pfft%local_no(3)-1 
          cube%pfft%data_out(index)= cube%pfft%data_out(index)*this%X(op)(kk, jj,ii)
          index=index+1
        end do
      end do
    end do
#else
    write(message(1),'(a)')'You have selected the PFFT for FFT, but it is not compiled.'
    call messages_fatal(1)
#endif
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
  call X(cube_function_free_FS)(cf)
  call profiling_out(prof)

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

  ASSERT(associated(cf%FS))

  cf%FS = M_z0

  do ip = 1, mesh%np_global
    ix = pad_feq(mesh%idx%lxyz(ip, 1), cube%n(1), .false.)
    if(ix > cube%nx) cycle ! negative frequencies are redundant
    iy = pad_feq(mesh%idx%lxyz(ip, 2), cube%n(2), .false.)
    iz = pad_feq(mesh%idx%lxyz(ip, 3), cube%n(3), .false.)

    cf%FS(ix, iy, iz) = mf(ip)
  end do
end subroutine X(mesh_to_fourier)


! ---------------------------------------------------------
subroutine X(fourier_to_mesh) (cube, cf, mesh, mf)
  type(cube_t),  intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  type(mesh_t),  intent(in)  :: mesh
  CMPLX,         intent(out) :: mf(:) ! mf(mesh%np_global)

  integer :: ip, ix, iy, iz

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

end subroutine X(fourier_to_mesh)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
