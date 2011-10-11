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

subroutine X(cube_function_alloc_FS)(cf)
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_alloc_FS))
  !Save memory if PFFT is used
  if (cf%fft_library /= PFFT_LIB) then
    ASSERT(.not.associated(cf%FS))
    ASSERT(associated(cf%fft))

    SAFE_ALLOCATE(cf%FS(1:cf%nx, 1:cf%n(2), 1:cf%n(3)))
  end if
  POP_SUB(X(cube_function_alloc_FS))
end subroutine X(cube_function_alloc_FS)


! ---------------------------------------------------------
subroutine X(cube_function_free_FS)(cf)
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_free_FS))
  !Save memory if PFFT is used
  if (cf%fft_library /= PFFT_LIB) then
    ASSERT(associated(cf%FS))
    SAFE_DEALLOCATE_P(cf%FS)
  end if

  POP_SUB(X(cube_function_free_FS))
end subroutine X(cube_function_free_FS)

! ---------------------------------------------------------
!> initializes the ffts. As the dimension of the fft may be adjusted, this
!! routine has to be called before allocating anything
subroutine X(cube_function_fft_init)(cf, sb)
  type(cube_function_t),     intent(inout) :: cf
  type(simul_box_t), intent(in)    :: sb

  PUSH_SUB(X(cube_function_fft_init))

  ASSERT(.not.associated(cf%fft))
  ASSERT(.not.associated(cf%X(RS)))
  ASSERT(.not.associated(cf%FS))
  SAFE_ALLOCATE(cf%fft)
  
  if (cf%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
    ASSERT(.not.associated(cf%pfft))
    
    SAFE_ALLOCATE(cf%pfft)
#ifdef R_TREAL
    call pfft_init(cf%n, sb%dim, fft_real, cf%pfft, optimize = .not.simul_box_is_periodic(sb))
    cf%nx = cf%n(1)/2 + 1
#else
    call pfft_init(cf%n, sb%dim, fft_complex, cf%pfft, .not.simul_box_is_periodic(sb))
    cf%nx = cf%n(1)
#endif 
#else 
    message(1) = "You have selected to use the PFFT library, "
    message(2) = "but it has not been linked."
    call messages_fatal(2)
#endif
    else
#ifdef R_TREAL
      call fft_init(cf%n, sb%dim, fft_real, cf%fft, optimize = .not.simul_box_is_periodic(sb))
      cf%nx = cf%n(1)/2 + 1
#else
      call fft_init(cf%n, sb%dim, fft_complex, cf%fft, optimize = .not.simul_box_is_periodic(sb))
      cf%nx = cf%n(1)
#endif  
    end if

  POP_SUB(X(cube_function_fft_init))
end subroutine X(cube_function_fft_init)

! ---------------------------------------------------------
!> The following routines convert the function between real space and Fourier space
!! Note that the dimensions of the function in FS are different depending on whether
!! f is real or complex, because the FFT representation is different (FFTW scheme).
subroutine X(cube_function_RS2FS)(cf)
  type(cube_function_t), intent(inout)  :: cf

  integer ::ii,jj,kk,index

  !Save memory if PFFT is used
  if (cf%fft_library /= PFFT_LIB) then
    ASSERT(associated(cf%X(RS)))
    if(.not.associated(cf%FS)) call X(cube_function_alloc_FS)(cf)
  end if
  
  if (cf%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
    call pfft_forward_3d(cf%pfft)
#endif
  else
    call X(fft_forward)(cf%fft, cf%X(RS), cf%FS)
  end if
   
end subroutine X(cube_function_RS2FS)


! ---------------------------------------------------------
subroutine X(cube_function_FS2RS)(cf)
  type(cube_function_t), intent(inout)  :: cf
  integer :: ii, jj, kk, index
  
  !Save memory if PFFT is used
  if (cf%fft_library /= PFFT_LIB) then
    ASSERT(associated(cf%FS))
    if(.not.associated(cf%X(RS))) call X(cube_function_alloc_RS)(cf)
  end if

  if (cf%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT 
    call pfft_backward_3d(cf%pfft)
#endif
  else
    call X(fft_backward)(cf%fft, cf%FS, cf%X(RS))
  end if
  
end subroutine X(cube_function_FS2RS)


! ---------------------------------------------------------
subroutine X(fourier_space_op_init)(this, cube, op)
  type(fourier_space_op_t), intent(out) :: this
  type(cube_function_t),            intent(in)  :: cube
  R_TYPE,                   intent(in)  :: op(:, :, :)

  integer :: ii, jj, kk

  nullify(this%dop)
  nullify(this%zop)
  SAFE_ALLOCATE(this%X(op)(1:cube%nx, 1:cube%n(2), 1:cube%n(3)))

  do kk = 1, cube%n(3)
    do jj = 1, cube%n(2)
      do ii = 1, cube%nx
        this%X(op)(ii, jj, kk) = op(ii, jj, kk)
      end do
    end do
  end do

end subroutine X(fourier_space_op_init)

! ---------------------------------------------------------
subroutine X(fourier_space_op_apply)(this, cube)
  type(fourier_space_op_t), intent(in)     :: this
  type(cube_function_t),    intent(inout)  :: cube
  
  integer :: ii, jj, kk, index

  type(profile_t), save :: prof_g, rs2fs_prof, fs2rs_prof, prof
  
  !Save memory if PFFT is used
  if (cube%fft_library /= PFFT_LIB) then
    call X(cube_function_alloc_FS)(cube)
  end if
  
  call profiling_in(prof, "OP_APPLY")
  call profiling_in(rs2fs_prof, "RS2FS")
  call X(cube_function_RS2FS)(cube)
  call profiling_out(rs2fs_prof)
  
  call profiling_in(prof_g,"G_APPLY")
  if (cube%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
    index = 1

    do kk =  cube%pfft%local_o_start(2),cube%pfft%local_o_start(2)+cube%pfft%local_no(2)-1
      do jj = cube%pfft%local_o_start(1), cube%pfft%local_o_start(1)+cube%pfft%local_no(1)-1
        do ii = cube%pfft%local_o_start(3), cube%pfft%local_o_start(3)+cube%pfft%local_no(3)-1 
          if (kk <= cube%nx)then
            cube%pfft%data_out(index)= cube%pfft%data_out(index)*this%X(op)(kk, jj, ii)
          else
            cube%pfft%data_out(index)= cube%pfft%data_out(index)*this%X(op)(cube%n(1)+2-kk, jj, ii)
          end if
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
          cube%FS(ii, jj, kk) = cube%FS(ii, jj, kk)*this%X(op)(ii, jj, kk)
        end do
      end do
    end do
  end if
  call profiling_out(prof_g)
  
  call profiling_in(fs2rs_prof, "FS2RS")
  call X(cube_function_FS2RS)(cube)
  call profiling_out(fs2rs_prof)
  call X(cube_function_free_FS)(cube)
  call profiling_out(prof)

end subroutine X(fourier_space_op_apply)

! ---------------------------------------------------------
!> The next two subroutines convert a function in Fourier space
!! between the normal mesh and the cube
!! Note that the function in the mesh should be defined
!! globally, not just in a partition (when running in
!! parallel in real-space domains).
! ---------------------------------------------------------
subroutine X(mesh_to_fourier) (mesh, mf, cf)
  type(mesh_t),  intent(in)    :: mesh
  CMPLX,         intent(in)    :: mf(:)   ! mf(mesh%np_global)
  type(cube_function_t), intent(inout) :: cf

  integer :: ip, ix, iy, iz

  ASSERT(associated(cf%FS))

  cf%FS = M_z0

  do ip = 1, mesh%np_global
    ix = pad_feq(mesh%idx%lxyz(ip, 1), cf%n(1), .false.)
    if(ix > cf%nx) cycle ! negative frequencies are redundant
    iy = pad_feq(mesh%idx%lxyz(ip, 2), cf%n(2), .false.)
    iz = pad_feq(mesh%idx%lxyz(ip, 3), cf%n(3), .false.)

    cf%FS(ix, iy, iz) = mf(ip)
  end do
end subroutine X(mesh_to_fourier)


! ---------------------------------------------------------
subroutine X(fourier_to_mesh) (mesh, cf, mf)
  type(mesh_t),  intent(in)  :: mesh
  type(cube_function_t), intent(in)  :: cf
  CMPLX,         intent(out) :: mf(:) ! mf(mesh%np_global)

  integer :: ip, ix, iy, iz

  ASSERT(associated(cf%FS))

  do ip = 1, mesh%np_global
    ix = pad_feq(mesh%idx%lxyz(ip, 1), cf%n(1), .false.)
    iy = pad_feq(mesh%idx%lxyz(ip, 2), cf%n(2), .false.)
    iz = pad_feq(mesh%idx%lxyz(ip, 3), cf%n(3), .false.)

#ifdef R_TREAL
    if(ix > cf%nx) then
      ix = pad_feq(-mesh%idx%lxyz(ip, 1), cf%n(1), .false.)
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
