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
subroutine X(cube_function_alloc_RS)(cf)
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_alloc_RS))
  !Safe memory if PFFT is used  
  if (cf%fft_library /= PFFT_LIB) then
    ASSERT(.not.associated(cf%X(RS)))
    SAFE_ALLOCATE(cf%X(RS)(1:cf%n(1), 1:cf%n(2), 1:cf%n(3)))
  end if
  	
  POP_SUB(X(cube_function_alloc_RS))
end subroutine X(cube_function_alloc_RS)


! ---------------------------------------------------------
subroutine X(cube_function_free_RS)(cf)
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_free_RS))
  !Safe memory if PFFT is used  
  if (cf%fft_library /= PFFT_LIB) then
    ASSERT(associated(cf%X(RS)))
    SAFE_DEALLOCATE_P(cf%X(RS))
  end if

  POP_SUB(X(cube_function_free_RS))
end subroutine X(cube_function_free_RS)

! ---------------------------------------------------------
!> The next two subroutines convert a function between the normal
!! mesh and the cube.
!! Note that the function in the mesh should be defined
!! globally, not just in a partition (when running in
!! parallel in real-space domains).
! ---------------------------------------------------------

subroutine X(mesh_to_cube)(mesh, mf, cf, local)
  type(mesh_t),          intent(in)    :: mesh
  R_TYPE,  target,       intent(in)    :: mf(:)  !< mf(mesh%np_global). It is the global rho
  type(cube_function_t), intent(inout) :: cf
  logical, optional,     intent(in)    :: local  !< If .true. the mf array is a local array. Considered .false. if not present.

  integer :: ip, ix, iy, iz, center(3), ixyz(1:3)
  integer :: im, ii, nn
  logical :: local_
  R_TYPE, pointer :: gmf(:)

  PUSH_SUB(X(mesh_to_cube))
  call profiling_in(prof_m2c, "MESH_TO_CUBE")

  local_ = optional_default(local, .false.) .and. mesh%parallel_in_domains
  
  if(local_) then
    ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)
    SAFE_ALLOCATE(gmf(1:mesh%np_global))
#ifdef HAVE_MPI
    call X(vec_allgather)(mesh%vp, gmf, mf)
#endif
  else
    ASSERT(ubound(mf, dim = 1) == mesh%np_global .or. ubound(mf, dim = 1) == mesh%np_part_global)
    gmf => mf
  end if

  ASSERT(associated(cf%X(RS)))

  center(1:3) = cf%n(1:3)/2 + 1

  cf%X(RS) = M_ZERO

  ASSERT(associated(mesh%cube_map%map))
  ASSERT(mesh%sb%dim <= 3)

  if(associated(mesh%idx%lxyz)) then

    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)
      
      ix = mesh%idx%lxyz(ip, 1) + center(1)
      iy = mesh%idx%lxyz(ip, 2) + center(2)
      iz = mesh%idx%lxyz(ip, 3) + center(3)
      forall(ii = 0:nn - 1) cf%X(RS)(ix, iy, iz + ii) = gmf(ip + ii)
    end do

  else

    ixyz = 0
    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      call index_to_coords(mesh%idx, mesh%sb%dim, ip, ixyz)
      ixyz = ixyz + center
      forall(ii = 0:nn - 1) cf%X(RS)(ixyz(1), ixyz(2), ixyz(3) + ii) = gmf(ip + ii)
    end do

  end if

  if(local_) then
    SAFE_DEALLOCATE_P(gmf)
  end if

  call profiling_count_transfers(mesh%np_global, mf(1))

  call profiling_out(prof_m2c)
  POP_SUB(X(mesh_to_cube))
end subroutine X(mesh_to_cube)

! ---------------------------------------------------------
subroutine X(cube_to_mesh) (mesh, cf, mf)
  type(mesh_t),  intent(in)  :: mesh
  type(cube_function_t), intent(in)  :: cf
  R_TYPE,        intent(out) :: mf(:)  ! mf(mesh%np_global)

  integer :: ip, ix, iy, iz, center(3), ixyz(1:3)
  integer :: im, ii, nn

  PUSH_SUB(X(cube_to_mesh))

  call profiling_in(prof_c2m, "CUBE_TO_MESH")

  ASSERT(associated(cf%X(RS)))

  center(1:3) = cf%n(1:3)/2 + 1

  ASSERT(associated(mesh%cube_map%map))

  if(associated(mesh%idx%lxyz)) then

    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)
      ix = mesh%idx%lxyz(ip, 1) + center(1)
      iy = mesh%idx%lxyz(ip, 2) + center(2)
      iz = mesh%idx%lxyz(ip, 3) + center(3)
      forall(ii = 0:nn - 1) mf(ip + ii) = cf%X(RS)(ix, iy, iz + ii)
    end do

  else

    ixyz = 0

    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      call index_to_coords(mesh%idx, mesh%sb%dim, ip, ixyz)
      ixyz = ixyz + center

      forall(ii = 0:nn - 1) mf(ip + ii) = cf%X(RS)(ixyz(1), ixyz(2), ixyz(3) + ii)
    end do

  end if

  call profiling_count_transfers(mesh%np_global, mf(1))

  call profiling_out(prof_c2m)

  POP_SUB(X(cube_to_mesh))

end subroutine X(cube_to_mesh)

#ifdef HAVE_PFFT
! ---------------------------------------------------------
!> The next two subroutines convert a function between the normal
!! mesh and the cube in parallel.
!! "Note that the function in the mesh should be defined
!! globally, not just in a partition (when running in
!! parallel in real-space domains)."
! ---------------------------------------------------------
subroutine X(mesh_to_cube_parallel)(mesh, mf, cf, local)
  type(mesh_t),          intent(in)    :: mesh
  R_TYPE,  target,       intent(in)    :: mf(:)  !< mf(mesh%np_global)
  type(cube_function_t), intent(inout) :: cf
  logical, optional,     intent(in)    :: local  !< If .true. the mf array is a local array. Considered .false. if not present.

  integer :: ip, ix, iy, iz, center(3)
  integer :: im, ii, jj, kk, nn
  integer :: min_x, min_y, min_z, max_x,max_y,max_z
  logical :: local_
  R_TYPE, pointer :: gmf(:)

  PUSH_SUB(X(mesh_to_cube_parallel))
  call profiling_in(prof_m2c, "MESH_TO_CUBE_PARALLEL")

  local_ = optional_default(local, .false.) .and. mesh%parallel_in_domains
  
  if(local_) then
    ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)
    SAFE_ALLOCATE(gmf(1:mesh%np_global))
#ifdef HAVE_MPI
    call X(vec_allgather)(mesh%vp, gmf, mf)
#endif
  else
    ASSERT(ubound(mf, dim = 1) == mesh%np_global .or. ubound(mf, dim = 1) == mesh%np_part_global)
    gmf => mf
  end if

  center(1:3) = cf%n(1:3)/2 + 1
  
  !The author of PFFT library (M. Pippig) says that  is not needed to initialize to zero the input matrix. 
  !Pippig, M. :An Efficient and Flexible Parallel FFT Implementation Based on FFTW (page 7)

  ! Save the limit values
  min_x = cf%pfft%local_i_start(1)
  min_y = cf%pfft%local_i_start(2)
  min_z = cf%pfft%local_i_start(3)
  max_x = cf%pfft%local_i_start(1)+cf%pfft%local_ni(1)
  max_y = cf%pfft%local_i_start(2)+cf%pfft%local_ni(2)
  max_z = cf%pfft%local_i_start(3)+cf%pfft%local_ni(3)
    
  ! Do the actual transform, only for the output values
  do im = 1, mesh%cube_map%nmap
    ip = mesh%cube_map%map(MCM_POINT, im)
    nn = mesh%cube_map%map(MCM_COUNT, im)
    
    ix = mesh%idx%lxyz(ip, 1) + center(1)
    if (ix >= min_x .and. ix < max_x) then
      iy = mesh%idx%lxyz(ip, 2) + center(2)
      if (iy >= min_y .and. iy < max_y) then
        iz = mesh%idx%lxyz(ip, 3) + center(3)
        do ii = 0, nn - 1
          if (iz+ii >= min_z .and. iz+ii < max_z) then      
            if (cf%pfft%is_real == 0) then
              cf%pfft%data_in(cube_get_pfft_index(cf%pfft,ix,iy,iz+ii)) = TOCMPLX(real(gmf(ip + ii)),M_ZERO)
            else 
              cf%pfft%data_in(cube_get_pfft_index(cf%pfft,ix,iy,iz+ii)) = gmf(ip + ii)
            end if
          else
            cycle
          end if
        end do
      end if
    end if
  end do
  
  if(local_) then
    SAFE_DEALLOCATE_P(gmf)
  end if

  call profiling_count_transfers(mesh%np_global, mf(1))

  call profiling_out(prof_m2c)
  POP_SUB(X(mesh_to_cube_parallel))
end subroutine X(mesh_to_cube_parallel)

! ---------------------------------------------------------
subroutine X(cube_to_mesh_parallel) (mesh, cf, mf)
  type(mesh_t),  intent(in)  :: mesh
  type(cube_function_t), intent(in)  :: cf
  R_TYPE,        intent(out) :: mf(:)  ! mf(mesh%np_global)

  integer :: ip, ix, iy, iz, center(3)
  integer :: im, ii, nn
  integer :: min_x, min_y, min_z, max_x,max_y,max_z,recv_count
  FLOAT :: scaling_fft_factor

  PUSH_SUB(X(cube_to_mesh_parallel))

  call profiling_in(prof_c2m, "CUBE_TO_MESH_PARALLEL")

  center(1:3) = cf%n(1:3)/2 + 1

  min_x = cf%pfft%local_i_start(1)
  min_y = cf%pfft%local_i_start(2)
  min_z = cf%pfft%local_i_start(3)
  max_x = cf%pfft%local_i_start(1)+cf%pfft%local_ni(1)
  max_y = cf%pfft%local_i_start(2)+cf%pfft%local_ni(2)
  max_z = cf%pfft%local_i_start(3)+cf%pfft%local_ni(3)
       
  scaling_fft_factor = real(cf%pfft%n(1)*cf%pfft%n(2)*cf%pfft%n(3))
  if (mpi_world%size == 1) then
    do im = 1, mesh%cube_map%nmap
      write(*,*)im,"do",scaling_fft_factor
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      ix = mesh%idx%lxyz(ip, 1) + center(1)
      iy = mesh%idx%lxyz(ip, 2) + center(2)
      iz = mesh%idx%lxyz(ip, 3) + center(3)
      forall(ii = 0:nn - 1) mf(ip + ii) = cf%pfft%global_data_in(cf%pfft%get_local_index(ix,iy,iz+ii))/scaling_fft_factor
    end do
  else
    !General execution
    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      ix = mesh%idx%lxyz(ip, 1) + center(1)
      iy = mesh%idx%lxyz(ip, 2) + center(2)
      iz = mesh%idx%lxyz(ip, 3) + center(3)

      do ii = 0, nn - 1
        !Only work with the points of the mesh of my mesh
        if (mesh%vp%part(ip+ii) == mesh%vp%partno) then
          mf(ip + ii) = cf%pfft%global_data_in(cf%pfft%get_local_index(ix,iy,iz+ii))/scaling_fft_factor
        end if
      end do
    end do
  end if

  call profiling_count_transfers(mesh%np_global, mf(1))
  call profiling_out(prof_c2m)

  POP_SUB(X(cube_to_mesh_parallel))

end subroutine X(cube_to_mesh_parallel)

#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
