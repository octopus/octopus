!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
subroutine X(cube_function_alloc_RS)(cube, cf)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_alloc_RS))

  !Save memory if PFFT is used  
  if (cube%fft_library /= PFFT_LIB) then
    ASSERT(.not.associated(cf%X(RS)))
    SAFE_ALLOCATE(cf%X(RS)(1:cube%n(1), 1:cube%n(2), 1:cube%n(3)))
  end if

  POP_SUB(X(cube_function_alloc_RS))
end subroutine X(cube_function_alloc_RS)


! ---------------------------------------------------------
subroutine X(cube_function_free_RS)(cf)
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_free_RS))

  if (associated(cf%X(RS))) then
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

subroutine X(mesh_to_cube)(mesh, mf, cube, cf, local)
  type(mesh_t),          intent(in)    :: mesh
  R_TYPE,  target,       intent(in)    :: mf(:) !< function defined on the mesh. Can be 
                                                !< mf(mesh%np) or mf(mesh%np_global), depending if it is a local or global function
  type(cube_t),          intent(in)    :: cube
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

  center(1:3) = cube%n(1:3)/2 + 1

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
subroutine X(cube_to_mesh) (cube, cf, mesh, mf, local)
  type(cube_t),          intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  type(mesh_t),          intent(in)  :: mesh
  R_TYPE,  target,       intent(out) :: mf(:) !< function defined on the mesh. Can be 
                                              !< mf(mesh%np) or mf(mesh%np_global), depending if it is a local or global function
  logical, optional,     intent(in)  :: local  !< If .true. the mf array is a local array. Considered .false. if not present.

  integer :: ip, ix, iy, iz, center(3), ixyz(1:3)
  integer :: im, ii, nn, last, first
  logical :: local_
  R_TYPE, pointer :: gmf(:)

  PUSH_SUB(X(cube_to_mesh))

  call profiling_in(prof_c2m, "CUBE_TO_MESH")

  local_ = optional_default(local, .false.) .and. mesh%parallel_in_domains
  
  if(local_) then
    ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)
    SAFE_ALLOCATE(gmf(1:mesh%np_global))
  else
    ASSERT(ubound(mf, dim = 1) == mesh%np_global .or. ubound(mf, dim = 1) == mesh%np_part_global)
    gmf => mf
  end if

  ASSERT(associated(cf%X(RS)))

  center(1:3) = cube%n(1:3)/2 + 1

  ASSERT(associated(mesh%cube_map%map))

  if(associated(mesh%idx%lxyz)) then

    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)
      ix = mesh%idx%lxyz(ip, 1) + center(1)
      iy = mesh%idx%lxyz(ip, 2) + center(2)
      iz = mesh%idx%lxyz(ip, 3) + center(3)
      forall(ii = 0:nn - 1) gmf(ip + ii) = cf%X(RS)(ix, iy, iz + ii)
    end do

  else

    ixyz = 0

    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      call index_to_coords(mesh%idx, mesh%sb%dim, ip, ixyz)
      ixyz = ixyz + center

      forall(ii = 0:nn - 1) gmf(ip + ii) = cf%X(RS)(ixyz(1), ixyz(2), ixyz(3) + ii)
    end do

  end if

  if(local_) then
#ifdef HAVE_MPI
    first = mesh%vp%xlocal(mesh%vp%partno)
    last = mesh%vp%np_local(mesh%vp%partno)   
    do ii = 1, last
      mf(ii) = gmf(mesh%vp%local(ii+first-1))
    end do
#endif
    SAFE_DEALLOCATE_P(gmf)
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
subroutine X(mesh_to_cube_parallel)(mesh, mf, cube, cf, local)
  type(mesh_t),          intent(in)    :: mesh
  R_TYPE,  target,       intent(in)    :: mf(:)  !< mf(mesh%np_global)
  type(cube_t),          intent(in)    :: cube
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

  center(1:3) = cube%n(1:3)/2 + 1
  

  ! Save the limit values
  min_x = cube%pfft%local_i_start(1)
  min_y = cube%pfft%local_i_start(2)
  min_z = cube%pfft%local_i_start(3)
  max_x = cube%pfft%local_i_start(1) + cube%pfft%local_ni(1)
  max_y = cube%pfft%local_i_start(2) + cube%pfft%local_ni(2)
  max_z = cube%pfft%local_i_start(3) + cube%pfft%local_ni(3)
  
  ! Initialize to zero the input matrix
  do ix = min_x, max_x - 1
    do iy = min_y, max_y - 1
      do iz = min_z, max_z - 1
        cube%pfft%data_in(cube_get_pfft_index(cube%pfft, ix, iy, iz)) = M_ZERO
      end do
    end do
  end do

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
            if (cube%pfft%is_real == 0) then
              cube%pfft%data_in(cube_get_pfft_index(cube%pfft,ix,iy,iz+ii)) = TOCMPLX(real(gmf(ip + ii)),M_ZERO)
            else 
              cube%pfft%data_in(cube_get_pfft_index(cube%pfft,ix,iy,iz+ii)) = gmf(ip + ii)
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
subroutine X(cube_to_mesh_parallel) (cube, cf, mesh, mf, local)
  type(cube_t),          intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  type(mesh_t),          intent(in)  :: mesh
  R_TYPE, target,        intent(out) :: mf(:)  ! mf(mesh%np_global)
  logical, optional,     intent(in)  :: local  !< If .true. the mf array is a local array. Considered .false. if not present.

  integer :: ip, ix, iy, iz, center(3)
  integer :: im, ii, nn, last, first
  integer :: min_x, min_y, min_z, max_x,max_y,max_z,recv_count
  FLOAT :: scaling_fft_factor
  logical :: local_
  R_TYPE, pointer :: gmf(:)

  PUSH_SUB(X(cube_to_mesh_parallel))

  call profiling_in(prof_c2m, "CUBE_TO_MESH_PARALLEL")

  local_ = optional_default(local, .false.) .and. mesh%parallel_in_domains
  
  if(local_) then
    ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)
    SAFE_ALLOCATE(gmf(1:mesh%np_global))
  else
    ASSERT(ubound(mf, dim = 1) == mesh%np_global .or. ubound(mf, dim = 1) == mesh%np_part_global)
    gmf => mf
  end if

  center(1:3) = cube%n(1:3)/2 + 1

  min_x = cube%pfft%local_i_start(1)
  min_y = cube%pfft%local_i_start(2)
  min_z = cube%pfft%local_i_start(3)
  max_x = cube%pfft%local_i_start(1) + cube%pfft%local_ni(1)
  max_y = cube%pfft%local_i_start(2) + cube%pfft%local_ni(2)
  max_z = cube%pfft%local_i_start(3) + cube%pfft%local_ni(3)
       
  scaling_fft_factor = real(cube%pfft%n(1)*cube%pfft%n(2)*cube%pfft%n(3))
  if (mpi_world%size == 1) then
    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      ix = mesh%idx%lxyz(ip, 1) + center(1)
      iy = mesh%idx%lxyz(ip, 2) + center(2)
      iz = mesh%idx%lxyz(ip, 3) + center(3)
      forall(ii = 0:nn - 1) gmf(ip + ii) = cube%pfft%global_data_in(cube%pfft%get_local_index(ix,iy,iz+ii))/scaling_fft_factor
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
          gmf(ip + ii) = cube%pfft%global_data_in(cube%pfft%get_local_index(ix,iy,iz+ii))/scaling_fft_factor
        end if
      end do
    end do
  end if

  if(local_) then
#ifdef HAVE_MPI
    first = mesh%vp%xlocal(mesh%vp%partno)
    last = mesh%vp%np_local(mesh%vp%partno)   
    do ii = 1, last
      mf(ii) = gmf(mesh%vp%local(ii+first-1))
    end do
#endif
    SAFE_DEALLOCATE_P(gmf)
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
