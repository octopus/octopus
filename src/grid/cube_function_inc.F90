!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira, J. Alberdi, P. Garcia Risueño
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

  if (cube%fft_library /= FFTLIB_PFFT) then
    ASSERT(.not.associated(cf%X(RS)))
    SAFE_ALLOCATE(cf%X(RS)(1:cube%n(1), 1:cube%n(2), 1:cube%n(3)))
#ifdef HAVE_PFFT
  else
    ASSERT(.not.associated(cf%pRS))
    cf%pRS => cube%pfft%rs_data
#endif
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
#ifdef HAVE_PFFT
  if (associated(cf%pRS)) then
    nullify(cf%pRS)
  end if
#endif

  POP_SUB(X(cube_function_free_RS))
end subroutine X(cube_function_free_RS)

! ---------------------------------------------------------
#ifdef HAVE_MPI
subroutine X(cube_function_allgather)(cube, cf, cf_local)
  type(cube_t),   intent(in) :: cube
  R_TYPE,         intent(out) :: cf(:,:,:)
  R_TYPE,         intent(in)  :: cf_local(:,:,:)

  integer :: ix, iy, iz, index
  R_TYPE, allocatable :: cf_tmp(:)
  type(profile_t), save :: prof_allgather

  PUSH_SUB(X(cube_function_allgather))
  call profiling_in(prof_allgather, "CF_ALLGATHER")

  SAFE_ALLOCATE(cf_tmp(cube%n(1)*cube%n(2)*cube%n(3)))

  call mpi_debug_in(cube%mpi_grp%comm, C_MPI_ALLGATHERV)
  call MPI_Allgatherv ( cf_local(1,1,1), cube%np_local(cube%mpi_grp%rank+1), R_MPITYPE, &
       cf_tmp(1), cube%np_local, cube%xlocal - 1, R_MPITYPE, &
       cube%mpi_grp%comm, mpi_err)
  call mpi_debug_out(cube%mpi_grp%comm, C_MPI_ALLGATHERV)

  ! Copy values to cf in the correct order
  do index = 1, cube%n(1)*cube%n(2)*cube%n(3)
    ix = cube%local(index, 1)
    iy = cube%local(index, 2)
    iz = cube%local(index, 3)
    cf(ix, iy, iz) = cf_tmp(index)
  end do

  SAFE_DEALLOCATE_A(cf_tmp)

  call profiling_out(prof_allgather)

  POP_SUB(X(cube_function_allgather))
end subroutine X(cube_function_allgather)
#endif

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

! ---------------------------------------------------------
!> The next two subroutines convert a function between the normal
!! mesh and the cube in parallel.
!! ``Note that the function in the mesh should be defined
!! globally, not just in a partition (when running in
!! parallel in real-space domains).``
! ---------------------------------------------------------
subroutine X(mesh_to_cube_parallel)(mesh, mf, cube, cf, local, pfft_part)
  type(mesh_t),          intent(in)    :: mesh
  R_TYPE,  target,       intent(in)    :: mf(:)  !< mf(mesh%np_global)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf
  logical, optional,     intent(in)    :: local  !< If .true. the mf array is a local array. Considered .false. if not present.
  logical, optional,     intent(in)    :: pfft_part !< If .true. the used partition is the pfft equal partition

  integer :: im, ip, ii, nn, index, center(3), ixyz(3), lxyz(3), ix, iy, iz
  logical :: local_, pfft_part_
  R_TYPE, pointer :: gmf(:)

  PUSH_SUB(X(mesh_to_cube_parallel))
  call profiling_in(prof_m2c, "MESH_TO_CUBE_PARALLEL")

  local_ = optional_default(local, .false.) .and. mesh%parallel_in_domains
  pfft_part_ = optional_default(pfft_part, .false.) 
  if (pfft_part_) then
    ASSERT(ubound(mf, dim=1) == mesh%np)
  else
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
  end if
  
  ! Initialize to zero the input matrix
  forall(iz = 1:cube%rs_n(3), iy = 1:cube%rs_n(2), ix = 1:cube%rs_n(1)) cf%pRS(ix, iy, iz) = M_Z0

  ! Do the actual transform, only for the output values
  center(1:3) = cube%n(1:3)/2 + 1
  do im = 1, mesh%cube_map%nmap
    ip = mesh%cube_map%map(MCM_POINT, im)
    nn = mesh%cube_map%map(MCM_COUNT, im)

    ixyz(1:3) = mesh%idx%lxyz(ip, 1:3) + center(1:3)
    do ii = 0, nn - 1
      ixyz(1:3) = mesh%idx%lxyz(ip + ii, 1:3) + center(1:3)
      if (cube_global2local(cube, ixyz, lxyz)) then

        if (pfft_part_) then
#ifdef HAVE_MPI
#ifdef R_TCOMPLEX
          cf%pRS(lxyz(1), lxyz(2), lxyz(3)) = mf(vec_global2local(mesh%vp,ip+ii, mesh%vp%partno))
#else
          cf%pRS(lxyz(1), lxyz(2), lxyz(3)) = TOCMPLX(real(mf(vec_global2local(mesh%vp,ip+ii, mesh%vp%partno))),M_ZERO)
#endif
#endif
        else
#ifdef R_TCOMPLEX
          cf%pRS(lxyz(1), lxyz(2), lxyz(3)) = gmf(ip + ii)
#else
          cf%pRS(lxyz(1), lxyz(2), lxyz(3)) = TOCMPLX(real(gmf(ip + ii)),M_ZERO)
#endif
        end if
      end if
      
    end do

  end do
  
  if(local_) then
    SAFE_DEALLOCATE_P(gmf)
  end if

  call profiling_count_transfers(mesh%np_global, mf(1))

  call profiling_out(prof_m2c)
  POP_SUB(X(mesh_to_cube_parallel))
end subroutine X(mesh_to_cube_parallel)

! ---------------------------------------------------------
subroutine X(cube_to_mesh_parallel) (cube, cf, mesh, mf, local, pfft_part)
  type(cube_t),          intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  type(mesh_t),          intent(in)  :: mesh
  R_TYPE, target,        intent(out) :: mf(:)  ! mf(mesh%np_global)
  logical, optional,     intent(in)  :: local  !< If .true. the mf array is a local array. Considered .false. if not present.
  logical, optional,     intent(in)  :: pfft_part !< If .true. the used partition is the pfft equal partition

  integer :: ip, im, ii, nn, index, center(3), ixyz(3), lxyz(3)
  integer :: last, first
  logical :: local_, pfft_part_
  R_TYPE, pointer :: gmf(:)
  CMPLX, allocatable :: gcf(:,:,:)

  PUSH_SUB(X(cube_to_mesh_parallel))

  call profiling_in(prof_c2m, "CUBE_TO_MESH_PARALLEL")

  local_ = optional_default(local, .false.) .and. mesh%parallel_in_domains
  pfft_part_ = optional_default(pfft_part, .false.) 
  if (pfft_part_) then
    ASSERT(ubound(mf, dim=1) == mesh%np)
  else
    if(local_) then
      ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)
      SAFE_ALLOCATE(gmf(1:mesh%np_global))
    else
      ASSERT(ubound(mf, dim = 1) == mesh%np_global .or. ubound(mf, dim = 1) == mesh%np_part_global)
      gmf => mf
    end if
  end if
  
  center(1:3) = cube%n(1:3)/2 + 1

  if (pfft_part_) then

    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      do ii = 0, nn - 1
        ixyz(1:3) = mesh%idx%lxyz(ip + ii, 1:3) + center(1:3)

        if (cube_global2local(cube, ixyz, lxyz)) then
#ifdef HAVE_MPI
#ifdef R_TCOMPLEX
          mf(vec_global2local(mesh%vp,ip+ii, mesh%vp%partno)) = cf%pRS(lxyz(1), lxyz(2), lxyz(3))
#else
          mf(vec_global2local(mesh%vp,ip+ii, mesh%vp%partno)) = real(cf%pRS(lxyz(1), lxyz(2), lxyz(3)))
#endif
#endif
        end if
      end do
    end do

  else

    !collect the data in all processes
    SAFE_ALLOCATE(gcf(cube%n(1), cube%n(2), cube%n(3)))
#ifdef HAVE_MPI
    call zcube_function_allgather(cube, gcf, cf%pRS)
#endif

    ! cube to mesh
    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      ixyz(1:3) = mesh%idx%lxyz(ip, 1:3) + center(1:3)
      do ii = 0, nn - 1
        if (mpi_world%size == 1) then
#ifdef R_TCOMPLEX
          gmf(ip + ii) = gcf(ixyz(1), ixyz(2), ixyz(3) + ii)
#else
          gmf(ip + ii) = real(gcf(ixyz(1), ixyz(2), ixyz(3) + ii))
#endif
        else
          if (mesh%vp%part(ip+ii) == mesh%vp%partno) then
#ifdef R_TCOMPLEX
            gmf(ip + ii) = gcf(ixyz(1), ixyz(2), ixyz(3) + ii)
#else
            gmf(ip + ii) = real(gcf(ixyz(1), ixyz(2), ixyz(3) + ii))
#endif

          end if
        end if
      end do
    end do

    SAFE_DEALLOCATE_A(gcf)

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
  end if

  call profiling_count_transfers(mesh%np_global, mf(1))
  call profiling_out(prof_c2m)

  POP_SUB(X(cube_to_mesh_parallel))

end subroutine X(cube_to_mesh_parallel)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
