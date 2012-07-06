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
!> Allocates locally the real space grid, if PFFT library is not used.
!! Otherwise, it assigns the PFFT real space grid to the cube real space grid,
!! via pointer.
subroutine X(cube_function_alloc_rs)(cube, cf, in_device, force_alloc)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf
  logical, optional,     intent(in)    :: in_device
  logical, optional,     intent(in)    :: force_alloc

  logical :: allocated
  
  

  PUSH_SUB(X(cube_function_alloc_rs))

  ASSERT(.not.associated(cf%X(rs)))

  allocated = .false.

  cf%forced_alloc = optional_default(force_alloc, .false.)

  if(associated(cube%fft)) then
    select case(cube%fft%library)
    case(FFTLIB_PFFT)

      ASSERT(associated(cube%fft))
      if(.not. cf%forced_alloc) then  
        allocated = .true.
        cf%X(rs) => cube%fft%X(rs_data)(1:cube%rs_n(1), 1:cube%rs_n(2), 1:cube%rs_n(3))
      end if
    case(FFTLIB_CLAMD)
      if(optional_default(in_device, .true.)) then
        allocated = .true.
        cf%in_device_memory = .true.
#ifdef HAVE_OPENCL
        call opencl_create_buffer(cf%real_space_buffer, CL_MEM_READ_WRITE, TYPE_CMPLX, product(cube%rs_n(1:3)))
#endif
      end if
    end select
  end if

  if(.not. allocated) then
    SAFE_ALLOCATE(cf%X(rs)(1:cube%rs_n(1), 1:cube%rs_n(2), 1:cube%rs_n(3)))
  end if

  POP_SUB(X(cube_function_alloc_rs))
end subroutine X(cube_function_alloc_rs)


! ---------------------------------------------------------
!> Deallocates the real space grid
subroutine X(cube_function_free_rs)(cube, cf)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf

  logical :: deallocated

  PUSH_SUB(X(cube_function_free_rs))

  deallocated = .false.

  if(associated(cube%fft)) then
    select case(cube%fft%library)
    case(FFTLIB_PFFT)
      if(.not. cf%forced_alloc) then
        deallocated = .true.
        nullify(cf%X(rs))
      end if
    case(FFTLIB_CLAMD)
#ifdef HAVE_OPENCL
      if(cf%in_device_memory) then
        deallocated = .true.
        ASSERT(cf%in_device_memory)
        call opencl_release_buffer(cf%real_space_buffer)
        cf%in_device_memory = .false.
      end if
#endif
    end select
  end if

  if(.not. deallocated) then
    SAFE_DEALLOCATE_P(cf%X(rs))
  end if

  POP_SUB(X(cube_function_free_rs))
end subroutine X(cube_function_free_rs)

! ---------------------------------------------------------
subroutine X(cube_function_allgather)(cube, cf, cf_local)
  type(cube_t),   intent(in) :: cube
  R_TYPE,         intent(out) :: cf(:,:,:)
  R_TYPE,         intent(in)  :: cf_local(:,:,:)

  integer :: ix, iy, iz, index
  R_TYPE, allocatable :: cf_tmp(:)
  type(profile_t), save :: prof_allgather

  #ifdef HAVE_MPI

  PUSH_SUB(X(cube_function_allgather))
  call profiling_in(prof_allgather, "CF_ALLGATHER")

  SAFE_ALLOCATE(cf_tmp(cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3)))

  call mpi_debug_in(cube%mpi_grp%comm, C_MPI_ALLGATHERV)
  ! Warning: in the next line we have to pass the full cf_local array, not just the first element.
  ! This is because cf_local might be a pointer to a subarray when using PFFT, such that
  ! memory will not be contiguous (see cube_function_alloc_rs). In that case the 
  ! Fortran compiler should do a temporary copy.
  call MPI_Allgatherv ( cf_local, cube%np_local(cube%mpi_grp%rank+1), R_MPITYPE, &
       cf_tmp(1), cube%np_local, cube%xlocal - 1, R_MPITYPE, &
       cube%mpi_grp%comm, mpi_err)
  call mpi_debug_out(cube%mpi_grp%comm, C_MPI_ALLGATHERV)

  ! Copy values to cf in the correct order
  do index = 1, cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3)
    ix = cube%local(index, 1)
    iy = cube%local(index, 2)
    iz = cube%local(index, 3)
    cf(ix, iy, iz) = cf_tmp(index)
  end do

  SAFE_DEALLOCATE_A(cf_tmp)

  call profiling_out(prof_allgather)
#endif
  POP_SUB(X(cube_function_allgather))
end subroutine X(cube_function_allgather)


! ---------------------------------------------------------
!> The next two subroutines convert a function between the normal
!! mesh and the cube.
!!
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

  integer :: ip, ix, iy, iz
  integer :: im, ii, nn, bsize
  logical :: local_
  R_TYPE, pointer :: gmf(:)
#ifdef HAVE_OPENCL
  type(opencl_mem_t)         :: mf_buffer
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)            :: kernel_ref
#endif

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

  if(.not. cf%in_device_memory) then


    ASSERT(associated(cf%X(rs)))

    cf%X(rs) = M_ZERO

    ASSERT(associated(mesh%cube_map%map))
    ASSERT(mesh%sb%dim <= 3)

    do im = 1, mesh%cube_map%nmap
      ix = mesh%cube_map%map(1, im) + cube%center(1)
      iy = mesh%cube_map%map(2, im) + cube%center(2)
      iz = mesh%cube_map%map(3, im) + cube%center(3)

      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)
      forall(ii = 0:nn - 1) cf%X(rs)(ix, iy, iz + ii) = gmf(ip + ii)
    end do

  else
#ifdef HAVE_OPENCL

    call opencl_set_buffer_to_zero(cf%real_space_buffer, TYPE_CMPLX, product(cube%rs_n(1:3)))

    call opencl_create_buffer(mf_buffer, CL_MEM_READ_ONLY, R_TYPE_VAL, mesh%np_global)
    call opencl_write_buffer(mf_buffer, mesh%np_global, gmf)

    call octcl_kernel_start_call(kernel, 'mesh_to_cube.cl', TOSTRING(X(mesh_to_cube)))
    kernel_ref = octcl_kernel_get_ref(kernel)
    
    call opencl_set_kernel_arg(kernel_ref, 0, mesh%cube_map%nmap)
    call opencl_set_kernel_arg(kernel_ref, 1, cube%rs_n(3)*cube%rs_n(2))
    call opencl_set_kernel_arg(kernel_ref, 2, cube%rs_n(3))
    call opencl_set_kernel_arg(kernel_ref, 3, 1)
    call opencl_set_kernel_arg(kernel_ref, 4, cube%center(1))
    call opencl_set_kernel_arg(kernel_ref, 5, cube%center(2))
    call opencl_set_kernel_arg(kernel_ref, 6, cube%center(3))
    call opencl_set_kernel_arg(kernel_ref, 7, mesh%cube_map%map_buffer)
    call opencl_set_kernel_arg(kernel_ref, 8, mf_buffer)
    call opencl_set_kernel_arg(kernel_ref, 9, cf%real_space_buffer)

    bsize = opencl_kernel_workgroup_size(kernel_ref)

    call opencl_kernel_run(kernel_ref, (/pad(mesh%cube_map%nmap, bsize)/), (/bsize/))
    call opencl_finish()
    call opencl_release_buffer(mf_buffer)

#endif
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

  integer :: ip, ix, iy, iz
  integer :: im, ii, nn
#ifdef HAVE_MPI
  integer :: first, last
#endif
  logical :: local_
  R_TYPE, pointer :: gmf(:)
#ifdef HAVE_OPENCL
  integer                    :: bsize
  type(opencl_mem_t)         :: mf_buffer
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)            :: kernel_ref
#endif

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

  if(.not. cf%in_device_memory) then

    ASSERT(associated(cf%X(rs)))
    ASSERT(associated(mesh%cube_map%map))

    do im = 1, mesh%cube_map%nmap
      ix = mesh%cube_map%map(1, im) + cube%center(1)
      iy = mesh%cube_map%map(2, im) + cube%center(2)
      iz = mesh%cube_map%map(3, im) + cube%center(3)

      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      forall(ii = 0:nn - 1) gmf(ip + ii) = cf%X(rs)(ix, iy, iz + ii)
    end do

  else

#ifdef HAVE_OPENCL    
    call opencl_create_buffer(mf_buffer, CL_MEM_READ_ONLY, R_TYPE_VAL, mesh%np_global)

    call octcl_kernel_start_call(kernel, 'mesh_to_cube.cl', TOSTRING(X(cube_to_mesh)))
    kernel_ref = octcl_kernel_get_ref(kernel)
    
    call opencl_set_kernel_arg(kernel_ref, 0, mesh%cube_map%nmap)
    call opencl_set_kernel_arg(kernel_ref, 1, cube%rs_n(3)*cube%rs_n(2))
    call opencl_set_kernel_arg(kernel_ref, 2, cube%rs_n(3))
    call opencl_set_kernel_arg(kernel_ref, 3, 1)
    call opencl_set_kernel_arg(kernel_ref, 4, cube%center(1))
    call opencl_set_kernel_arg(kernel_ref, 5, cube%center(2))
    call opencl_set_kernel_arg(kernel_ref, 6, cube%center(3))
    call opencl_set_kernel_arg(kernel_ref, 7, mesh%cube_map%map_buffer)
    call opencl_set_kernel_arg(kernel_ref, 8, cf%real_space_buffer)
    call opencl_set_kernel_arg(kernel_ref, 9, mf_buffer)

    bsize = opencl_kernel_workgroup_size(kernel_ref)

    call opencl_kernel_run(kernel_ref, (/pad(mesh%cube_map%nmap, bsize)/), (/bsize/))
    call opencl_finish()

    call opencl_read_buffer(mf_buffer, mesh%np_global, gmf)
    call opencl_release_buffer(mf_buffer)
#endif

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
! ---------------------------------------------------------
subroutine X(mesh_to_cube_parallel)(mesh, mf, cube, cf, map)
  type(mesh_t),          intent(in)    :: mesh
  R_TYPE,  target,       intent(in)    :: mf(:)  !< mf(mesh%np)
  type(cube_t),          intent(in)    :: cube
  type(cube_function_t), intent(inout) :: cf
  type(mesh_cube_parallel_map_t), intent(in) :: map

  integer :: ip, ix, iy, iz
  integer :: im, ii, nn
  integer :: min_x, min_y, min_z, max_x, max_y, max_z
  R_TYPE, allocatable :: in(:), out(:)

  PUSH_SUB(X(mesh_to_cube_parallel))
  call profiling_in(prof_m2c, "MESH_TO_CUBE_PARALLEL")

  ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)

  if (mesh%parallel_in_domains) then
    SAFE_ALLOCATE(in(map%m2c_nsend))
    SAFE_ALLOCATE(out(map%m2c_nrec))

    !Put the mesh function data in correct order for transfer
    do ip = 1, map%m2c_nsend
      in(ip) = mf(map%m2c_mf_order(ip))
    end do

    !Transfer the mesh function from the mesh partition to the cube partition
    call X(partition_transfer)(map%m2c, in, out)

    !Put the transfered mesh function in the cube
    cf%X(rs) = R_TOTYPE(M_ZERO)
    do ip = 1, map%m2c_nrec
      ix = map%m2c_cf_order(ip, 1)
      iy = map%m2c_cf_order(ip, 2)
      iz = map%m2c_cf_order(ip, 3)
      cf%X(rs)(ix, iy, iz) = out(ip)
    end do

    SAFE_DEALLOCATE_A(in)
    SAFE_DEALLOCATE_A(out)

  else
    ! Save the limit values
    min_x = cube%rs_istart(1)
    min_y = cube%rs_istart(2)
    min_z = cube%rs_istart(3)
    max_x = cube%rs_istart(1) + cube%rs_n(1)
    max_y = cube%rs_istart(2) + cube%rs_n(2)
    max_z = cube%rs_istart(3) + cube%rs_n(3)
  
    ! Initialize to zero the input matrix
    forall(iz = 1:cube%rs_n(3), iy = 1:cube%rs_n(2), ix = 1:cube%rs_n(1)) cf%X(rs)(ix, iy, iz) = M_ZERO

    ! Do the actual transform, only for the output values
    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)
    
      ix = mesh%cube_map%map(1, im) + cube%center(1)
      if (ix >= min_x .and. ix < max_x) then
        iy = mesh%cube_map%map(2, im) + cube%center(2)
        if (iy >= min_y .and. iy < max_y) then
          iz = mesh%cube_map%map(3, im) + cube%center(3)
          do ii = 0, nn - 1
            if (iz+ii >= min_z .and. iz+ii < max_z) then
              cf%X(rs)(ix-min_x+1, iy-min_y+1, iz+ii-min_z+1) = mf(ip + ii)
            end if
          end do
        end if
      end if

    end do

  end if


  call profiling_count_transfers(mesh%np_global, mf(1))
  call profiling_out(prof_m2c)

  POP_SUB(X(mesh_to_cube_parallel))
end subroutine X(mesh_to_cube_parallel)

! ---------------------------------------------------------
subroutine X(cube_to_mesh_parallel) (cube, cf, mesh, mf, map)
  type(cube_t),          intent(in)  :: cube
  type(cube_function_t), intent(in)  :: cf
  type(mesh_t),          intent(in)  :: mesh
  R_TYPE,                intent(out) :: mf(:)  !< mf(mesh%np)
  type(mesh_cube_parallel_map_t), intent(in) :: map

  integer :: ip, ix, iy, iz, im, ii, nn, ixyz(3)
  R_TYPE, allocatable :: gcf(:,:,:)
  R_TYPE, allocatable :: in(:), out(:)

  PUSH_SUB(X(cube_to_mesh_parallel))
  call profiling_in(prof_c2m, "CUBE_TO_MESH_PARALLEL")

  ASSERT(.not. cf%in_device_memory)
  ASSERT(ubound(mf, dim = 1) == mesh%np .or. ubound(mf, dim = 1) == mesh%np_part)

  if(mesh%parallel_in_domains) then
    SAFE_ALLOCATE(in(map%c2m_nsend))
    SAFE_ALLOCATE(out(map%c2m_nrec))

    !Put the cube function in the mesh in the correct order for transfer
    do ip = 1, map%c2m_nsend
      ix = map%c2m_cf_order(ip, 1)
      iy = map%c2m_cf_order(ip, 2)
      iz = map%c2m_cf_order(ip, 3)
      in(ip) = cf%X(rs)(ix, iy, iz)
    end do

    !Transfer the mesh function from the cube partition to the mesh partition
    call X(partition_transfer)(map%c2m, in, out)

    !Put the transfered mesh function data in correct order
    do ip = 1, map%c2m_nrec
      mf(map%c2m_mf_order(ip)) = out(ip)
    end do

    SAFE_DEALLOCATE_A(in)
    SAFE_DEALLOCATE_A(out)

  else
    !collect the data in all processes
    SAFE_ALLOCATE(gcf(cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3)))
#ifdef HAVE_MPI
    call X(cube_function_allgather)(cube, gcf, cf%X(rs))
#endif

    ! cube to mesh
    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      ixyz(1:3) = mesh%cube_map%map(1:3, im) + cube%center(1:3)

      forall(ii = 0:nn - 1) mf(ip + ii) = gcf(ixyz(1), ixyz(2), ixyz(3) + ii)
    end do

    SAFE_DEALLOCATE_A(gcf)
  end if

  call profiling_count_transfers(mesh%np_global, mf(1))
  call profiling_out(prof_c2m)

  POP_SUB(X(cube_to_mesh_parallel))
end subroutine X(cube_to_mesh_parallel)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
