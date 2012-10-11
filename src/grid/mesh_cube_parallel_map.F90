!! Copyright (C) 2012 M. Oliveira
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

#include "global.h"

module mesh_cube_parallel_map_m
  use cube_m
  use datasets_m
  use global_m
  use index_m
  use mesh_m
  use mesh_cube_map_m
  use messages_m
  use mpi_m
  use partition_transfer_m
  use par_vec_m
  use profiling_m
  use simul_box_m

  implicit none

  private
  public ::                       &
    mesh_cube_parallel_map_t,     &
    mesh_cube_parallel_map_init,  &
    mesh_cube_parallel_map_end

  type mesh_cube_parallel_map_t
    ! Mesh to cube:
    type(partition_transfer_t) :: m2c
    integer :: m2c_nsend !< How many points will this process send to the cube partition
    integer :: m2c_nrec  !< How many points will this process receive from the mesh partition
    integer, pointer :: m2c_mf_order(:)   !< How the points of the mesh function should be ordered to execute 
                                          !! the transfer from the mesh partition to the cube partition. 
                                          !! These are *local* indexes.
    integer, pointer :: m2c_cf_order(:,:) !< How the points of the mesh function are ordered after executing 
                                          !! the transfer from the mesh partition to the cube partition with 
                                          !! respect to the *local* cube indexes.

    ! Cube to mesh:
    type(partition_transfer_t) :: c2m
    integer :: c2m_nsend !< How many points will this process send to the mesh partition
    integer :: c2m_nrec  !< How many points will this process receive from the cube partition
    integer, pointer :: c2m_cf_order(:,:) !< How the points of the mesh function should be ordered to execute
                                          !! the transfer from the cube partition to the mesh partition with 
                                          !! respect to the *local* cube indexes.
    integer, pointer :: c2m_mf_order(:)   !< How the points of the mesh function are ordered after executing 
                                          !! the transfer from the cube partition to the mesh partition. 
                                          !! These are *local* indexes.
  end type mesh_cube_parallel_map_t

contains

  ! ---------------------------------------------------------
  subroutine mesh_cube_parallel_map_init(this, mesh, cube)
    type(mesh_cube_parallel_map_t), intent(out) :: this
    type(mesh_t),                   intent(in)  :: mesh
    type(cube_t),                   intent(in)  :: cube

    integer :: im, ip, nn, ixyz(3), lxyz(3), ii
    integer, allocatable :: cube_part(:)
    integer, pointer :: mf_order(:), cf_order(:)
    type(dimensions_t), allocatable :: part(:)

    integer :: last_found_proc
    PUSH_SUB(mesh_cube_parallel_map_init)

    !Get the cube partition on the mesh
    SAFE_ALLOCATE(part(1:cube%mpi_grp%size))
    call cube_partition(cube, part)

    SAFE_ALLOCATE(cube_part(1:mesh%np_global))
    
    ixyz = 0
    last_found_proc = 1
    do im = 1, mesh%cube_map%nmap
      ip = mesh%cube_map%map(MCM_POINT, im)
      nn = mesh%cube_map%map(MCM_COUNT, im)

      call index_to_coords(mesh%idx, mesh%sb%dim, ip, ixyz)
      ixyz = ixyz + cube%center

      forall(ii = 0:nn - 1) cube_part(ip + ii) = cube_point_to_process(ixyz(1), ixyz(2), ixyz(3) + ii, part, last_found_proc)
    end do

    SAFE_DEALLOCATE_A(part)

    !Init mesh to cube partition transfer
    call partition_transfer_init(this%m2c, mesh%np_global, &
                                 mesh%mpi_grp, mesh%vp%part, cube%mpi_grp, cube_part, &
                                 this%m2c_nsend, this%m2c_nrec, mf_order, cf_order)

    SAFE_ALLOCATE(this%m2c_mf_order(1:this%m2c_nsend))
    SAFE_ALLOCATE(this%m2c_cf_order(1:this%m2c_nrec, 1:3))    
    do ip = 1, this%m2c_nsend
#ifdef HAVE_MPI
      this%m2c_mf_order(ip) = vec_global2local(mesh%vp, mf_order(ip), mesh%vp%partno)
#endif

      if (this%m2c_mf_order(ip) == 0) then
        write(message(1),'(a,i3,a,i3)') "Error in mesh_cube_parallel_map_init: m2c point ", &
             mf_order(ip), " is not stored in partition ", mesh%vp%partno
        call messages_fatal(1)
      end if

    end do
    do ip = 1, this%m2c_nrec
      call index_to_coords(mesh%idx, mesh%sb%dim, cf_order(ip), ixyz)
      ixyz = ixyz + cube%center

      if (.not. cube_global2local(cube, ixyz, lxyz)) then
        message(1) = "Error in mesh_cube_parallel_map_init"
        call messages_fatal(1)
      end if

      this%m2c_cf_order(ip, 1:3) = lxyz(1:3)
    end do
    SAFE_DEALLOCATE_P(mf_order)
    SAFE_DEALLOCATE_P(cf_order)

    !Init cube to mesh partition transfer
    call partition_transfer_init(this%c2m, mesh%np_global, &
                                 cube%mpi_grp, cube_part, mesh%mpi_grp, mesh%vp%part, &
                                 this%c2m_nsend, this%c2m_nrec, cf_order, mf_order)

    SAFE_ALLOCATE(this%c2m_cf_order(1:this%c2m_nsend, 1:3))
    SAFE_ALLOCATE(this%c2m_mf_order(1:this%c2m_nrec))
    do ip = 1, this%c2m_nrec
#ifdef HAVE_MPI
      this%c2m_mf_order(ip) = vec_global2local(mesh%vp, mf_order(ip), mesh%vp%partno)
#endif
      if (this%c2m_mf_order(ip) == 0) then
        write(message(1),'(a,i3,a,i3)') "Error in mesh_cube_parallel_map_init: c2m point ", &
             mf_order(ip), " is not stored in partition ", mesh%vp%partno
      end if

    end do
    do ip = 1, this%c2m_nsend
      call index_to_coords(mesh%idx, mesh%sb%dim, cf_order(ip), ixyz)
      ixyz = ixyz + cube%center

      if (.not. cube_global2local(cube, ixyz, lxyz)) then
        message(1) = "Error in mesh_cube_parallel_map_init"
        call messages_fatal(1)
      end if

      this%c2m_cf_order(ip, 1:3) = lxyz(1:3)
    end do
    SAFE_DEALLOCATE_P(mf_order)
    SAFE_DEALLOCATE_P(cf_order)

    SAFE_DEALLOCATE_A(cube_part)

    POP_SUB(mesh_cube_parallel_map_init)
  end subroutine mesh_cube_parallel_map_init

  ! ---------------------------------------------------------
  subroutine mesh_cube_parallel_map_end(this)
    type(mesh_cube_parallel_map_t), intent(inout) :: this

    PUSH_SUB(mesh_cube_parallel_map_end)

    SAFE_DEALLOCATE_P(this%m2c_mf_order)
    SAFE_DEALLOCATE_P(this%m2c_cf_order)
    SAFE_DEALLOCATE_P(this%c2m_mf_order)
    SAFE_DEALLOCATE_P(this%c2m_cf_order)
    call partition_transfer_end(this%m2c)
    call partition_transfer_end(this%c2m)

    POP_SUB(mesh_cube_parallel_map_end)
  end subroutine mesh_cube_parallel_map_end

end module mesh_cube_parallel_map_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
