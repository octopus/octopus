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

  ASSERT(.not.associated(cf%X(RS)))
  SAFE_ALLOCATE(cf%X(RS)(1:cf%n(1), 1:cf%n(2), 1:cf%n(3)))

  POP_SUB(X(cube_function_alloc_RS))
end subroutine X(cube_function_alloc_RS)


! ---------------------------------------------------------
subroutine X(cube_function_free_RS)(cf)
  type(cube_function_t), intent(inout) :: cf

  PUSH_SUB(X(cube_function_free_RS))

  ASSERT(associated(cf%X(RS)))
  SAFE_DEALLOCATE_P(cf%X(RS))

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
  R_TYPE,  target,       intent(in)    :: mf(:)  !< mf(mesh%np_global)
  type(cube_function_t), intent(inout) :: cf
  logical, optional,     intent(in)    :: local  !< If .true. the mf array is a local array. Considered .false. if not present.

  integer :: ip, ix, iy, iz, center(3)
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

  do im = 1, mesh%cube_map%nmap
    ip = mesh%cube_map%map(MCM_POINT, im)
    nn = mesh%cube_map%map(MCM_COUNT, im)

    ix = mesh%idx%lxyz(ip, 1) + center(1)
    iy = mesh%idx%lxyz(ip, 2) + center(2)
    iz = mesh%idx%lxyz(ip, 3) + center(3)
    forall(ii = 0:nn - 1) cf%X(RS)(ix, iy, iz + ii) = gmf(ip + ii)
  end do

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

  integer :: ip, ix, iy, iz, center(3)
  integer :: im, ii, nn

  PUSH_SUB(X(cube_to_mesh))

  call profiling_in(prof_c2m, "CUBE_TO_MESH")

  ASSERT(associated(cf%X(RS)))

  center(1:3) = cf%n(1:3)/2 + 1

  ASSERT(associated(mesh%cube_map%map))

  do im = 1, mesh%cube_map%nmap
    ip = mesh%cube_map%map(MCM_POINT, im)
    nn = mesh%cube_map%map(MCM_COUNT, im)
    ix = mesh%idx%lxyz(ip, 1) + center(1)
    iy = mesh%idx%lxyz(ip, 2) + center(2)
    iz = mesh%idx%lxyz(ip, 3) + center(3)
    forall(ii = 0:nn - 1) mf(ip + ii) = cf%X(RS)(ix, iy, iz + ii)
  end do
  
  call profiling_count_transfers(mesh%np_global, mf(1))

  call profiling_out(prof_c2m)

  POP_SUB(X(cube_to_mesh))

end subroutine X(cube_to_mesh)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
