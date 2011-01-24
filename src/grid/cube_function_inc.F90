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
! The following routines handle creation/destruction of the cube
subroutine X(cf_new)(n, cf)
  integer, intent(in) :: n(3)
  type(X(cf_t)), intent(out) :: cf

  PUSH_SUB(X(cf_new))

  ASSERT(all(n>0))

  nullify(cf%RS)
  nullify(cf%FS)
  cf%n = n

  nullify(cf%fft)
  POP_SUB(X(cf_new))
end subroutine X(cf_new)


! ---------------------------------------------------------
subroutine X(cf_new_from)(cf, cf_i)
  type(X(cf_t)), intent(out) :: cf
  type(X(cf_t)), intent( in) :: cf_i

  PUSH_SUB(X(cf_new_from))
  ASSERT(all(cf_i%n>0))

  nullify(cf%RS)
  nullify(cf%FS)
  cf%n = cf_i%n

  if(associated(cf_i%fft)) then
    SAFE_ALLOCATE(cf%fft)
    call fft_copy(cf_i%fft, cf%fft)
    cf%nx = cf_i%nx
  else
    nullify(cf%fft)
  end if
  POP_SUB(X(cf_new_from))
end subroutine X(cf_new_from)


! ---------------------------------------------------------
subroutine X(cf_alloc_RS)(cf)
  type(X(cf_t)), intent(inout) :: cf

  PUSH_SUB(X(cf_alloc_RS))

  ASSERT(.not.associated(cf%RS))
  SAFE_ALLOCATE(cf%RS(1:cf%n(1), 1:cf%n(2), 1:cf%n(3)))

  POP_SUB(X(cf_alloc_RS))
end subroutine X(cf_alloc_RS)


! ---------------------------------------------------------
subroutine X(cf_free_RS)(cf)
  type(X(cf_t)), intent(inout) :: cf

  PUSH_SUB(X(cf_free_RS))

  ASSERT(associated(cf%RS))
  SAFE_DEALLOCATE_P(cf%RS)

  POP_SUB(X(cf_free_RS))
end subroutine X(cf_free_RS)

! ---------------------------------------------------------
subroutine X(cf_free)(cf)
  type(X(cf_t)), intent(inout) :: cf

  PUSH_SUB(X(cf_free))

  SAFE_DEALLOCATE_P(cf%RS)
  SAFE_DEALLOCATE_P(cf%FS)

  if(associated(cf%fft)) then
    call fft_end(cf%fft)
    SAFE_DEALLOCATE_P(cf%fft)
  end if

  POP_SUB(X(cf_free))
end subroutine X(cf_free)

! ---------------------------------------------------------
! The next two subroutines convert a function between the normal
! mesh and the cube.
! Note that the function in the mesh should be defined
! globally, not just in a partition (when running in
! parallel in real-space domains).
! ---------------------------------------------------------

subroutine X(mesh_to_cube) (mesh, mf, cf)
  type(mesh_t),  intent(in)    :: mesh
  R_TYPE,        intent(in)    :: mf(:)  ! mf(mesh%np_global)
  type(X(cf_t)), intent(inout) :: cf

  integer :: ip, ix, iy, iz, center(3)
  integer :: im, ii, nn

  PUSH_SUB(X(mesh_to_cube))
  call profiling_in(prof_m2c, "MESH_TO_CUBE")

  ASSERT(associated(cf%RS))

  center(1:3) = cf%n(1:3)/2 + 1

  cf%RS = M_ZERO

  do im = 1, mesh%cube_map%nmap
    ip = mesh%cube_map%map(POINT, im)
    nn = mesh%cube_map%map(COUNT, im)

    ix = mesh%idx%Lxyz(ip, 1) + center(1)
    iy = mesh%idx%Lxyz(ip, 2) + center(2)
    iz = mesh%idx%Lxyz(ip, 3) + center(3)
    forall(ii = 0:nn - 1) cf%RS(ix, iy, iz + ii) = mf(ip + ii)
  end do
  
  call profiling_count_transfers(mesh%np_global, mf(1))

  call profiling_out(prof_m2c)
  POP_SUB(X(mesh_to_cube))
end subroutine X(mesh_to_cube)

! ---------------------------------------------------------
subroutine X(cube_to_mesh) (mesh, cf, mf)
  type(mesh_t),  intent(in)  :: mesh
  type(X(cf_t)), intent(in)  :: cf
  R_TYPE,        intent(out) :: mf(:)  ! mf(mesh%np_global)

  integer :: ip, ix, iy, iz, center(3)
  integer :: im, ii, nn

  PUSH_SUB(X(cube_to_mesh))

  call profiling_in(prof_c2m, "CUBE_TO_MESH")

  ASSERT(associated(cf%RS))

  center(1:3) = cf%n(1:3)/2 + 1

  do im = 1, mesh%cube_map%nmap
    ip = mesh%cube_map%map(POINT, im)
    nn = mesh%cube_map%map(COUNT, im)
    ix = mesh%idx%Lxyz(ip, 1) + center(1)
    iy = mesh%idx%Lxyz(ip, 2) + center(2)
    iz = mesh%idx%Lxyz(ip, 3) + center(3)
    forall(ii = 0:nn - 1) mf(ip + ii) = cf%RS(ix, iy, iz + ii)
  end do
  
  call profiling_count_transfers(mesh%np_global, mf(1))

  call profiling_out(prof_c2m)

  POP_SUB(X(cube_to_mesh))

end subroutine X(cube_to_mesh)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
