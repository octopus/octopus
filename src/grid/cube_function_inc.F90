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

  call push_sub('cf_inc.Xcf_new')

  ASSERT(all(n>0))

  nullify(cf%RS)
  nullify(cf%FS)
  cf%n = n

  nullify(cf%fft)
  call pop_sub()
end subroutine X(cf_new)


! ---------------------------------------------------------
subroutine X(cf_new_from)(cf, cf_i)
  type(X(cf_t)), intent(out) :: cf
  type(X(cf_t)), intent( in) :: cf_i

  call push_sub('cf_inc.Xcf_new_from')
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
  call pop_sub()
end subroutine X(cf_new_from)


! ---------------------------------------------------------
subroutine X(cf_alloc_RS)(cf)
  type(X(cf_t)), intent(inout) :: cf

  call push_sub('cf_inc.Xcf_alloc_RS')

  ASSERT(.not.associated(cf%RS))
  SAFE_ALLOCATE(cf%RS(1:cf%n(1), 1:cf%n(2), 1:cf%n(3)))

  call pop_sub()
end subroutine X(cf_alloc_RS)


! ---------------------------------------------------------
subroutine X(cf_free_RS)(cf)
  type(X(cf_t)), intent(inout) :: cf

  call push_sub('cf_inc.Xcf_free_RS')

  ASSERT(associated(cf%RS))
  SAFE_DEALLOCATE_P(cf%RS)
  nullify(cf%RS)

  call pop_sub()
end subroutine X(cf_free_RS)

! ---------------------------------------------------------
subroutine X(cf_free)(cf)
  type(X(cf_t)), intent(inout) :: cf

  call push_sub('cf_inc.Xcf_free')

  if(associated(cf%RS)) then
    SAFE_DEALLOCATE_P(cf%RS)
    nullify(cf%RS)
  end if

  if(associated(cf%FS)) then
    SAFE_DEALLOCATE_P(cf%FS)
    nullify(cf%FS)
  end if

  if(associated(cf%fft)) then
    call fft_end(cf%fft)
    SAFE_DEALLOCATE_P(cf%fft)
    nullify(cf%fft)
  end if

  call pop_sub()
end subroutine X(cf_free)

! ---------------------------------------------------------
! The next two subroutines convert a function between the normal
! mesh and the cube.
! Note that the function in the mesh should be defined
! globally, not just in a partition (when running in
! parallel in real-space domains).
! ---------------------------------------------------------

subroutine X(mesh_to_cube) (m, mf, cf)
  type(mesh_t),  intent(in)    :: m
  R_TYPE,        intent(in)    :: mf(:)  ! mf(m%np_global)
  type(X(cf_t)), intent(inout) :: cf

  integer :: i, ix, iy, iz, c(3)

  call push_sub('cf_inc.Xmesh_to_cube')

  ASSERT(associated(cf%RS))

  c(1:3)  =  cf%n(1:3)/2 + 1

  cf%RS =  M_ZERO

  do i = 1, m%np_global
    ix = m%idx%Lxyz(i, 1) + c(1)
    iy = m%idx%Lxyz(i, 2) + c(2)
    iz = m%idx%Lxyz(i, 3) + c(3)

    cf%RS(ix, iy, iz) = mf(i)
  end do

  call pop_sub()
end subroutine X(mesh_to_cube)

! ---------------------------------------------------------
subroutine X(cube_to_mesh) (m, cf, mf)
  type(mesh_t),  intent(in)  :: m
  type(X(cf_t)), intent(in)  :: cf
  R_TYPE,        intent(out) :: mf(:)  ! mf(m%np_global)

  integer :: i, ix, iy, iz, c(3)

  call push_sub('cf_inc.Xcube_to_mesh')

  ASSERT(associated(cf%RS))

  c(1:3) =  cf%n(1:3)/2 + 1

  do i = 1, m%np_global
    ix = m%idx%Lxyz(i, 1) + c(1)
    iy = m%idx%Lxyz(i, 2) + c(2)
    iz = m%idx%Lxyz(i, 3) + c(3)
    mf(i) = cf%RS(ix, iy, iz)
  end do

  call pop_sub()

end subroutine X(cube_to_mesh)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
