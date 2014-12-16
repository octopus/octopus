!! Copyright (C) 2014 X. Andrade
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

R_TYPE function X(mesh_interpolation_evaluate)(this, values, position) result(interpolated_value)
  type(mesh_interpolation_t), target, intent(in)    :: this
  R_TYPE,                             intent(in)    :: values(:)
  FLOAT,                              intent(in)    :: position(:)

  type(mesh_t), pointer :: mesh
  integer :: nm(1:MAX_DIM)
  FLOAT :: xd(1:MAX_DIM), posrel(1:MAX_DIM)
  R_TYPE :: c00, c10, c01, c11, c0, c1

  PUSH_SUB(X(mesh_interpolation_evaluate))

  mesh => this%mesh

  ASSERT(mesh%sb%dim <= 3)

  posrel(1:mesh%sb%dim) = (position(1:mesh%sb%dim) - mesh%sb%box_offset(1:mesh%sb%dim))/mesh%spacing(1:mesh%sb%dim)

  nm(1:mesh%sb%dim) = floor(posrel(1:mesh%sb%dim))

  xd(1:mesh%sb%dim) = posrel(1:mesh%sb%dim) - nm(1:mesh%sb%dim)

  ASSERT(all(xd(1:mesh%sb%dim) >= CNST(0.0)))
  ASSERT(all(xd(1:mesh%sb%dim) <= CNST(1.0)))

  select case(mesh%sb%dim)
  case(3)

    ! trilinear interpolation : http://en.wikipedia.org/wiki/Trilinear_interpolation
    c00 = values(mesh%idx%lxyz_inv(0 + nm(1), 0 + nm(2), 0 + nm(3)))*(CNST(1.0) - xd(1)) + &
      values(mesh%idx%lxyz_inv(1 + nm(1), 0 + nm(2), 0 + nm(3)))*xd(1)
    c10 = values(mesh%idx%lxyz_inv(0 + nm(1), 1 + nm(2), 0 + nm(3)))*(CNST(1.0) - xd(1)) + &
      values(mesh%idx%lxyz_inv(1 + nm(1), 1 + nm(2), 0 + nm(3)))*xd(1)
    c01 = values(mesh%idx%lxyz_inv(0 + nm(1), 0 + nm(2), 1 + nm(3)))*(CNST(1.0) - xd(1)) + &
      values(mesh%idx%lxyz_inv(1 + nm(1), 0 + nm(2), 1 + nm(3)))*xd(1)
    c11 = values(mesh%idx%lxyz_inv(0 + nm(1), 1 + nm(2), 1 + nm(3)))*(CNST(1.0) - xd(1)) + &
      values(mesh%idx%lxyz_inv(1 + nm(1), 1 + nm(2), 1 + nm(3)))*xd(1)
    c0 = c00*(CNST(1.0) - xd(2)) + c10*xd(2)
    c1 = c01*(CNST(1.0) - xd(2)) + c11*xd(2)
    interpolated_value = c0*(CNST(1.0) - xd(3)) + c1*xd(3)

  case(2)

    ! bilinear interpolation: http://en.wikipedia.org/wiki/Bilinear_interpolation
    interpolated_value = &
      values(mesh%idx%lxyz_inv(0 + nm(1), 0 + nm(2), 0))*(CNST(1.0) - xd(1))*(CNST(1.0) - xd(2)) + &    
      values(mesh%idx%lxyz_inv(1 + nm(1), 0 + nm(2), 0))*xd(1)*(CNST(1.0) - xd(2)) + &    
      values(mesh%idx%lxyz_inv(0 + nm(1), 1 + nm(2), 0))*(CNST(1.0) - xd(1))*xd(2) + &    
      values(mesh%idx%lxyz_inv(1 + nm(1), 1 + nm(2), 0))*xd(1)*xd(2)

  case(1)
    
    ! linear interpolation
    interpolated_value = &
      values(mesh%idx%lxyz_inv(0 + nm(1), 0, 0))*(CNST(1.0) - xd(1)) + &
      values(mesh%idx%lxyz_inv(1 + nm(1), 0, 0))*xd(1)

  end select
  
  POP_SUB(X(interpolation_point_evaluate))

end function X(mesh_interpolation_evaluate)

! --------------------------------------------------------------

subroutine X(mesh_interpolation_test)(mesh)
  type(mesh_t), intent(in) :: mesh

  integer :: ip, iunit, idir, itest
  FLOAT :: xx(1:MAX_DIM)
  R_TYPE :: coeff(1:MAX_DIM), calculated, interpolated
  R_TYPE, allocatable :: ff(:)
  type(c_ptr)  :: random_gen_pointer
  type(mesh_interpolation_t) :: interp
  integer, parameter :: ntest_points = 20

  SAFE_ALLOCATE(ff(1:mesh%np_part))

  call loct_ran_init(random_gen_pointer)

  do idir = 1, mesh%sb%dim
    coeff(idir) = loct_ran_gaussian(random_gen_pointer, CNST(100.0)) + M_ZI*loct_ran_gaussian(random_gen_pointer, CNST(100.0))
  end do
  
  do ip = 1, mesh%np_part
    ff(ip) = sum(coeff(1:mesh%sb%dim)*mesh%x(ip, 1:mesh%sb%dim))
  end do

  call mesh_interpolation_init(interp, mesh)

  do itest = 1, ntest_points

    ip = 1 + nint(loct_ran_flat(random_gen_pointer, CNST(0.0), CNST(1.0))*(mesh%np - 1))

    do idir = 1, mesh%sb%dim
      xx(idir) = mesh%x(ip, idir) + mesh%spacing(idir)*loct_ran_flat(random_gen_pointer, CNST(-1.0), CNST(1.0))
    end do

    calculated = sum(coeff(1:mesh%sb%dim)*xx(1:mesh%sb%dim))
    interpolated = mesh_interpolation_evaluate(interp, ff, xx)
    call messages_write('Random point ')
    call messages_write(itest, fmt = '(i3)')
    call messages_write(' error:')
    call messages_write(abs(calculated - interpolated), fmt = '(e8.2)', align_left = .true.)
    call messages_info()

  end do
  
  call loct_ran_end(random_gen_pointer)

  call mesh_interpolation_end(interp)

  SAFE_DEALLOCATE_A(ff)

end subroutine X(mesh_interpolation_test)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:


