!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "config_F90.h"
module mesh
use global
use units
use fft
use atom
  
implicit none

integer, parameter ::     &
     REAL_SPACE = 0,      &
     RECIPROCAL_SPACE = 1

integer, parameter :: &
     SPHERE   = 1,    &
     CYLINDER = 2,    &
     MINIMUM  = 3,    &
     PARALLELEPIPED = 4

type mesh_derivatives_type
  integer :: space ! 0 for FD, 1 for reciprocal-space

  integer :: norder ! order on the discretization of the gradient/laplace operator
  real(r8), pointer :: dgidfj(:) ! the coefficients for the gradient
  real(r8), pointer :: dlidfj(:) ! the coefficients for the laplacian
  
end type mesh_derivatives_type

type mesh_type
  integer  :: box_shape ! 1->sphere, 2->cylinder, 3->sphere around each atom,
                        ! 4->parallelpiped (orthonormal, up to now).
  real(r8) :: h(3)      ! the (constant) spacing between the points
  logical  :: iso
  
  real(r8) :: rsize     ! the radius of the sphere or of the cylinder
  real(r8) :: zsize     ! the length of the cylinder in the z direction
  real(r8) :: lsize(3)  ! the lengths of the parallelepipeds in each direction.
  
  integer  :: np        ! number of points in inner mesh
  real(r8) :: vol_pp    ! element of volume for integrations
  
  ! return x, y and z for each point
  integer, pointer :: Lx(:), Ly(:), Lz(:)
  ! return points # for each xyz
  integer, pointer :: Lxyz_inv(:,:,:)
  integer, pointer :: ind(:, :, :)

  ! mesh elargement... this is currently only used in solving
  ! the poisson equation with conjugated gradients...
  integer :: nk        ! number of points in outer mesh
  integer, pointer :: Kx(:), Ky(:), Kz(:)

  ! some other vars
  integer :: nr(3)  ! number of points per radius inside
  integer :: nx(3)  ! the # of points in the whole radius

  type(mesh_derivatives_type) :: d

  integer :: fft_n(3), hfft_n ! we always want to perform FFTs ;)
  integer(POINTER_SIZE) :: dplanf, dplanb, zplanf, zplanb
  real(r8) :: fft_alpha
  integer :: fft_n2(3), hfft_n2
  integer(POINTER_SIZE) :: dplanf2, dplanb2
end type mesh_type

contains

#if defined(ONE_D)
#  include "mesh1D.F90"
#elif defined(THREE_D)
#  include "mesh3D.F90"
#endif
#include "undef.F90"
#include "real.F90"
#include "mesh_inc.F90"
#include "undef.F90"
#include "complex.F90"
#include "mesh_inc.F90"
#include "undef.F90"

subroutine mesh_write_info(m, unit)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: unit

  character(len=15), parameter :: bs(3) = (/ &
      'sphere       ', &
      'cylinder     ', &
      'around nuclei'/)

  sub_name = 'mesh_write_info'; call push_sub()

  if(m%box_shape == SPHERE .or. m%box_shape == CYLINDER  .or. m%box_shape == MINIMUM) then
    write(message(1), '(a,a,1x,3a,f7.3)') '  Type = ', bs(m%box_shape), &
       ' Radius [', trim(units_out%length%abbrev), '] = ', m%rsize/units_out%length%factor
  endif
  if(m%box_shape == CYLINDER) then
    write(message(1), '(a,3a,f7.3)') trim(message(1)), ', zlength [', &
         trim(units_out%length%abbrev), '] = ', m%zsize/units_out%length%factor
  end if
  if(m%box_shape == PARALLELEPIPED) then
    write(message(1),'(3a, a, f6.3, a, f6.3, a, f6.3, a)') &
       '  Lengths [', trim(units_out%length%abbrev), '] = ',         &
       '(', m%lsize(1)/units_out%length%factor, ',',                     &
            m%lsize(2)/units_out%length%factor, ',',                     &
            m%lsize(3)/units_out%length%factor, ')'
  endif
  write(message(2),'(3a, a, f6.3, a, f6.3, a, f6.3, a, 1x, 3a, f8.5)') &
       '  Spacing [', trim(units_out%length%abbrev), '] = ',         &
       '(', m%h(1)/units_out%length%factor, ',',                     &
            m%h(2)/units_out%length%factor, ',',                     &
            m%h(3)/units_out%length%factor, ')',                     &
       '   volume/point [', trim(units_out%length%abbrev), '^3] = ', &
       (m%h(1)*m%h(2)*m%h(3))/units_out%length%factor**3
  write(message(3),'(a, i6, a, i6)') '  # inner mesh = ', m%np, &
      '   # outer mesh = ', m%nk

  call write_info(3, unit)

  call pop_sub()
  return
end subroutine mesh_write_info

! subroutines to get xyzr
subroutine mesh_xyz(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x(3)

  x(1) = m%Lx(i)*m%H(1)
  x(2) = m%Ly(i)*m%H(2)
  x(3) = m%Lz(i)*m%H(3)
end subroutine mesh_xyz

subroutine mesh_x(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x

  x = m%Lx(i)*m%H(1)
end subroutine mesh_x

subroutine mesh_y(m, i, y)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: y

  y = m%Ly(i)*m%H(2)
end subroutine mesh_y

subroutine mesh_z(m, i, z)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: z

  z = m%Lz(i)*m%H(3)
end subroutine mesh_z

subroutine mesh_r(m, i, r, a, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: r
  real(r8), intent(in), optional :: a(3)
  real(r8), intent(out), optional :: x(3)

  real(r8) :: xx(3)

  call mesh_xyz(m, i, xx)
  if(present(a)) xx = xx - a
  r = sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)

  if(present(x)) x = xx
end subroutine mesh_r

subroutine mesh_init_derivatives_coeff(m)
  use math, only: weights
  type(mesh_type), intent(inout) :: m

  integer :: k, i, j, morder
  real(r8), allocatable :: cc(:,:,:)

  call oct_parse_int(C_string('OrderDerivatives'), 4, k)
  m%d%norder = k
  if (k < 1) then
    write(message(1), '(a,i4,a)') "Input: '", k, "' is not a valid OrderDerivatives"
    message(2) = '(1 <= OrderDerivatives)'
    call write_fatal(2)
  end if
  allocate(m%d%dlidfj(-k:k), m%d%dgidfj(-k:k))
  morder = 2*k
  allocate(cc(0:morder, 0:morder, 0:2))
  call weights(2, morder, cc)
  m%d%dgidfj(0) = cc(0, morder, 1)
  m%d%dlidfj(0) = cc(0, morder, 2)
  
  j = 1
  do i = 1, k
    m%d%dgidfj(-i) = cc(j, morder, 1)
    m%d%dlidfj(-i) = cc(j, morder, 2)
    j = j + 1
    m%d%dgidfj( i) = cc(j, morder, 1)
    m%d%dlidfj( i) = cc(j, morder, 2)
    j = j + 1
  end do
  deallocate(cc)

end subroutine mesh_init_derivatives_coeff

! Deallocates what has to be deallocated ;)
subroutine mesh_end(m)
  type(mesh_type), intent(inout) :: m

  sub_name = 'mesh_end'; call push_sub()

  if(associated(m%Lx)) then
    deallocate(m%Lx, m%Ly, m%Lz, m%Kx, m%Ky, m%Kz)
    nullify(m%Lx, m%Ly, m%Lz, m%Kx, m%Ky, m%Kz)
  end if
  if(associated(m%Lxyz_inv)) then
    deallocate(m%Lxyz_inv, m%ind)
    nullify(m%Lxyz_inv, m%ind)
  end if
  if(m%d%space == 0 .and. associated(m%d%dgidfj)) then
    deallocate(m%d%dgidfj, m%d%dlidfj)
    nullify(m%d%dgidfj, m%d%dlidfj)
  end if

  if(m%dplanf.ne.int(-1, POINTER_SIZE)) then
    call fftw_f77_destroy_plan(m%dplanf)
    call fftw_f77_destroy_plan(m%dplanb)
    call fftw_f77_destroy_plan(m%zplanf)
    call fftw_f77_destroy_plan(m%zplanb)
  end if

  if(m%dplanf2.ne.int(-1, POINTER_SIZE)) then
    call fftw_f77_destroy_plan(m%dplanf2)
    call fftw_f77_destroy_plan(m%dplanb2)
  end if

  call pop_sub()
  return
end subroutine mesh_end

end module mesh
