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
use oct_parser
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

type der_lookup_type
  integer :: lapl_n, grad_n(3)
  integer,  pointer :: lapl_i(:), grad_i(:,:)
  real(r8), pointer :: lapl_w(:), grad_w(:,:)
end type der_lookup_type

type mesh_derivatives_type
  integer :: space ! 0 for FD, 1 for reciprocal-space (TAKE AWAY)

  integer :: norder ! order on the discretization of the gradient/laplace operator
  real(r8), pointer :: dgidfj(:) ! the coefficients for the gradient
  real(r8), pointer :: dlidfj(:) ! the coefficients for the laplacian
  
end type mesh_derivatives_type

type mesh_type
  integer  :: box_shape ! 1->sphere, 2->cylinder, 3->sphere around each atom,
                        ! 4->parallelpiped (orthonormal, up to now).
  real(r8) :: h(3)      ! the (constant) spacing between the points
  
  real(r8) :: rsize     ! the radius of the sphere or of the cylinder
  real(r8) :: xsize     ! the length of the cylinder in the x direction
  real(r8) :: lsize(3)  ! half of the length of the parallelepiped in each direction.
  real(r8) :: rlat(3,3)   ! lattice primitive vectors
  real(r8) :: klat(3,3)   ! reciprocal lattice primitive vectors
  real(r8) :: shift(27,3) ! shift to equivalent positions in nearest neighbour primitive cells
  
  integer  :: np        ! number of points in mesh
  integer  :: np_in     ! number of points without the enlargement

  integer, pointer :: Lxyz(:,:)       ! return x, y and z for each point
  integer, pointer :: Lxyz_inv(:,:,:) ! return points # for each xyz

  ! some other vars
  integer :: nr(2,3)  ! dimensions of the box where the points are contained
  integer :: l(3)     ! literally n(2,:) - n(1,:) + 1

  type(mesh_derivatives_type), pointer :: d

  type(der_lookup_type), pointer :: der_lookup(:)   ! lookup tables for derivatives
  
  type(fft_type) :: dfft, zfft ! to calculate derivatives using ffts
  real(r8) :: fft_alpha ! enlargement factor for double box

  real(r8) :: vol_pp    ! element of volume for integrations
end type mesh_type

contains

subroutine mesh_init(m, natoms, atom)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)

  integer :: i, k, morder
  logical :: fft_optimize

  call push_sub('mesh_init')

  allocate(m%d)
  call oct_parse_int('OrderDerivatives', 4, m%d%norder)
  call mesh_init_derivatives_coeff(m%d)

  call mesh_create(m, natoms, atom)

  call oct_parse_double('DoubleFFTParameter', 2.0_r8, m%fft_alpha)
  if (m%fft_alpha < 1.0_r8 .or. m%fft_alpha > 3.0_r8 ) then
    write(message(1), '(a,f12.5,a)') "Input: '", m%fft_alpha, &
         "' is not a valid DoubleFFTParameter"
    message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
    call write_fatal(2)
  end if

  call oct_parse_int('DerivativesSpace', REAL_SPACE, m%d%space)
  if(m%d%space < 0 .or. m%d%space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", m%d%space, &
         "' is not a valid DerivativesSpace"
    message(2) = '(0 <= DerivativesSpace <=1)'
    call write_fatal(2)
  end if

  ! order is very important
  ! init rs -> create -> init fft
  if(m%d%space == REAL_SPACE) then 
    message(1) = 'Info: Derivatives calculated in real-space'
  else
    call fft_init(m%l, fft_real,    m%dfft)
    call fft_init(m%l, fft_complex, m%zfft)
    message(1) = 'Info: Derivatives calculated in reciprocal-space'
  end if
  call write_info(1)

  call pop_sub()
end subroutine mesh_init

! finds the dimension of a box doubled in the non-periodic dimensions
subroutine mesh_double_box(m, db)
  type(mesh_type), intent(in) :: m
  integer, intent(out) :: db(3)

  integer :: i

  db = 1
  do i = conf%periodic_dim+1, conf%dim
    db(i) = nint(m%fft_alpha*(m%nr(2,i)-m%nr(1,i))) + 1
  end do
  do i = 1, conf%periodic_dim ! for periodic systems
    db(i) = m%nr(2,i) - m%nr(1,i) + 1
  end do

end subroutine mesh_double_box

subroutine mesh_write_info(m, unit)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: unit

  character(len=15), parameter :: bs(4) = (/ &
      'sphere        ', &
      'cylinder      ', &
      'around nuclei ', &
      'parallelepiped'/)

  call push_sub('mesh_write_info')
  
  write(message(1), '(a,a,1x)') '  Type = ', bs(m%box_shape)
  if(m%box_shape == SPHERE .or. m%box_shape == CYLINDER  .or. m%box_shape == MINIMUM) then
    write(message(2), '(3a,f7.3)') '  Radius  [', trim(units_out%length%abbrev), '] = ', &
                                   m%rsize/units_out%length%factor
  endif
  if(m%box_shape == CYLINDER) then
    write(message(2), '(a,3a,f7.3)') trim(message(2)), ', xlength [', &
         trim(units_out%length%abbrev), '] = ', m%xsize/units_out%length%factor
  end if
  if(m%box_shape == PARALLELEPIPED) then
    write(message(2),'(3a, a, f6.3, a, f6.3, a, f6.3, a)') &
       '  Lengths [', trim(units_out%length%abbrev), '] = ',         &
       '(', m%lsize(1)/units_out%length%factor, ',',                     &
            m%lsize(2)/units_out%length%factor, ',',                     &
            m%lsize(3)/units_out%length%factor, ')'
  endif
  write(message(3),'(3a, a, f6.3, a, f6.3, a, f6.3, a, 1x, 3a, f8.5)') &
       '  Spacing [', trim(units_out%length%abbrev), '] = ',         &
       '(', m%h(1)/units_out%length%factor, ',',                     &
            m%h(2)/units_out%length%factor, ',',                     &
            m%h(3)/units_out%length%factor, ')',                     &
       '   volume/point [', trim(units_out%length%abbrev), '^3] = ', &
       m%vol_pp/units_out%length%factor**3
  write(message(4),'(a, i6)') '  # inner mesh = ', m%np
  call write_info(4, unit)

  write(message(1),'(a,f9.3,a)') '  Grid Cutoff ['+trim(units_out%energy%abbrev)+'] = ', &
                                 (M_PI**2/(M_TWO*maxval(m%h)**2))/units_out%energy%factor
  call write_info(1, unit)
      
  if (conf%periodic_dim > 0) then
    write(message(1),'(1x)')
    write(message(2),'(a,3a,a)')'Lattice Primitive Vectors [', trim(units_out%length%abbrev), ']'
    write(message(3),'(a,f8.3)')'    x axis ', &
                    m%rlat(1,1)/units_out%length%factor
    write(message(4),'(a,f8.3)')'    y axis ', &
                    m%rlat(2,2)/units_out%length%factor
    write(message(5),'(a,f8.3)')'    z axis ', &
                    m%rlat(3,3) /units_out%length%factor
    write(message(6),'(a,3a,a)') 'Reciprocal Lattice Primitive Vectors [', &
                    trim(units_out%length%abbrev), '^-1]'
    write(message(7),'(a,f8.3)')'  k_x axis ', &
                    m%klat(1,1)*units_out%length%factor
    write(message(8),'(a,f8.3)')'  k_y axis ', &
                    m%klat(2,2)*units_out%length%factor
    write(message(9),'(a,f8.3)')'  k_z axis ', &
                    m%klat(3,3)*units_out%length%factor
  call write_info(9, unit)
  end if


  call pop_sub()
  return
end subroutine mesh_write_info

! subroutines to get xyzr
subroutine mesh_xyz(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x(conf%dim)

  x(1:conf%dim) = m%Lxyz(1:conf%dim, i)*m%h(1:conf%dim)
end subroutine mesh_xyz

subroutine mesh_x(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x

  x = m%Lxyz(1, i)*m%h(1)
end subroutine mesh_x

subroutine mesh_y(m, i, y)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: y

  y = m%Lxyz(2, i)*m%h(2)
end subroutine mesh_y

subroutine mesh_z(m, i, z)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: z

  z = m%Lxyz(3, i)*m%h(3)
end subroutine mesh_z

subroutine mesh_r(m, i, r, a, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: r
  real(r8), intent(in), optional :: a(conf%dim)
  real(r8), intent(out), optional :: x(conf%dim)

  real(r8) :: xx(conf%dim)

  call mesh_xyz(m, i, xx)
  if(present(a)) xx = xx - a
  r = sqrt(sum(xx(:)**2))

  if(present(x)) x = xx
end subroutine mesh_r

subroutine mesh_init_derivatives_coeff(d)
  use math, only: weights
  type(mesh_derivatives_type), intent(inout) :: d

  integer :: k, i, j, morder
  real(r8), allocatable :: cc(:,:,:)

  k = d%norder
  if (k < 1) then
    write(message(1), '(a,i4,a)') "Input: '", k, "' is not a valid OrderDerivatives"
    message(2) = '(1 <= OrderDerivatives)'
    call write_fatal(2)
  end if
  allocate(d%dlidfj(-k:k), d%dgidfj(-k:k))
  morder = 2*k
  allocate(cc(0:morder, 0:morder, 0:2))
  call weights(2, morder, cc)
  d%dgidfj(0) = cc(0, morder, 1)
  d%dlidfj(0) = cc(0, morder, 2)
  
  j = 1
  do i = 1, k
    d%dgidfj(-i) = cc(j, morder, 1)
    d%dlidfj(-i) = cc(j, morder, 2)
    j = j + 1
    d%dgidfj( i) = cc(j, morder, 1)
    d%dlidfj( i) = cc(j, morder, 2)
    j = j + 1
  end do
  deallocate(cc)

end subroutine mesh_init_derivatives_coeff

! this actually adds to outp
subroutine phase_factor(m, n, vec, inp, outp)
  implicit none
  type(mesh_type), intent(IN) :: m
  integer, intent(in)         :: n(3)
  real(r8), intent(IN)        :: vec(3)
  complex(r8), intent(IN)     :: inp (n(1)/2+1, n(2), n(3))
  complex(r8), intent(inout)  :: outp(n(1)/2+1, n(2), n(3))
  
  complex(r8) :: k(3)
  integer     :: ix, iy, iz, ixx, iyy, izz
  
  k = M_z0
  k(1:conf%dim) = M_zI * ((2.0_r8*M_Pi)/(n(1:conf%dim)*m%h(1:conf%dim)))
  do iz = 1, n(3)
    izz = pad_feq(iz, n(3), .true.)
    do iy = 1, n(2)
      iyy = pad_feq(iy, n(2), .true.)
      do ix = 1, n(1)/2 + 1
        ixx = pad_feq(ix, n(1), .true.)
        outp(ix, iy, iz) = outp(ix, iy, iz) + &
             exp( -(k(1)*vec(1)*ixx + k(2)*vec(2)*iyy + k(3)*vec(3)*izz) ) * inp(ix, iy, iz)
      end do
    end do
  end do
end subroutine phase_factor

! Deallocates what has to be deallocated ;)
subroutine mesh_end(m)
  type(mesh_type), intent(inout) :: m

  integer :: i

  call push_sub('mesh_end')

  if(associated(m%Lxyz)) then
    deallocate(m%Lxyz, m%Lxyz_inv)
    nullify(m%Lxyz, m%Lxyz_inv)
  end if
  if(associated(m%d)) then
    if(m%d%space == 0) then
      deallocate(m%d%dgidfj, m%d%dlidfj)
      nullify(m%d%dgidfj, m%d%dlidfj)
    end if
    deallocate(m%d); nullify(m%d)
  end if

  if(associated(m%der_lookup)) then
    do i = 1, m%np
      deallocate(m%der_lookup(i)%lapl_i, m%der_lookup(i)%grad_i)
      deallocate(m%der_lookup(i)%lapl_w, m%der_lookup(i)%grad_w)
    end do
    deallocate(m%der_lookup)
    nullify(m%der_lookup)
  end if

  call pop_sub()
  return
end subroutine mesh_end

#include "mesh_create.F90"

#include "undef.F90"
#include "real.F90"
#include "mesh_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_inc.F90"

end module mesh
