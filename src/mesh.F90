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

#include "global.h"

module mesh
use global
use oct_parser
use units
use math
use atom
  
implicit none

integer, parameter :: &
     SPHERE   = 1,    &
     CYLINDER = 2,    &
     MINIMUM  = 3,    &
     PARALLELEPIPED = 4

type lookup_type
  integer :: n
  integer, pointer :: i(:)
  FLOAT, pointer :: w(:)
end type lookup_type

type derivatives_type
  integer :: norder
  type(lookup_type), pointer :: lookup(:)
end type derivatives_type

type mesh_type
  integer  :: box_shape ! 1->sphere, 2->cylinder, 3->sphere around each atom,
                        ! 4->parallelpiped (orthonormal, up to now).
  FLOAT :: h(3)      ! the (constant) spacing between the points
  
  FLOAT :: rsize     ! the radius of the sphere or of the cylinder
  FLOAT :: xsize     ! the length of the cylinder in the x direction
  FLOAT :: lsize(3)  ! half of the length of the parallelepiped in each direction.
  FLOAT :: rlat(3,3)   ! lattice primitive vectors
  FLOAT :: klat(3,3)   ! reciprocal lattice primitive vectors
  FLOAT :: shift(27,3) ! shift to equivalent positions in nearest neighbour primitive cells
  
  integer  :: np        ! number of points in mesh
  integer  :: np_in     ! number of points without the enlargement

  integer, pointer :: Lxyz(:,:)       ! return x, y and z for each point
  integer, pointer :: Lxyz_inv(:,:,:) ! return points # for each xyz

  ! some other vars
  integer :: nr(2,3)  ! dimensions of the box where the points are contained
  integer :: l(3)     ! literally n(2,:) - n(1,:) + 1

  type(derivatives_type) :: laplacian
  type(derivatives_type), pointer :: grad(:)

  FLOAT :: fft_alpha ! enlargement factor for double box

  FLOAT :: vol_pp    ! element of volume for integrations
end type mesh_type

contains

subroutine mesh_init(m, natoms, atom)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)

  integer :: i, k, morder
  logical :: fft_optimize

  call push_sub('mesh_init')

  call mesh_create(m, natoms, atom)

  call oct_parse_float('DoubleFFTParameter', M_TWO, m%fft_alpha)
  if (m%fft_alpha < M_ONE .or. m%fft_alpha > M_THREE ) then
    write(message(1), '(a,f12.5,a)') "Input: '", m%fft_alpha, &
         "' is not a valid DoubleFFTParameter"
    message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
    call write_fatal(2)
  end if

  nullify(m%grad)

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

  write(message(1),'(3a,f9.3,a)') '  Grid Cutoff [',trim(units_out%energy%abbrev),'] = ', &
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
  FLOAT, intent(out) :: x(conf%dim)

  x(1:conf%dim) = m%Lxyz(1:conf%dim, i)*m%h(1:conf%dim)
end subroutine mesh_xyz

subroutine mesh_x(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  FLOAT, intent(out) :: x

  x = m%Lxyz(1, i)*m%h(1)
end subroutine mesh_x

subroutine mesh_y(m, i, y)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  FLOAT, intent(out) :: y

  y = m%Lxyz(2, i)*m%h(2)
end subroutine mesh_y

subroutine mesh_z(m, i, z)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  FLOAT, intent(out) :: z

  z = m%Lxyz(3, i)*m%h(3)
end subroutine mesh_z

subroutine mesh_r(m, i, r, a, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  FLOAT, intent(out) :: r
  FLOAT, intent(in), optional :: a(conf%dim)
  FLOAT, intent(out), optional :: x(conf%dim)

  FLOAT :: xx(conf%dim)

  call mesh_xyz(m, i, xx)
  if(present(a)) xx = xx - a
  r = sqrt(dot_product(xx))

  if(present(x)) x = xx
end subroutine mesh_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! /* Finds out if a given point of a mesh belongs to the "border" of the mesh.
! A point belongs to the border of the mesh if it is too close to any of the
! walls of the mesh. The criterium is set by input parameter "width".
!
! m     : the mesh.
! i     : the point in the mesh.
! n     : on output, the number (0<=n<=3) of "walls" of the mesh that the point is
!         too close to, in order to consider it belonging to a mesh.
! d     : the distances of the point to the walls, for each of the walls that the
!         point is too close to.
! width : the width of the border.
!
! So, if n>0, the point is in the border.  */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mesh_inborder(m, i, n, d, width)
  type(mesh_type), intent(in) :: m
  integer, intent(in)   :: i
  FLOAT, intent(in)  :: width
  integer, intent(out)  :: n
  FLOAT, intent(out) :: d(3)

  integer :: j
  FLOAT :: x(3), r, dd

  call mesh_r(m, i, r, x=x)
  n = 0
  select case(m%box_shape)
    case(SPHERE)
      dd = r - (m%rsize - width)
      if(dd.gt.M_ZERO) then
        n = 1; d(1) = dd
      end if
    case(CYLINDER)
      dd = sqrt(x(2)**2 + x(3)**2) - (m%rsize - width)
      if(dd.gt.M_ZERO) then
        n = 1; d(1) = dd
      endif
      dd = abs(x(1)) - (m%xsize - width)
      if(dd.gt.M_ZERO) then
        n = n + 1; d(n) = dd
      endif
    case(PARALLELEPIPED)
      do j = 1, conf%dim
         dd = abs(x(j)) - (m%lsize(j) - width)
         if(dd.gt.M_ZERO) then
           n = n + 1; d(n) = dd
         end if
      end do
  end select

end subroutine mesh_inborder

subroutine mesh_end(m)
  type(mesh_type), intent(inout) :: m

  integer :: i, j

  call push_sub('mesh_end')

  if(associated(m%Lxyz)) then
    deallocate(m%Lxyz, m%Lxyz_inv)
    nullify(m%Lxyz, m%Lxyz_inv)
  end if

  call pop_sub()
  return
end subroutine mesh_end

#include "mesh_create.F90"

end module mesh
