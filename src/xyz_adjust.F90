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

module xyz_adjust
  use geometry

  implicit none

  private
  public :: xyz_adjust_it

contains

  subroutine xyz_adjust_it(geo)
    type(geometry_type), intent(inout) :: geo

    integer, parameter :: &
         INERTIA = 1, &
         PSEUDO  = 2, &
         LARGE   = 3

    FLOAT :: center(3), x1(3), x2(3), to(3)
    character(len=80) :: str
    integer :: axis_type

    ! is there something to do
    if(geo%natoms <= 1) return
  
    ! get to axis
    str = "MainAxis"
    if(loct_parse_isdef(str) .ne. 0) then
      call loct_parse_block_float(str, 0, 0, to(1))
      call loct_parse_block_float(str, 0, 1, to(2))
      call loct_parse_block_float(str, 0, 2, to(3))
    else
      to(1) = M_ONE; to(2) = M_ZERO; to(3) = M_ZERO
    end if
    to = to / sqrt(sum(to**2))
    
    call loct_parse_int("AxisType", INERTIA, axis_type)

    select case(axis_type)
    case(INERTIA, PSEUDO)
      call find_center_of_mass(geo, center, (axis_type==PSEUDO))
      call translate(geo, -center)
      call axis_inertia(geo, x1, x2, (axis_type==PSEUDO))
    case(LARGE)
      call find_center(geo, center)
      call translate(geo, -center)
      call axis_large(geo, x1, x2)
    case default
      write(message(1), '(a,i2,a)') 'AxisType = ', axis_type, ' not known by octopus'
      call write_fatal(1)
    end select

    ! check if the axis are OK
    if(sum(x2**2) == M_ZERO) then ! linear molecule
      if(x1(1) == M_ZERO) then
        x2(1) = x1(1); x2(2) = -x1(3); x2(3) = x1(2)
      else if(x1(2) == M_ZERO) then
        x2(2) = x1(2); x2(1) = -x1(3); x2(3) = x1(1)
      else
        x2(3) = x1(3); x2(1) = -x1(2); x2(2) = x1(1)
      end if
    end if
    x2 = x2/sqrt(sum(x2**2))

    ! rotate to main axis
    call rotate(geo, x1, x2, to)
    
    ! recenter
    call find_center(geo, center)
    call translate(geo, -center)
  end subroutine xyz_adjust_it
    
  subroutine find_center(geo, x)
    type(geometry_type), intent(IN) :: geo
    FLOAT, intent(out) :: x(3)

    FLOAT :: xmin(3), xmax(3)
    integer  :: i, j

    xmin =  CNST(1e10)
    xmax = -CNST(1e10)
    do i = 1, geo%natoms
      do j = 1, 3
        if(geo%atom(i)%x(j) > xmax(j)) xmax(j) = geo%atom(i)%x(j)
        if(geo%atom(i)%x(j) < xmin(j)) xmin(j) = geo%atom(i)%x(j)
      end do
    end do
    
    x = (xmax + xmin)/M_TWO
  end subroutine find_center

  subroutine find_center_of_mass(geo, x, pseudo)
    type(geometry_type), intent(IN) :: geo
    FLOAT, intent(out) :: x(3)
    logical, intent(in) :: pseudo

    FLOAT :: sm, m
    integer  :: i

    x = M_ZERO
    m = M_ONE
    do i = 1, geo%natoms
      if(.not.pseudo) m = geo%atom(i)%spec%weight
      sm = sm + m
      x = x + m * geo%atom(i)%x
    end do
    
    x = x / sm
  end subroutine find_center_of_mass

  subroutine axis_large(geo, x, x2)
    type(geometry_type), intent(IN) :: geo
    FLOAT, intent(out) :: x(3), x2(3)

    integer  :: i, j
    FLOAT :: rmax, r, r2
    
    ! first get the further apart atoms
    rmax = -CNST(1e10)
    do i = 1, geo%natoms
      do j = 1, geo%natoms/2 + 1
        r = sqrt(sum((geo%atom(i)%x-geo%atom(j)%x)**2))
        if(r > rmax) then
          rmax = r
          x = geo%atom(i)%x - geo%atom(j)%x
        end if
      end do
    end do
    x  = x /sqrt(sum(x**2))
    
    ! now let us find out what is the second most important axis
    rmax = -CNST(1e10)
    do i = 1, geo%natoms
      r2 = sum(x * geo%atom(i)%x)
      r = sqrt(sum((geo%atom(i)%x - r2*x)**2))
      if(r > rmax) then
        rmax = r
        x2 = geo%atom(i)%x - r2*x
      end if
    end do    
  end subroutine axis_large

  ! This subroutine assumes that the origin of the coordinates is the
  ! center of mass of the system
  subroutine axis_inertia(geo, x, x2, pseudo)
    type(geometry_type), intent(IN) :: geo
    FLOAT, intent(out) :: x(3), x2(3)
    logical, intent(in) :: pseudo

    FLOAT :: m, tinertia(3,3), vec(3,3), e(3)
    integer :: i, j, iatom

    ! first calculate the inertia tensor
    tinertia = M_ZERO
    m = M_ONE
    do iatom = 1, geo%natoms
      if(.not.pseudo) m = geo%atom(iatom)%spec%weight
      do i = 1, 3
        do j = 1, 3
          tinertia(i, j) = tinertia(i, j) - m*geo%atom(iatom)%x(i)*geo%atom(iatom)%x(j)
        end do
        tinertia(i, i) = tinertia(i, i) + M_THREE*m*sqrt(sum(geo%atom(iatom)%x(:)**2))
      end do
    end do

    call lalg_eigensolve(3, tinertia, vec, e)

    x  = vec(:,1) / sqrt(sum(vec(:,1)**2))
    x2 = vec(:,2) / sqrt(sum(vec(:,2)**2))
  end subroutine axis_inertia
  
  subroutine translate(geo, x)
    type(geometry_type), intent(inout) :: geo
    FLOAT, intent(in) :: x(3)
    
    integer  :: i
    
    do i = 1, geo%natoms
      geo%atom(i)%x = geo%atom(i)%x + x
    end do
    do i = 1, geo%ncatoms
      geo%catom(i)%x = geo%catom(i)%x + x
    end do
  end subroutine translate

  subroutine rotate(geo, from, from2, to)
    type(geometry_type), intent(inout) :: geo
    FLOAT, intent(in) :: from(3), from2(3), to(3) ! assumed to be normalize
    
    integer :: i
    FLOAT :: m1(3,3), m2(3,3), m3(3,3), f2(3), per(3)
    FLOAT :: alpha, r
    
    ! initialize matrices
    m1 = M_ZERO; m1(1,1) = M_ONE; m1(2,2) = M_ONE; m1(3,3) = M_ONE
    
    ! rotate the to axis to the z axis
    if(to(2).ne.M_ZERO) then
      alpha = atan2(to(2), to(1))
      call rotate_z(m1, alpha)
    end if
    alpha = atan2(sqrt(to(1)**2 + to(2)**2), to(3))
    call rotate_y(m1, -alpha)

    ! get perpendicular to z and from
    f2 = matmul(m1, from)
    per(1) = -f2(2)
    per(2) =  f2(1)
    per(3) = M_ZERO
    r = sqrt(sum(per**2))
    if(r > M_ZERO) then
      per = per/r
    else
      per(2) = M_ONE
    end if

    ! rotate perpendicular axis to the y axis
    m2 = M_ZERO; m2(1,1) = M_ONE; m2(2,2) = M_ONE; m2(3,3) = M_ONE
    alpha = atan2(per(1), per(2))
    call rotate_z(m2, -alpha)

    ! rotate from => to (around the y axis)
    m3 = M_ZERO; m3(1,1) = M_ONE; m3(2,2) = M_ONE; m3(3,3) = M_ONE
    alpha = acos(sum(from*to))
    call rotate_y(m3, -alpha)
    
    ! join matrices
    m2 = matmul(transpose(m2), matmul(m3, m2))

    ! rotate around the z axis to get the second axis
    per = matmul(m2, matmul(m1, from2))
    alpha = atan2(per(1), per(2))
    call rotate_z(m2, -alpha) ! second axis is now y

    ! get combined transformation
    m1 = matmul(transpose(m1), matmul(m2, m1))

    ! now transform the coordinates
    ! it is written in this way to avoid what I consider a bug in the Intel compiler
    do i = 1, geo%natoms
      f2 = geo%atom(i)%x
      geo%atom(i)%x = matmul(m1, f2)
    end do

    do i = 1, geo%ncatoms
      f2 = geo%catom(i)%x
      geo%catom(i)%x = matmul(m1, f2)
    end do

  end subroutine rotate

  subroutine rotate_x(m, angle)
    FLOAT, intent(inout) :: m(3,3)
    FLOAT, intent(in)    :: angle
    
    FLOAT :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(1, 1) = M_ONE
    aux(2, 2) = ca
    aux(3, 3) = ca
    aux(2, 3) = sa
    aux(3, 2) = -sa

    m = matmul(aux, m)
  end subroutine rotate_x

  subroutine rotate_y(m, angle)
    FLOAT, intent(inout) :: m(3,3)
    FLOAT, intent(in)    :: angle

    FLOAT :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(2, 2) = M_ONE
    aux(1, 1) = ca
    aux(3, 3) = ca
    aux(1, 3) = sa
    aux(3, 1) = -sa

    m = matmul(aux, m)
  end subroutine rotate_y

  subroutine rotate_z(m, angle)
    FLOAT, intent(inout) :: m(3,3)
    FLOAT, intent(in) :: angle

    FLOAT :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(3, 3) = M_ONE
    aux(1, 1) = ca
    aux(2, 2) = ca
    aux(1, 2) = sa
    aux(2, 1) = -sa

    m = matmul(aux, m)
  end subroutine rotate_z
  
end module xyz_adjust
