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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module xyz_adjust_m
  use datasets_m
  use global_m
  use geometry_m
  use lalg_adv_m
  use parser_m
  use messages_m
  use species_m
  use unit_m
  use unit_system_m
  use utils_m

  implicit none

  private
  public :: &
    xyz_adjust_it,  &
    geometry_rotate

contains

  ! ---------------------------------------------------------
  subroutine xyz_adjust_it(geo, rotate)
    type(geometry_t),           intent(inout) :: geo
    logical,          optional, intent(in)    :: rotate

    integer, parameter :: &
      NONE    = 0,        &
      INERTIA = 1,        &
      PSEUDO  = 2,        &
      LARGE   = 3

    FLOAT :: center(MAX_DIM), x1(MAX_DIM), x2(MAX_DIM), to(MAX_DIM)
    integer :: axis_type, idir, default
    type(block_t) :: blk

    PUSH_SUB(xyz_adjust_it)

    ! is there something to do
    if(geo%natoms <= 1) return

    if(optional_default(rotate, .true.)) then

      ! get to axis
      if(parse_block(datasets_check('MainAxis'), blk)==0) then
        do idir = 1, geo%space%dim
          call parse_block_float(blk, 0, idir - 1, to(idir))
        enddo
        call parse_block_end(blk)
      else
        to(:) = M_ZERO
        to(1) = M_ONE
      end if
      to = to / sqrt(sum(to(1:geo%space%dim)**2))

      write(message(1),'(a,6f15.6)') 'Using main axis ', to(1:geo%space%dim)
      call messages_info(1)

      !%Variable MainAxis
      !%Type block
      !%Section Utilities::oct-center-geom
      !%Description 
      !% A vector of reals defining the axis to which the molecule
      !% should be aligned. If not present, the default value will
      !% be the x-axis. For example in 3D:
      !% <tt>%MainAxis
      !% <br> 1 | 0 | 0 
      !% <br>%</tt>
      !%End

      !%Variable AxisType
      !%Type integer
      !%Default inertia
      !%Section Utilities::oct-center-geom
      !%Description
      !% After the structure is centered, it is also aligned to a set of orthogonal axes.
      !% This variable decides which set of axes to use. Only implemented for 3D; otherwise
      !% 'none' is default and the only legal value.
      !%Option none 0
      !% Do not rotate. Will still give output regarding center of mass and moment of inertia.
      !%Option inertia 1
      !% The axis of inertia.
      !%Option pseudo_inertia 2
      !% Pseudo-axis of inertia, calculated considering all species to have equal mass.
      !%Option large_axis 3
      !% The larger axis of the molecule.
      !%End
      if(geo%space%dim == 3) then
        default = INERTIA
      else
        default = NONE
      endif
      call parse_integer(datasets_check('AxisType'), default, axis_type)
      call messages_print_var_option(stdout, "AxisType", axis_type)

      if(geo%space%dim /= 3 .and. axis_type /= NONE) then
        call messages_not_implemented("alignment to axes (AxisType /= none) in other than 3 dimensions")
      endif

      select case(axis_type)
      case(NONE, INERTIA, PSEUDO)
        call find_center_of_mass(geo, center, pseudo = (axis_type==PSEUDO))

        write(message(1),'(3a,99f15.6)') 'Center of mass [', trim(units_abbrev(units_out%length)), '] = ', &
          (units_from_atomic(units_out%length, center(idir)), idir = 1, geo%space%dim)
        call messages_info(1)

        call translate(geo, center)
        call axis_inertia(geo, x1, x2, pseudo = (axis_type==PSEUDO))
      case(LARGE)
        call find_center(geo, center)

        write(message(1),'(3a,99f15.6)') 'Center [', trim(units_abbrev(units_out%length)), '] = ', &
          (units_from_atomic(units_out%length, center(idir)), idir = 1, geo%space%dim)
        call messages_info(1)

        call translate(geo, center)
        call axis_large(geo, x1, x2)
      case default
        write(message(1), '(a,i2,a)') 'AxisType = ', axis_type, ' not known by Octopus.'
        call messages_fatal(1)
      end select

      write(message(1),'(a,99f15.6)') 'Found primary   axis ', x1(1:geo%space%dim)
      write(message(2),'(a,99f15.6)') 'Found secondary axis ', x2(1:geo%space%dim)
      call messages_info(2)

      if(axis_type /= NONE) call geometry_rotate(geo, x1, x2, to)

    end if

    ! recenter
    call find_center(geo, center)
    call translate(geo, center)

    POP_SUB(xyz_adjust_it)
  end subroutine xyz_adjust_it


  ! ---------------------------------------------------------
  subroutine find_center(geo, x)
    type(geometry_t), intent(in) :: geo
    FLOAT, intent(out) :: x(MAX_DIM)

    FLOAT :: xmin(MAX_DIM), xmax(MAX_DIM)
    integer  :: i, j

    PUSH_SUB(find_center)

    xmin =  CNST(1e10)
    xmax = -CNST(1e10)
    do i = 1, geo%natoms
      do j = 1, geo%space%dim
        if(geo%atom(i)%x(j) > xmax(j)) xmax(j) = geo%atom(i)%x(j)
        if(geo%atom(i)%x(j) < xmin(j)) xmin(j) = geo%atom(i)%x(j)
      end do
    end do

    x = M_ZERO
    x(1:geo%space%dim) = (xmax(1:geo%space%dim) + xmin(1:geo%space%dim))/M_TWO

    POP_SUB(find_center)
  end subroutine find_center


  ! ---------------------------------------------------------
  subroutine find_center_of_mass(geo, x, pseudo)
    type(geometry_t), intent(in) :: geo
    FLOAT,  intent(out) :: x(MAX_DIM)
    logical, intent(in) :: pseudo

    FLOAT :: sm, m
    integer  :: i

    PUSH_SUB(find_center_of_mass)

    x = M_ZERO
    m = M_ONE
    sm = M_ZERO
    do i = 1, geo%natoms
      if(.not.pseudo) m = species_weight(geo%atom(i)%spec)
      sm = sm + m
      x = x + m * geo%atom(i)%x
    end do

    x = x / sm

    POP_SUB(find_center_of_mass)
  end subroutine find_center_of_mass


  ! ---------------------------------------------------------
  subroutine axis_large(geo, x, x2)
    type(geometry_t), intent(in) :: geo
    FLOAT, intent(out) :: x(MAX_DIM), x2(MAX_DIM)

    integer  :: i, j
    FLOAT :: rmax, r, r2

    PUSH_SUB(axis_large)

    ! first get the further apart atoms
    rmax = -CNST(1e10)
    do i = 1, geo%natoms
      do j = 1, geo%natoms/2 + 1
        r = sqrt(sum((geo%atom(i)%x(1:geo%space%dim)-geo%atom(j)%x(1:geo%space%dim))**2))
        if(r > rmax) then
          rmax = r
          x = geo%atom(i)%x - geo%atom(j)%x
        end if
      end do
    end do
    x  = x /sqrt(sum(x(1:geo%space%dim)**2))

    ! now let us find out what is the second most important axis
    rmax = -CNST(1e10)
    do i = 1, geo%natoms
      r2 = sum(x(1:geo%space%dim) * geo%atom(i)%x(1:geo%space%dim))
      r = sqrt(sum((geo%atom(i)%x(1:geo%space%dim) - r2*x(1:geo%space%dim))**2))
      if(r > rmax) then
        rmax = r
        x2 = geo%atom(i)%x - r2*x
      end if
    end do

    POP_SUB(axis_large)
  end subroutine axis_large


  ! ---------------------------------------------------------
  !> This subroutine assumes that the origin of the coordinates is the
  !! center of mass of the system
  subroutine axis_inertia(geo, x, x2, pseudo)
    type(geometry_t), intent(in) :: geo
    FLOAT,  intent(out) :: x(MAX_DIM), x2(MAX_DIM)
    logical, intent(in) :: pseudo

    FLOAT :: m, tinertia(MAX_DIM, MAX_DIM), eigenvalues(MAX_DIM)
    integer :: ii, jj, iatom
    type(unit_t) :: unit

    PUSH_SUB(axis_inertia)

    ! first calculate the inertia tensor
    tinertia = M_ZERO
    m = M_ONE
    do iatom = 1, geo%natoms
      if(.not.pseudo) m = species_weight(geo%atom(iatom)%spec)
      do ii = 1, geo%space%dim
        do jj = 1, geo%space%dim
          tinertia(ii, jj) = tinertia(ii, jj) - m*geo%atom(iatom)%x(ii)*geo%atom(iatom)%x(jj)
        end do
        tinertia(ii, ii) = tinertia(ii, ii) + m*sum(geo%atom(iatom)%x(:)**2)
      end do
    end do

    unit = units_out%length**2
    ! note: we always use amu for atomic masses, so no unit conversion to/from atomic is needed.
    if(pseudo) then
      write(message(1),'(a)') 'Moment of pseudo-inertia tensor [' // trim(units_abbrev(unit)) // ']'
    else
      write(message(1),'(a)') 'Moment of inertia tensor [amu*' // trim(units_abbrev(unit)) // ']'
    endif
    call messages_info(1)
    call output_tensor(stdout, tinertia, geo%space%dim, unit, write_average = .true.)

    call lalg_eigensolve(geo%space%dim, tinertia, eigenvalues)

    write(message(1),'(a,6f25.6)') 'Eigenvalues: ', &
      (units_from_atomic(unit, eigenvalues(jj)), jj = 1, geo%space%dim)
    call messages_info(1)

    ! make a choice to fix the sign of the axis.
    do ii = 1, 2
      jj = maxloc(abs(tinertia(:,ii)), dim = 1)
      if(tinertia(jj,ii) < M_ZERO) tinertia(:,ii) = -tinertia(:,ii)
    enddo
    x  = tinertia(:,1)
    x2 = tinertia(:,2)

    POP_SUB(axis_inertia)
  end subroutine axis_inertia


  ! ---------------------------------------------------------
  subroutine translate(geo, x)
    type(geometry_t), intent(inout) :: geo
    FLOAT,            intent(in)    :: x(MAX_DIM)

    integer  :: iatom

    PUSH_SUB(translate)

    do iatom = 1, geo%natoms
      geo%atom(iatom)%x = geo%atom(iatom)%x - x
    end do
    do iatom = 1, geo%ncatoms
      geo%catom(iatom)%x = geo%catom(iatom)%x - x
    end do

    POP_SUB(translate)
  end subroutine translate


  ! ---------------------------------------------------------
  subroutine geometry_rotate(geo, from, from2, to)
    type(geometry_t), intent(inout) :: geo
    FLOAT,            intent(in)    :: from(MAX_DIM)   !< assumed to be normalized
    FLOAT,            intent(in)    :: from2(MAX_DIM)  !< assumed to be normalized
    FLOAT,            intent(in)    :: to(MAX_DIM)     !< assumed to be normalized

    integer :: iatom, idim
    FLOAT :: m1(MAX_DIM, MAX_DIM), m2(MAX_DIM, MAX_DIM)
    FLOAT :: m3(MAX_DIM, MAX_DIM), f2(MAX_DIM), per(MAX_DIM)
    FLOAT :: alpha, r

    PUSH_SUB(geometry_rotate)

    if(geo%space%dim /= 3) &
      call messages_not_implemented("geometry_rotate in other than 3 dimensions")

    ! initialize matrices
    m1 = M_ZERO
    do idim = 1, MAX_DIM
      m1(idim, idim) = M_ONE
    end do

    ! rotate the to-axis to the z-axis
    if(to(2) /= M_ZERO) then
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

    ! rotate perpendicular axis to the y-axis
    m2 = M_ZERO; m2(1,1) = M_ONE; m2(2,2) = M_ONE; m2(3,3) = M_ONE
    alpha = atan2(per(1), per(2))
    call rotate_z(m2, -alpha)

    ! rotate from => to (around the y-axis)
    m3 = M_ZERO; m3(1,1) = M_ONE; m3(2,2) = M_ONE; m3(3,3) = M_ONE
    alpha = acos(sum(from*to))
    call rotate_y(m3, -alpha)

    ! join matrices
    m2 = matmul(transpose(m2), matmul(m3, m2))

    ! rotate around the z-axis to get the second axis
    per = matmul(m2, matmul(m1, from2))
    alpha = atan2(per(1), per(2))
    call rotate_z(m2, -alpha) ! second axis is now y

    ! get combined transformation
    m1 = matmul(transpose(m1), matmul(m2, m1))

    ! now transform the coordinates
    ! it is written in this way to avoid what I consider a bug in the Intel compiler
    do iatom = 1, geo%natoms
      f2 = geo%atom(iatom)%x
      geo%atom(iatom)%x = matmul(m1, f2)
    end do

    do iatom = 1, geo%ncatoms
      f2 = geo%catom(iatom)%x
      geo%catom(iatom)%x = matmul(m1, f2)
    end do

    POP_SUB(geometry_rotate)
  end subroutine geometry_rotate


  ! ---------------------------------------------------------
  subroutine rotate_x(m, angle)
    FLOAT, intent(inout) :: m(MAX_DIM, MAX_DIM)
    FLOAT, intent(in)    :: angle

    FLOAT :: aux(MAX_DIM, MAX_DIM), ca, sa

    PUSH_SUB(rotate_x)

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(1, 1) = M_ONE
    aux(2, 2) = ca
    aux(3, 3) = ca
    aux(2, 3) = sa
    aux(3, 2) = -sa

    m = matmul(aux, m)

    POP_SUB(rotate_x)
  end subroutine rotate_x


  ! ---------------------------------------------------------
  subroutine rotate_y(m, angle)
    FLOAT, intent(inout) :: m(MAX_DIM, MAX_DIM)
    FLOAT, intent(in)    :: angle

    FLOAT :: aux(MAX_DIM, MAX_DIM), ca, sa

    PUSH_SUB(rotate_y)

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(2, 2) = M_ONE
    aux(1, 1) = ca
    aux(3, 3) = ca
    aux(1, 3) = sa
    aux(3, 1) = -sa

    m = matmul(aux, m)

    POP_SUB(rotate_y)
  end subroutine rotate_y


  ! ---------------------------------------------------------
  subroutine rotate_z(m, angle)
    FLOAT, intent(inout) :: m(MAX_DIM, MAX_DIM)
    FLOAT,    intent(in) :: angle

    FLOAT :: aux(MAX_DIM, MAX_DIM), ca, sa

    PUSH_SUB(rotate_z)

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(3, 3) = M_ONE
    aux(1, 1) = ca
    aux(2, 2) = ca
    aux(1, 2) = sa
    aux(2, 1) = -sa

    m = matmul(aux, m)

    POP_SUB(rotate_z)
  end subroutine rotate_z

end module xyz_adjust_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
