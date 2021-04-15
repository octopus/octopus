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

#include "global.h"

program centergeom
  use command_line_oct_m
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  integer :: ierr

  call global_init(is_serial = .true.)

  call getopt_init(ierr)
  if(ierr == 0) call getopt_center_geom()
  call getopt_end()

  call parser_init()

  call messages_init()

  call io_init()
  call unit_system_init(global_namespace)

  call center_geo()

  call io_end()
  call messages_end()

  call parser_end()

  call global_end()

contains

  subroutine center_geo()

    integer, parameter :: &
      NONE    = 0,        &
      INERTIA = 1,        &
      PSEUDO  = 2,        &
      LARGE   = 3

    type(geometry_t)  :: geo
    type(space_t)     :: space
    FLOAT :: center(MAX_DIM), x1(MAX_DIM), x2(MAX_DIM), to(MAX_DIM)
    integer :: axis_type, idir, default
    type(block_t) :: blk

    call space_init(space, global_namespace)
    call geometry_init(geo, global_namespace, space)

    ! is there something to do
    if (geo%natoms > 1) then

      !%Variable MainAxis
      !%Type block
      !%Section Utilities::oct-center-geom
      !%Description 
      !% A vector of reals defining the axis to which the molecule
      !% should be aligned. If not present, the default value will
      !% be the x-axis. For example in 3D:
      !% <tt>
      !% <br>%MainAxis
      !% <br> 1 | 0 | 0 
      !% <br>%</tt>
      !%End
      if(parse_block(global_namespace, 'MainAxis', blk)==0) then
        do idir = 1, geo%space%dim
          call parse_block_float(blk, 0, idir - 1, to(idir))
        end do
        call parse_block_end(blk)
      else
        to(:) = M_ZERO
        to(1) = M_ONE
      end if
      to = to / sqrt(sum(to(1:geo%space%dim)**2))

      write(message(1),'(a,6f15.6)') 'Using main axis ', to(1:geo%space%dim)
      call messages_info(1)

      !%Variable AxisType
      !%Type integer
      !%Default inertia
      !%Section Utilities::oct-center-geom
      !%Description
      !% After the structure is centered, it is also aligned to a set of orthogonal axes.
      !% This variable decides which set of axes to use. Only implemented for 3D, in which case
      !% the default is <tt>inertia</tt>; otherwise <tt>none</tt> is default and the only legal value.
      !%Option none 0
      !% Do not rotate. Will still give output regarding center of mass and moment of inertia.
      !%Option inertia 1
      !% The axis of inertia.
      !%Option pseudo_inertia 2
      !% Pseudo-axis of inertia, calculated considering all species to have equal mass.
      !%Option large_axis 3
      !% The larger axis of the molecule.
      !%End
      if (geo%space%dim == 3) then
        default = INERTIA
      else
        default = NONE
      end if
      call parse_variable(global_namespace, 'AxisType', default, axis_type)
      call messages_print_var_option(stdout, "AxisType", axis_type)

      if (geo%space%dim /= 3 .and. axis_type /= NONE) then
        call messages_not_implemented("alignment to axes (AxisType /= none) in other than 3 dimensions", namespace=global_namespace)
      end if

      select case (axis_type)
      case (NONE, INERTIA, PSEUDO)
        center = geometry_center_of_mass(geo, pseudo = (axis_type==PSEUDO))

        write(message(1),'(3a,99f15.6)') 'Center of mass [', trim(units_abbrev(units_out%length)), '] = ', &
          (units_from_atomic(units_out%length, center(idir)), idir = 1, geo%space%dim)
        call messages_info(1)

        call geometry_translate(geo, center)
        call geometry_axis_inertia(geo, x1, x2, pseudo = (axis_type==PSEUDO))
      case (LARGE)
        center = geometry_center(geo)

        write(message(1),'(3a,99f15.6)') 'Center [', trim(units_abbrev(units_out%length)), '] = ', &
          (units_from_atomic(units_out%length, center(idir)), idir = 1, geo%space%dim)
        call messages_info(1)

        call geometry_translate(geo, center)
        call geometry_axis_large(geo, x1, x2)
      case default
        write(message(1), '(a,i2,a)') 'AxisType = ', axis_type, ' not known by Octopus.'
        call messages_fatal(1, namespace=global_namespace)
      end select

      write(message(1),'(a,99f15.6)') 'Found primary   axis ', x1(1:geo%space%dim)
      write(message(2),'(a,99f15.6)') 'Found secondary axis ', x2(1:geo%space%dim)
      call messages_info(2)

      if (axis_type /= NONE) call geometry_rotate(geo, global_namespace, x1, x2, to)

    end if

    ! recenter
    center = geometry_center(geo)
    call geometry_translate(geo, center)

    ! write adjusted geometry
    call geometry_write_xyz(geo, './adjusted', global_namespace)

    call geometry_end(geo)

  end subroutine center_geo

end program centergeom

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
