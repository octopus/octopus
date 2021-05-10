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
  use global_oct_m
  use io_oct_m
  use ions_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
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

    type(ions_t),     pointer :: ions
    FLOAT, allocatable :: center(:), x1(:), x2(:), to(:)
    integer :: axis_type, idir, default
    type(block_t) :: blk

    ions => ions_t(global_namespace)

    SAFE_ALLOCATE(center(1:ions%space%dim))
    SAFE_ALLOCATE(x1(1:ions%space%dim))
    SAFE_ALLOCATE(x2(1:ions%space%dim))
    SAFE_ALLOCATE(to(1:ions%space%dim))

    ! is there something to do
    if (ions%natoms > 1) then

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
      if(parse_block(ions%namespace, 'MainAxis', blk)==0) then
        do idir = 1, ions%space%dim
          call parse_block_float(blk, 0, idir - 1, to(idir))
        end do
        call parse_block_end(blk)
      else
        to(:) = M_ZERO
        to(1) = M_ONE
      end if
      to = to / norm2(to)

      write(message(1),'(a,6f15.6)') 'Using main axis ', to
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
      if (ions%space%dim == 3) then
        default = INERTIA
      else
        default = NONE
      end if
      call parse_variable(ions%namespace, 'AxisType', default, axis_type)
      call messages_print_var_option(stdout, "AxisType", axis_type)

      if (ions%space%dim /= 3 .and. axis_type /= NONE) then
        call messages_not_implemented("alignment to axes (AxisType /= none) in other than 3 dimensions", namespace=ions%namespace)
      end if

      select case (axis_type)
      case (NONE, INERTIA, PSEUDO)
        center = ions%center_of_mass(pseudo = (axis_type==PSEUDO))

        write(message(1),'(3a,99f15.6)') 'Center of mass [', trim(units_abbrev(units_out%length)), '] = ', &
          units_from_atomic(units_out%length, center)
        call messages_info(1)

        call ions%translate(center)
        call ions%axis_inertia(x1, x2, pseudo = (axis_type==PSEUDO))
      case (LARGE)
        center = ions%center()

        write(message(1),'(3a,99f15.6)') 'Center [', trim(units_abbrev(units_out%length)), '] = ', &
          units_from_atomic(units_out%length, center)
        call messages_info(1)

        call ions%translate(center)
        call ions%axis_large(x1, x2)
      case default
        write(message(1), '(a,i2,a)') 'AxisType = ', axis_type, ' not known by Octopus.'
        call messages_fatal(1, namespace=ions%namespace)
      end select

      write(message(1),'(a,99f15.6)') 'Found primary   axis ', x1
      write(message(2),'(a,99f15.6)') 'Found secondary axis ', x2
      call messages_info(2)

      if (axis_type /= NONE) then
        call ions%rotate(x1, x2, to)
      end if

    end if

    ! recenter
    center = ions%center()
    call ions%translate(center)

    ! write adjusted geometry
    call ions%write_xyz('./adjusted')

    SAFE_DEALLOCATE_A(center)
    SAFE_DEALLOCATE_A(x1)
    SAFE_DEALLOCATE_A(x2)
    SAFE_DEALLOCATE_A(to)
    SAFE_DEALLOCATE_P(ions)

  end subroutine center_geo

end program centergeom

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
