!! Copyright (C) 2010 X. Andrade
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
!! $Id: liquid.F90 6346 2010-03-08 22:05:37Z acastro $

#include "global.h"

program centergeom
  use command_line_m
  use c_pointer_m
  use datasets_m
  use geometry_m
  use global_m
  use io_m
  use messages_m
  use loct_math_m
  use parser_m
  use profiling_m
  use unit_m
  use unit_system_m
  use xyz_adjust_m

  implicit none

  integer :: ierr
  type(geometry_t) :: geo
  integer, parameter :: dim = 3

  ! this is a hack that must be fixed
  calc_dim = 3

  call global_init()                       ! initialize
  call getopt_init(ierr)
  if(ierr.eq.0) call getopt_center_geom
  call parser_init()
  call parse_integer('DebugLevel', 0, conf%debug_level)
  if(conf%debug_level>0) then
    in_debug_mode = .true.
  else
    in_debug_mode = .false.
  end if

  call datasets_init(1)
  call io_init()

  call unit_system_init()

  call geometry_init(geo)

  call generate_liquid()

  call geometry_end(geo)

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()

contains
  
  subroutine generate_liquid()
    FLOAT        :: density, cell_size(1:dim), center(1:dim)
    integer      :: nmolecules
    type(unit_t) :: gr_per_cm3
    type(c_ptr)  :: random_gen_pointer
    integer      :: iatom, jatom, imolecule, idir, iunit
    FLOAT, allocatable  :: coordinates(:, :)

    ! we read the density in sensible units
    gr_per_cm3%factor = CNST(10.541009)
    gr_per_cm3%abbrev = 'gr/cm^3'
    gr_per_cm3%name   = 'gr/cm^3'

    !%Variable LiquidDensity
    !%Type float
    !%Section Utilities::oct-liquid
    !%Description
    !% This variable specifies the density of the liquid to be
    !% generated. It has to be given in units of grammes per cubic
    !% centimeters <math>gr/cm^3</math>.
    !%End
    call parse_float('LiquidDensity', CNST(-1.0), density, gr_per_cm3)

    if(density <= M_ZERO) call input_error('LiquidDensity')

    write(message(1), '(a, f10.3, a)') 'Info: Liquid density = ', units_from_atomic(gr_per_cm3, density), ' '//units_abbrev(gr_per_cm3)
    call write_info(1)

    call xyz_adjust_it(geo)

    !%Variable LiquidNumberOfMolecules
    !%Type float
    !%Section Utilities::oct-liquid
    !%Description
    !% This variable specifies the number of molecules the will be
    !% included in the liquid.
    !%End
    call parse_integer('LiquidNumberOfMolecules', -1, nmolecules)
    if(nmolecules <= 0) call input_error('LiquidNumberOfMolecules')

    write(message(1), '(a, i5, a)') 'Info: ', nmolecules, ' molecules will be included in the liquid.'
    call write_info(1)

    ! we assume a cubic cell
    cell_size = (geometry_mass(geo)*nmolecules/density)**(1.0/3.0)

    write(message(1), '(a, f10.5, a)') 'Info: cubic cell size ', units_from_atomic(units_out%length, cell_size(1)), ' '//units_abbrev(units_out%length)
    call write_info(1)    

    ! Now comes the real part, the calculation of the coordinates. For
    ! the moment this we don't avoid molecule colisions. We don't
    ! rotate the molecules either.

    SAFE_ALLOCATE(coordinates(1:dim, nmolecules*geo%natoms)) 

    call loct_ran_init(random_gen_pointer)

    jatom = 1
    do imolecule = 1, nmolecules
      
      do idir = 1, dim
        center(idir) = loct_ran_flat(random_gen_pointer, -M_HALF*cell_size(idir), M_HALF*cell_size(idir))
      end do
      do iatom = 1, geo%natoms
        do idir = 1, dim
          coordinates(idir, jatom) = center(idir) + geo%atom(iatom)%x(idir)
        end do
        INCR(jatom, 1)
      end do

    end do
    
    call loct_ran_end(random_gen_pointer)

    ! now print the coordinates to a xyz file

    iunit = io_open('liquid.xyz', action='write')

    write(iunit, '(i4)') nmolecules*geo%natoms
    write(iunit, '(3a,3f20.5)') '# units: ', trim(units_abbrev(units_out%length)), '    cell size:',  cell_size(1:3)
    jatom = 1
    do imolecule = 1, nmolecules
      do iatom = 1, geo%natoms
        write(iunit, '(6x,a,2x,99f12.6)') geo%atom(iatom)%label, &
          (units_from_atomic(units_out%length, coordinates(idir, jatom)), idir = 1, dim)
        INCR(jatom, 1)
      end do
    end do
    call io_close(iunit)

    write(message(1), '(a)') 'Info: the atomic coordinates have written to a file called liquid.xyz'
    call write_info(1)   

  end subroutine generate_liquid
end program centergeom

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
