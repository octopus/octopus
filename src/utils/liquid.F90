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

program liquid
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
  use space_m
  use unit_m
  use unit_system_m
  use xyz_adjust_m

  implicit none

  integer :: ierr
  type(geometry_t) :: geo
  type(space_t)    :: space
  integer, parameter :: dim = 3

  call global_init()                       ! initialize
  call getopt_init(ierr)
  if(ierr.eq.0) call getopt_center_geom
  call parser_init()
  call messages_init()

  call datasets_init(1)
  call io_init()

  call unit_system_init()

  call space_init(space)
  call geometry_init(geo, space)

  call generate_liquid()

  call geometry_end(geo)

  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()

contains
  
  subroutine generate_liquid()
    FLOAT        :: density, cell_size(1:dim), axis(1:MAX_DIM, 1:3)
    FLOAT        :: ang1, ang2, molecule_min_dist, min_dist, scale, radius, dist
    integer      :: nmolecules, dx, dy, dz
    type(unit_t) :: gr_per_cm3
    type(c_ptr)  :: random_gen_pointer
    integer      :: iatom, jatom, imolecule, jmolecule, idir, iunit, iaxis, iter, niter
    FLOAT, allocatable  :: coordinates(:, :), center(:, :)
    logical      :: too_close

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

    write(message(1), '(a, f10.3, a)') 'Info: liquid density = ', &
      units_from_atomic(gr_per_cm3, density), ' '//units_abbrev(gr_per_cm3)
    call messages_info(1)

    call xyz_adjust_it(geo)

    !%Variable LiquidNumberOfMolecules
    !%Type integer
    !%Section Utilities::oct-liquid
    !%Description
    !% This variable specifies the number of molecules that will be
    !% included in the liquid.
    !%End
    call parse_integer('LiquidNumberOfMolecules', -1, nmolecules)
    if(nmolecules <= 0) call input_error('LiquidNumberOfMolecules')

    write(message(1), '(a, i5, a)') 'Info: ', nmolecules, ' molecules will be included in the liquid.'
    call messages_info(1)

    !%Variable LiquidMoleculeScale
    !%Type float
    !%Section Utilities::oct-liquid
    !%Description
    !% This value specifies how much the original molecule will be
    !% scaled in the liquid. The default is 1.0.
    !%End
    call parse_float('LiquidMoleculeScale', CNST(1.0), scale)

    radius = M_ZERO
    do iatom = 1, geo%natoms
      do idir = 1, dim
        geo%atom(iatom)%x(idir) = scale*geo%atom(iatom)%x(idir)
        radius = max(radius, geo%atom(iatom)%x(idir))
      end do
    end do

    molecule_min_dist = geometry_min_distance(geo)

    radius = radius + CNST(0.5)*molecule_min_dist

    ! we assume a cubic cell
    cell_size = (geometry_mass(geo)*nmolecules/density)**(1.0/3.0)

    write(message(1), '(a, f10.5, a)') 'Info: cubic cell size ', &
      units_from_atomic(units_out%length, cell_size(1)), ' '//units_abbrev(units_out%length)
    call messages_info(1)    

    ! Now comes the real part, the calculation of the coordinates.

    SAFE_ALLOCATE(coordinates(1:dim, nmolecules*geo%natoms))
    SAFE_ALLOCATE(center(1:dim, nmolecules)) 

    call loct_ran_init(random_gen_pointer)

    axis = M_ZERO
    niter = 10
    jatom = 1
    do imolecule = 1, nmolecules

      ! randomly select a point where the molecule will be placed
      do iter = 1, niter
        do idir = 1, dim
          center(idir, imolecule) = loct_ran_flat(random_gen_pointer, -M_HALF*cell_size(idir), M_HALF*cell_size(idir))
        end do

        ! we avoid points too close to the previous ones
        too_close = .false.

        do jmolecule = 1, imolecule - 1

          do dx = -1, 1
            do dy = -1, 1
              do dz = -1, 1

                dist = sum((center(:, imolecule) - center(:, jmolecule) + (/dx, dy, dz/)*cell_size(1:dim))**2)
                if(dist <= M_FOUR*radius**2) then
                  too_close = .true.
                  exit
                end if

              end do
            end do
          end do

        end do

        if(.not. too_close) exit
      end do

      ! and some rotation axes
      do iaxis = 1, 3
        ang1 = loct_ran_flat(random_gen_pointer, M_ZERO, M_PI)
        ang2 = loct_ran_flat(random_gen_pointer, M_ZERO, CNST(2.0)*M_PI)
        axis(1, iaxis) = sin(ang1)*cos(ang2)
        axis(2, iaxis) = sin(ang1)*sin(ang2)
        axis(3, iaxis) = cos(ang1)
      end do

      call geometry_rotate(geo, axis(:, 1), axis(:, 2), axis(:, 3))

      do iatom = 1, geo%natoms
        do idir = 1, dim
          coordinates(idir, jatom) = center(idir, imolecule) + geo%atom(iatom)%x(idir)

          ! now force the atoms to be in the cell
          if(coordinates(idir, jatom) < -M_HALF*cell_size(idir)) &
            coordinates(idir, jatom) = coordinates(idir, jatom) + cell_size(idir)

          if(coordinates(idir, jatom) >= M_HALF*cell_size(idir)) &
            coordinates(idir, jatom) = coordinates(idir, jatom) - cell_size(idir)

        end do
        INCR(jatom, 1)
      end do

    end do

    !check the distance between atoms
    min_dist = HUGE(min_dist)

    do iatom = 1, geo%natoms*nmolecules
      do jatom = 1, geo%natoms*nmolecules
        if(iatom == jatom) cycle

        do dx = -1, 1
          do dy = -1, 1
            do dz = -1, 1
              dist = sum((coordinates(1:dim, iatom) - coordinates(1:dim, jatom) + (/dx, dy, dz/)*cell_size(1:dim))**2)
              min_dist = min(min_dist, dist)
            end do
          end do
        end do

      end do
    end do

    min_dist = sqrt(min_dist)

    write(message(1), '(a, f10.5, a)') 'Info: minimum distance between atoms in the liquid :', &
      units_from_atomic(units_out%length, min_dist), ' '//trim(units_abbrev(units_out%length))
    call messages_info(1)

    SAFE_DEALLOCATE_A(center)

    call loct_ran_end(random_gen_pointer)

    ! now print the coordinates to a xyz file

    iunit = io_open('liquid.xyz', action='write')

    write(iunit, '(i4)') nmolecules*geo%natoms
    write(iunit, '(3a,3f20.5)') '# units: ', trim(units_abbrev(units_out%length)), &
      '    cell size:',  (units_from_atomic(units_out%length, cell_size(idir)), idir = 1, 3)
    jatom = 1
    do imolecule = 1, nmolecules
      do iatom = 1, geo%natoms
        write(iunit, '(6x,a,2x,99f12.6)') geo%atom(iatom)%label, &
          (units_from_atomic(units_out%length, coordinates(idir, jatom)), idir = 1, dim)
        INCR(jatom, 1)
      end do
    end do
    call io_close(iunit)

    write(message(1), '(a)') "Info: the atomic coordinates have been written to 'liquid.xyz'."
    call messages_info(1)   

    SAFE_DEALLOCATE_A(coordinates)

  end subroutine generate_liquid
end program liquid

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
