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

program atoms_magnet
  use global
  use lib_oct_parser
  use lib_oct
  use units
  use mesh
  use mesh_function
  use functions
  use geometry
  use output
  use states
  use restart

  implicit none

  type(mesh_type), pointer      :: m          ! the mesh
  type(f_der_type)              :: f_der      ! manages the calculation of derivatives of functions
  type(geometry_type), pointer  :: geo        ! the geometry
  type(states_type), pointer    :: st         ! the states
  integer :: i, ia, j, iunit, err
  FLOAT :: val_charge, def_h, def_rsize, r, rmin
  FLOAT, allocatable :: rho(:,:), magnet(:,:)

  call init_()

  ! load information
  call X(restart_read)('tmp/restart_gs', st, m, err)
  if(err.ne.0) then
    message(1) = "Error opening states in 'tmp/restart_gs'"
    call write_fatal(1)
  end if

  ! get default radius for the spheres
  call geometry_min_distance(geo, rmin)
  call loct_parse_float("MagnetizationSpheresRadius", rmin*M_HALF/units_inp%length%factor, r)
  r = r * units_inp%length%factor

  ! open file for output
  iunit = io_open('static/atoms_magnet')

  write(iunit, '(a)') geo%sysname
  write(iunit,'(/a)') "Magnetization"
  select case (st%d%ispin)
  case (SPIN_POLARIZED)
    write(iunit,'(a,7x,13x,a)') ' Ion','mz'
  case (SPINORS)
    write(iunit,'(a,7x,13x,a,13x,a,13x,a)') ' Ion','mx','my','mz'
  end select

  !get density
  allocate(rho(m%np, st%d%nspin))
  call X(states_calc_dens)(st, m%np, rho, reduce=.true.)

  !get atomic magnetizations
  allocate(magnet(3, geo%natoms))
  call states_local_magnetic_moments(m, st, geo, rho, r, magnet)
  deallocate(rho)

  do ia = 1, geo%natoms
    select case (st%d%ispin)
    case (SPIN_POLARIZED)
      write(iunit,'(i4,a10,f15.6)') ia, trim(geo%atom(ia)%spec%label), magnet(3, ia)
    case (SPINORS)
      write(iunit,'(i4,a10,3f15.6)') ia, trim(geo%atom(ia)%spec%label), magnet(1, ia), magnet(2, ia), magnet(3, ia)
    end select
  end do

  deallocate(magnet)
  call io_close(iunit)
  call end_()

contains
  ! Initialize stuff
  subroutine init_()
    call global_init()
    call units_init()

    allocate(geo, m)
    
    call geometry_init_xyz(geo)
    call geometry_init_species(geo, val_charge, def_h, def_rsize)
    
    call mesh_init(m, geo, def_h, def_rsize)
    call f_der_init(m, f_der)
    call mesh_create_xyz(m, f_der%n_ghost(1))
    call f_der_build(f_der)
    call mesh_write_info(m, stdout)
    
    allocate(st)
    call states_init(st, m, geo, val_charge, geo%nlcc)
    if (st%d%ispin == UNPOLARIZED) then
      message(1) = "You cannot use this utility with spin-unpolarized calculations"
      call write_fatal(1)
    end if
    allocate(st%X(psi)(m%np, st%d%dim, st%nst, st%d%nik))

  end subroutine init_

  subroutine end_()
    call states_end(st)
    deallocate(st); nullify(st)
      
    call geometry_end(geo)
    call f_der_end(f_der)
    call mesh_end(m)

    deallocate(geo, m)
    
    call global_end()
    
  end subroutine end_

end program atoms_magnet
