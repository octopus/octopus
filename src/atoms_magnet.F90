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
  integer :: ik, ist, i, ia, j, iunit, err
  R_TYPE :: c
  FLOAT :: val_charge, def_h, def_rsize, mg(3), r, rmin, ri

  ! Initialize stuff
  call global_init()
  call units_init()

  allocate(geo)
  call geometry_init_xyz(geo)
  call geometry_init_species(geo, val_charge, def_h, def_rsize)
  allocate(m)
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
  allocate(st%X(psi)(m%np, st%dim, st%nst, st%nik))

  ! load information
  call X(restart_read)("tmp/restart_gs", st, m, err)
  if(err.ne.0) then
    message(1) = "Error opening states in 'tmp/restart_gs'"
    call write_fatal(1)
  end if

  ! get default radius for the spheres
  rmin = huge(PRECISION)
  do i = 1, geo%natoms
    do j = i+1, geo%natoms
      r = sqrt(sum((geo%atom(i)%x-geo%atom(j)%x)**2))
      if(r < rmin) then
        rmin = r
      end if
    end do
  end do
  rmin = M_HALF * rmin

  call loct_parse_float("MagnetizationSpheresRadius", rmin/units_inp%length%factor, r)
  r = r * units_inp%length%factor

  ! open file for output
  call io_assign(iunit)
  open(iunit, status='unknown', file="atoms_magnet")
  write(iunit, '(a)') geo%sysname
  write(iunit,'(/a)') "Magnetization"
  select case (st%d%ispin)
  case (SPIN_POLARIZED)
    write(iunit,'(a,7x,13x,a)') ' Ion','mx'
  case (SPINORS)
    write(iunit,'(a,7x,13x,a,13x,a,13x,a)') ' Ion','mx','my','mz'
  end select

  ! compute the local magnetization around each atom 
  do ia = 1, geo%natoms
    mg = M_ZERO
    do ik = 1, st%nik
      do ist = 1, st%nst
        do i = 1, m%np
          call mesh_r(m, i, ri, a=geo%atom(ia)%x)
          if (ri > r) cycle
          select case (st%d%ispin)
          case (SPIN_POLARIZED)
            c = R_CONJ(st%X(psi)(i, 1, ist, ik)) * st%X(psi)(i, 1, ist, ik) * m%vol_pp(i)
            if(modulo(ik+1, 2) == 0) then
              mg(1) = mg(1) - st%d%kweights(ik)*st%occ(ist, ik)* M_TWO*R_REAL(c)
            else
              mg(1) = mg(1) + st%d%kweights(ik)*st%occ(ist, ik)* M_TWO*R_REAL(c)
            end if
          case (SPINORS)
            c = st%X(psi)(i, 1, ist, ik) * R_CONJ(st%X(psi)(i, 2, ist, ik)) * m%vol_pp(i)
            mg(1) = mg(1) + st%d%kweights(ik)*st%occ(ist, ik)* M_TWO*R_REAL(c)
            mg(2) = mg(2) - st%d%kweights(ik)*st%occ(ist, ik)* M_TWO*R_AIMAG(c)
            c = R_ABS(st%X(psi)(i, 1, ist, ik))**2 - R_ABS(st%X(psi)(i, 2, ist, ik))**2 * m%vol_pp(i)
            mg(3) = mg(3) + st%d%kweights(ik)*st%occ(ist, ik)* R_REAL(c)
          end select
        end do
      end do
    end do

    select case (st%d%ispin)
    case (SPIN_POLARIZED)
      write(iunit,'(i4,a10,f15.6)') ia, trim(geo%atom(ia)%spec%label), mg(1)
    case (SPINORS)
      write(iunit,'(i4,a10,3f15.6)') ia, trim(geo%atom(ia)%spec%label), mg(1), mg(2), mg(3)
    end select
  end do

  ! end stuff
  call io_close(iunit)
  if(associated(st)) then
    call states_end(st)
    deallocate(st); nullify(st)
  end if
  call geometry_end(geo)
  call mesh_end(m)
  deallocate(geo)

end program atoms_magnet
