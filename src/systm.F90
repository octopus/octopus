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

module system
use io
use specie
use atom
use mesh
use states
use output

implicit none

type system_type
  character(len=20) :: sysname ! the name of the system we are running

  integer :: natoms
  type(atom_type), pointer :: atom(:)
  integer :: ncatoms  ! For QM+MM calculations
  type(atom_classical_type), pointer :: catom(:)

  real(r8) :: eii, kinetic_energy ! the ion-ion energy

  integer :: nspecies
  type(specie_type), pointer :: specie(:)
  logical :: nlpp    ! is any species having non-local pp
  logical :: nlcc    ! is any species having non-local core corrections?

  ! the charge of the valence electrons (necessary to initialize states)
  real(r8) :: val_charge 

  type(mesh_type) :: m
  type(states_type), pointer :: st

  type(output_type) :: outp
end type system_type

contains

subroutine system_init(s)
  type(system_type), intent(out) :: s

  integer :: i

  call push_sub('system_init')

  call oct_parse_string('SystemName', 'system', s%sysname)
  s%nspecies = specie_init(s%specie)

  call atom_init(s%natoms, s%atom, s%ncatoms, s%catom, s%nspecies, s%specie)

  call output_init(s%outp)

  call mesh_init(s%m, s%natoms, s%atom)
  call functions_init(s%m)

  !  find total charge of the system
  s%val_charge = 0
  do i = 1, s%natoms
    s%val_charge = s%val_charge - s%atom(i)%spec%Z_val
  enddo

  allocate(s%st)
  call states_init(s%st, s%m, s%val_charge)

  ! find out if we need non-local core corrections
  s%nlcc = .false.
  s%nlpp = .false.
  do i = 1, s%nspecies
    s%nlcc = (s%nlcc.or.s%specie(i)%nlcc)
    s%nlpp = (s%nlcc.or.(.not.s%specie(i)%local))
  end do
  if(s%nlcc) allocate(s%st%rho_core(s%m%np))

  call pop_sub()
end subroutine system_init

subroutine system_end(s)
  type(system_type), intent(inout) :: s

  call push_sub('system_end')

  if(associated(s%st)) then
    call states_end(s%st)
    deallocate(s%st); nullify(s%st)
  end if
  call functions_end(s%m)
  call mesh_end(s%m)
  call atom_end(s%natoms, s%atom, s%ncatoms, s%catom)
  call specie_end(s%nspecies, s%specie)

  call pop_sub()
end subroutine system_end

end module system
