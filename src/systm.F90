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

  type(mesh_type) :: m
  type(states_type), pointer :: st

  type(output_type) :: outp
end type system_type

contains

subroutine system_init(s)
  type(system_type), intent(out) :: s

  integer :: i
  real(r8) :: val_charge

  sub_name = 'system_init'; call push_sub()

  call oct_parse_str('SystemName', 'system', s%sysname)
  s%nspecies = specie_init(s%specie)

  call atom_init(s%natoms, s%atom, s%ncatoms, s%catom, s%nspecies, s%specie)

  call output_init(s%outp)

  call mesh_init(s%m, s%natoms, s%atom)

  !  find total charge of the system
  val_charge = 0
  do i = 1, s%natoms
    val_charge = val_charge - s%atom(i)%spec%Z_val
  enddo

  allocate(s%st)
  call states_init(s%st, s%m, val_charge)

  ! find out if we need non-local core corrections
  s%st%nlcc = .false.
  do i = 1, s%nspecies
    s%st%nlcc = (s%st%nlcc.or.s%specie(i)%nlcc)
  end do
  if(s%st%nlcc) allocate(s%st%rho_core(s%m%np))

  call pop_sub()
end subroutine system_init

subroutine system_end(s)
  type(system_type), intent(inout) :: s

  sub_name = 'system_end'; call push_sub()

  if(associated(s%st)) then
    call states_end(s%st)
    deallocate(s%st); nullify(s%st)
  end if
  call mesh_end(s%m)
  call atom_end(s%natoms, s%atom, s%ncatoms, s%catom)
  call specie_end(s%nspecies, s%specie)

  call pop_sub()
end subroutine system_end

subroutine system_write_xyz(dir, fname, s)
  character(len=*), intent(in)  :: dir, fname
  type(system_type), intent(IN) :: s

  integer i, iunit
  
  ! xyz format, for easy plot in rasmol
#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    call oct_mkdir(C_string(trim(dir)))

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+'.xyz', status='unknown')
    write(iunit, '(i4)') s%natoms
    write(iunit, '(1x)')
    do i = 1, s%natoms
      write(iunit, '(6x,a,2x,3f12.6)') &
           s%atom(i)%spec%label, s%atom(i)%x(:)/units_out%length%factor
    end do
    call io_close(iunit)
    
    if(s%ncatoms > 0) then
      call io_assign(iunit)
      open(iunit, file=trim(dir)+"/"+trim(fname)+'_classical.xyz', status='unknown')
      write(iunit, '(i4)') s%ncatoms
      write(iunit, '(1x)')
      do i = 1, s%ncatoms
        write(iunit, '(6x,a1,2x,3f12.6,a,f12.6)') &
             s%catom(i)%label(1:1), s%catom(i)%x(:)/units_out%length%factor, &
             " # ", s%catom(i)%charge
      end do
      call io_close(iunit)
    end if

#ifdef HAVE_MPI
  end if
#endif
  
  return
end subroutine system_write_xyz

end module system
