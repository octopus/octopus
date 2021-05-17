!! Copyright (C) 2015 X. Andrade
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

module partial_charges_oct_m
  use global_oct_m
  use hirshfeld_oct_m
  use ions_oct_m
  use mesh_oct_m
  use messages_oct_m
  use profiling_oct_m
  use states_elec_oct_m

  implicit none

  private

  public ::                             &
    partial_charges_calculate

contains

  !----------------------------------------------
  subroutine partial_charges_calculate(mesh, st, ions, hirshfeld_charges)
    type(mesh_t),            intent(in)    :: mesh
    type(states_elec_t),     intent(in)    :: st
    type(ions_t),            intent(in)    :: ions
    FLOAT, optional,         intent(out)   :: hirshfeld_charges(:)

    integer :: iatom
    type(profile_t), save :: prof
    type(hirshfeld_t) :: hirshfeld
    
    PUSH_SUB(partial_charges_calculate)
    call profiling_in(prof, 'PARTIAL_CHARGES')

    if(present(hirshfeld_charges)) then

      call hirshfeld_init(hirshfeld, mesh, ions, st)
      
      do iatom = 1, ions%natoms
        call hirshfeld_charge(hirshfeld, iatom, st%rho, hirshfeld_charges(iatom))
      end do
      
      call hirshfeld_end(hirshfeld)
    end if
    
    call profiling_out(prof)
    POP_SUB(partial_charges_calculate)

  end subroutine partial_charges_calculate

end module partial_charges_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
