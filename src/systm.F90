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

module system
use global
use mesh
use geometry
use states

implicit none

type system_type
  FLOAT                      :: val_charge ! the charge of the valence electrons (necessary to initialize states)
  type(mesh_type)            :: m          ! the mesh
  type(geometry_type)        :: geo        ! the geometry
  type(states_type), pointer :: st         ! the states
  type(output_type)          :: outp
end type system_type

contains

subroutine system_init(s)
  type(system_type), intent(out) :: s

  call push_sub('system_init')

  ! initialize the stuff related to the mesh
  call mesh_init(s%m)
  call functions_init(s%m)

  ! initialize the other stuff
  call geometry_init(s%geo, s%val_charge)
  allocate(s%st)
  call states_init(s%st, s%m, s%val_charge, s%geo%nlcc)
  call output_init(s%outp)

  call pop_sub()
end subroutine system_init

subroutine system_end(s)
  type(system_type), intent(inout) :: s

  call push_sub('system_end')

  if(associated(s%st)) then
    call states_end(s%st)
    deallocate(s%st); nullify(s%st)
  end if
  call geometry_end(s%geo)
  call functions_end(s%m)
  call mesh_end(s%m)

  call pop_sub()
end subroutine system_end

end module system
