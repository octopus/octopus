!! Copyright (C) 2016 N. Tancogne-Dejean 
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
!! $Id$

#include "global.h"

module orbital_oct_m
  use global_oct_m
  use messages_oct_m  
  use profiling_oct_m
  use types_oct_m
 
  implicit none

  private

  public ::                        &
             orbital_t,            &
             orbital_nullify,      &
             orbital_init,         &
             orbital_end

  type orbital_t
    FLOAT, pointer      :: dorb(:) !> The orbital, if real, on the submesh
    CMPLX, pointer      :: zorb(:) !> The orbital, if complex, on the submesh
    CMPLX, pointer      :: eorb(:,:) !> Orbitals with its phase factor
  end type orbital_t

contains

 subroutine orbital_nullify(this)
  type(orbital_t),             intent(inout) :: this

  PUSH_SUB(orbital_nullify)

  nullify(this%dorb)
  nullify(this%zorb)
  nullify(this%eorb)

  POP_SUB(orbital_nullify)

 end subroutine orbital_nullify

 subroutine orbital_init(this)
  type(orbital_t),             intent(inout) :: this

  PUSH_SUB(orbital_init)


  POP_SUB(orbital_init)
 end subroutine orbital_init


 subroutine orbital_end(this)
   type(orbital_t), intent(inout) :: this

   PUSH_SUB(orbital_end)  

   SAFE_DEALLOCATE_P(this%dorb)
   SAFE_DEALLOCATE_P(this%zorb)
   SAFE_DEALLOCATE_P(this%eorb)

   POP_SUB(orbital_end)
 end subroutine orbital_end

end module orbital_oct_m
