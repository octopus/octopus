!! Copyright (C) 2017 N. Tancogne-Dejean 
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


  !At the present time this routine can only return atomic orbitals, but this could be generalized
  subroutine X(orbitalset_utils_getorbitals)(os, ions, mesh, use_mesh, normalize)
    type(orbitalset_t),      intent(inout) :: os
    type(ions_t),            intent(in)    :: ions
    type(mesh_t),            intent(in)    :: mesh
    logical,                 intent(in)    :: use_mesh
    logical,                 intent(in)    :: normalize

    integer :: iorb

    PUSH_SUB(X(orbitalset_utils_getorbitals))

    do iorb = 1, os%norbs
      if (debug%info) then
        write(message(1),'(a,i3,1x,i1,1x,i1,1x,i1,1x,f3.1)')  'get_atomic_orbital ', os%iatom, &
                                                                    iorb, os%ii, os%ll, os%jj
        call messages_info(1)
      end if
      ! We obtain the orbital
      call X(get_atomic_orbital)(ions, mesh, os%sphere, os%iatom, os%ii, os%ll, os%jj, &
                                              os, iorb, os%radius, os%ndim, use_mesh, normalize)
    end do !iorb

    POP_SUB(X(orbitalset_utils_getorbitals))

   end subroutine X(orbitalset_utils_getorbitals)

