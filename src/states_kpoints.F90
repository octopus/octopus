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

! We assume time-reversal symmetry, i.e. k = -k
subroutine states_choose_kpoints(st, m)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m

  integer :: ik, k
  real(r8) :: l

  allocate(st%kpoints(3, st%nik), st%kweights(st%nik))
  st%kpoints = 0._r8

  if(st%nik == 1.or.(st%nik == 2.and.st%ispin == 2)) then ! just return the Gamma point
    st%kweights(:) = 1._r8
    return
  end if

  l = M_PI/m%h(3)
  do ik = 1, st%nik
    if(st%ispin == 2) then
      k = (ik-1)/2
    else
      k = ik - 1
    end if

    st%kpoints(3, ik) = -l*k/(st%nik - 1)
    st%kweights(ik) = 1._r8/real(st%nik, r8)
  end do
end subroutine states_choose_kpoints
