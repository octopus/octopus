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

  integer :: ik, i, k
  real(r8) :: l(3)

  allocate(st%kpoints(3, st%nik), st%kweights(st%nik))
  st%kpoints = M_ZERO

  if(st%nik == 1.or.(st%nik == 2.and.st%ispin == 2)) then ! just return the Gamma point
    st%kweights(:) = M_ONE
    return
  end if

  l=M_ZERO
  do ik = 1, st%nik
    if(st%ispin == 2) then
      k = (ik-1)/2
    else
      k = ik - 1
    end if
    do i=1,conf%periodic_dim
        l(i) = M_PI/m%lsize(i)
        st%kpoints(i, ik) = -l(i)*k/(st%nik - 1)
    end do
! this does not take into accout the symmetries
! it is correct only for simple cubic with center of inversion in point group
    st%kweights(ik) = M_ONE/real(st%nik, r8)
  end do

end subroutine states_choose_kpoints

subroutine kpoints_write_info(st,iunit)
  
  type(states_type), intent(IN) :: st
  integer, intent(IN) :: iunit
  
  
  write(message(1),'(a,i4)') 'Number of K points in each direction = ',st%nik
  call write_info(1,iunit)

end subroutine kpoints_write_info
