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

subroutine states_choose_kpoints(st, m)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m

  logical  :: coi
  integer  :: i, ik, jk, kk, j, k, skip, tmp_nik_axis(3)
  FLOAT :: l(3), total_weight

  allocate(st%kpoints(3, st%nik), st%kweights(st%nik))
  st%kpoints = M_ZERO
  
  ! just return the Gamma point
  if(st%nik == 1.or.(st%nik == 2.and.st%ispin == 2)) then 
    st%kweights(:) = M_ONE
    if(conf%periodic_dim>0) then 
      message(1) = 'Info: Only the Gamma k point is used'
      call write_info(1)
    end if
    return
  end if

  call loct_parse_logical('CenterOfInversion', .false., coi)
  if (conf%periodic_dim>1) then
    message(1) = 'Symmetries for periodic_dim > 1 not implemented yet'
    message(2) = 'K points are generated for the bull BZ'
    call write_warning(2)
    coi=.false.
  end if

  if(st%ispin == 2) then
    skip=M_TWO
  else
    skip=M_ONE
  end if

  ! if there is a center of inversion k is in [0.0,0.5]G, else in [0.0,1.0]G
  l = M_ZERO
  if (coi) then
    l = M_HALF 
  else
    l = M_ONE 
  end if
  
  do i=1,3
  if (st%nik_axis(i)==1) then
    tmp_nik_axis(i)=st%nik_axis(i)+1
  else
    tmp_nik_axis(i)=st%nik_axis(i)
  end if
  end do

  total_weight = M_ZERO
  k=M_ONE
  do kk = 1, st%nik_axis(3)
  do jk = 1, st%nik_axis(2)
  do ik = 1, st%nik_axis(1)
    st%kpoints(1, k) = l(1)*m%klat(1,1)*(ik-1)/(tmp_nik_axis(1)-1)
    st%kpoints(2, k) = l(2)*m%klat(2,2)*(jk-1)/(tmp_nik_axis(2)-1)
    st%kpoints(3, k) = l(3)*m%klat(3,3)*(kk-1)/(tmp_nik_axis(3)-1)
    if (coi) then
      if (conf%periodic_dim==M_ONE .and. (k==0 .or. k==st%nik_axis(i)-1)) then
        st%kweights(k) = M_ONE
      else
        st%kweights(k) = M_TWO
      end if
    else
      st%kweights(k) = M_ONE
    end if
    total_weight=total_weight+st%kweights(k)
    k = (k+1)/skip
  end do
  end do
  end do
  st%kweights=st%kweights/total_weight

end subroutine states_choose_kpoints

subroutine kpoints_write_info(st,iunit)
  
  type(states_type), intent(IN) :: st
  integer, intent(IN) :: iunit
  integer :: ik
    
  write(message(1),'(a,i4)') 'Number of k points in each direction = ',st%nik
  call write_info(1,iunit)
  
  write(message(1),'(3x,a,7x,a,9x,a,9x,a,8x,a)')'ik','K_x','K_y','K_z','Weight'
  message(2)='       --------------------------------------------------'
  call write_info(2,iunit,verbose_limit=80)
  do ik=1, st%nik
    write(message(1),'(i4,1x,4f12.4)') ik,st%kpoints(:,ik)*units_out%length%factor, st%kweights(ik)
    call write_info(1,iunit,verbose_limit=80)
  end do

end subroutine kpoints_write_info
