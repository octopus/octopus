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

  integer  :: coi, i, ik, jk, kk, j, k, nkmax, skip
  FLOAT :: total_weight, kmax
  
! local variables for the crystal_init call
  character(len=80) :: str, label
  integer :: nspecies, natoms
  integer :: is, nk
  integer, allocatable :: natom(:)
  FLOAT :: kshifts(3)   
  FLOAT, allocatable :: coorat(:,:,:)   
  FLOAT, allocatable :: kp(:,:),kw(:)
  character(len=80),allocatable :: splabel(:)    
  
! nspecies is the total number of species
! natoms is the total number of atoms
! natom(i) is the number of atoms of specie i
! coorat(i,j,k) is the k-th component (lattice coordinates) 
!     of the position of the j-th atom of type i.
! splabel(i) is the label of each specie

  ! if not periodic just return the Gamma point
  if (conf%periodic_dim == 0) then
    select case(st%d%ispin)
    case(1)
      st%d%nik = 1
      allocate(st%d%kpoints(3,st%d%nik),st%d%kweights(st%d%nik))
    case(2)
      st%d%nik = 2
      allocate(st%d%kpoints(3,st%d%nik),st%d%kweights(st%d%nik))
    case default    
      message(1) = 'Input: invalid SpinComponents'
      call write_fatal(1)
    end select
    st%d%kpoints = M_ZERO
    st%d%kweights = M_ONE
    return
  else
    if(loct_parse_block_n('NumberKPoints')<1) then
      message(1) = 'Block "NumberKPoints" not found in input file.'
      call write_fatal(1)
    end if
    st%d%nik_axis = 1
    do i = 1, conf%periodic_dim
      call loct_parse_block_int('NumberKPoints', 0, i-1, st%d%nik_axis(i))
    end do
    if (any(st%d%nik_axis < 1)) then
      message(1) = 'Input: NumberKPoints is not valid'
      message(2) = '(NumberKPoints >= 1)'
      call write_fatal(2)
    end if
    nkmax = PRODUCT(st%d%nik_axis)
    if(loct_parse_block_n('ShiftKPoints')<1) then
      kshifts = M_ZERO
    else
      do i = 1, conf%periodic_dim
        call loct_parse_block_float('ShiftKPoints', 0, i-1, kshifts(i))
      end do
    end if
    if (conf%periodic_dim == 1) then
      call loct_parse_int('CenterOfInversion', 0, coi)
      st%d%nik = st%d%nik_axis(1)/(1 + coi)+1
      allocate(st%d%kpoints(3,st%d%nik),st%d%kweights(st%d%nik))
      st%d%kpoints = M_ZERO
      st%d%kweights = M_ONE
      kmax = (M_ONE - coi*M_HALF)*m%klat(1,1)
      total_weight=M_ZERO
      do i = 1, st%d%nik
        st%d%kpoints(1,i) = (i-1)*kmax/real(st%d%nik-1)
        if (i /= 1 .and. i /= st%d%nik) then
          st%d%kweights(i) = st%d%kweights(i)+coi
        end if
        total_weight=total_weight+st%d%kweights(i)
      end do
      st%d%kweights=st%d%kweights/total_weight
      if (st%d%ispin == 2) then
        message(1) = 'Not implemented yet.'
        call write_fatal(1)
      end if
      return
    end if
  end if
    
! Read the atomic species and coordinates

  str = "Species"
  nspecies = loct_parse_block_n(str)
  if (nspecies < 1) then
    message(1) = "Input: Species block not specified"
    message(2) = '% Species'
    message(3) = '   specie <params>'
    message(4) = '%'
    call write_fatal(4)    
  end if

  allocate(natom(nspecies),splabel(nspecies))
  
  do i = 1, nspecies
    call loct_parse_block_string(str, i-1, 0, splabel(i))
  enddo  

  ! read the positions
  str = "Coordinates"
  natoms = loct_parse_block_n(str)
  if(natoms <= 0) then
    message(1) = "Input: Coordinates block not specified"
    message(2) = '% Coordinates'
    message(3) = '  specie  x  y  z  move'
    message(4) = '%'
    call write_fatal(4)
  end if

  allocate(coorat(nspecies,natoms,3))
  
  natom = 0
  do i = 1, natoms
    call loct_parse_block_string(str, i-1, 0, label)
      is = get_specie_number(label,splabel,nspecies,natom)
    call loct_parse_block_float (str, i-1, 1, coorat(is,natom(is),1))
    call loct_parse_block_float (str, i-1, 2, coorat(is,natom(is),2))
    call loct_parse_block_float (str, i-1, 3, coorat(is,natom(is),3))
  end do
 
  do i=1,3
    coorat(:,:,i) = coorat(:,:,i) / m%rlat(i,i) + M_HALF
  end do

  allocate(kp(3,nkmax),kw(nkmax))

  call init_crystal(m%rlat,nspecies,natom,natoms,coorat,st%d%nik_axis, &
                    kshifts,nk,kp,kw)


  ! double st%d%nik and copy points for spin polarized calc
  select case(st%d%ispin)
  case(1)
    st%d%nik = nk
    allocate(st%d%kpoints(3, st%d%nik), st%d%kweights(st%d%nik))  
    do i = 1,3
      st%d%kpoints(i,:) = kp(i,:)*m%klat(i,i)
    end do
    st%d%kweights = kw
  case(2)
    st%d%nik = 2 * nk
    allocate(st%d%kpoints(3, st%d%nik), st%d%kweights(st%d%nik))  
    do i = 1,3
      st%d%kpoints(i,::2) = kp(i,:)*m%klat(i,i)
      st%d%kpoints(i,2::2) = kp(i,:)*m%klat(i,i)
    end do
    st%d%kweights(::2) = kw(:)
    st%d%kweights(2::2) = kw(:)
  end select
  
  deallocate(natom,splabel,coorat)
  deallocate(kp,kw)

contains

  integer function get_specie_number(label,splabel,nspecies,natom)
    
    character(len=*),intent(in) :: label
    character(len=*), intent(in) :: splabel(:)
    integer, intent(in) :: nspecies
    integer, intent(inout) :: natom(:)

    integer :: j
    logical :: flag

    flag = .false.
    do j = 1, nspecies
      if( trim(label) == trim(splabel(j)) ) then
        flag = .true.
        get_specie_number = j
        natom(j) = natom(j) + 1
        return
      end if
    end do
    
    if(.not. flag) then
      message(1) = "Specie '"//trim(label)//"' not found"
      call write_fatal(1)
    end if
    
  end function get_specie_number
  
end subroutine states_choose_kpoints

subroutine kpoints_write_info(st,iunit)
  
  type(states_type), intent(IN) :: st
  integer, intent(IN) :: iunit
  integer :: ik
    
  write(message(1),'(a,3(i3,1x))') 'Number of k points in each direction = ',st%d%nik_axis(:)
  call write_info(1,iunit)
  
  write(message(1),'(3x,a,7x,a,9x,a,9x,a,8x,a)')'ik','K_x','K_y','K_z','Weight'
  message(2)='       --------------------------------------------------'
  call write_info(2,iunit,verbose_limit=80)
  do ik=1, st%d%nik
    write(message(1),'(i4,1x,4f12.4)') ik,st%d%kpoints(:,ik)*units_out%length%factor, st%d%kweights(ik)
    call write_info(1,iunit,verbose_limit=80)
  end do

end subroutine kpoints_write_info
