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
!!
!! $Id$

subroutine states_choose_kpoints(d, m, sb, geo)
  type(states_dim_type), intent(inout) :: d
  type(mesh_type),       intent(in)    :: m
  type(simul_box_type),  intent(in)    :: sb
  type(geometry_type),   intent(in)    :: geo

  integer  :: coi, i, nkmax
  integer(POINTER_SIZE) :: blk
  FLOAT :: total_weight, kmax
  
  ! local variables for the crystal_init call
  integer :: is, nk
  integer, allocatable :: natom(:)      ! natom(i) is the number of atoms of specie i
  FLOAT :: kshifts(3)   
  FLOAT, allocatable :: coorat(:,:,:)   ! coorat(i,j,k) is the k-th component (lattice coordinates) 
                                        ! of the position of the j-th atom of type i.
  FLOAT, allocatable :: kp(:,:),kw(:)

  call push_sub('states_choose_kpoints')
  
  ! if not periodic just return the Gamma point
  if (sb%periodic_dim == 0) then
    select case(d%ispin)
    case(1,3)
      d%nik = 1
    case(2)
      d%nik = 2
    case default    
      message(1) = 'Input: invalid SpinComponents'
      call write_fatal(1)
    end select

    allocate(d%kpoints(3,d%nik),d%kweights(d%nik))
    d%kpoints  = M_ZERO
    d%kweights = M_ONE

    call pop_sub()
    return
  end if

  if(loct_parse_block(check_inp('NumberKPoints'), blk) .ne. 0) then
    message(1) = 'Block "NumberKPoints" not found in input file.'
    call write_fatal(1)
  end if

  d%nik_axis = 1
  do i = 1, sb%periodic_dim
    call loct_parse_block_int(blk, 0, i-1, d%nik_axis(i))
  end do
  call loct_parse_block_end(blk)
  if (any(d%nik_axis < 1)) then
    message(1) = 'Input: NumberKPoints is not valid'
    message(2) = '(NumberKPoints >= 1)'
    call write_fatal(2)
  end if
  nkmax = PRODUCT(d%nik_axis)

  if(loct_parse_block(check_inp('ShiftKPoints'), blk) .ne. 0) then
    kshifts = M_ZERO
  else
    do i = 1, sb%periodic_dim
      call loct_parse_block_float(blk, 0, i-1, kshifts(i))
    end do
    call loct_parse_block_end(blk)
  end if

  if(sb%periodic_dim == 1) then
    call loct_parse_int(check_inp('CenterOfInversion'), 0, coi)
    d%nik = d%nik_axis(1)/(1 + coi) + 1
    allocate(d%kpoints(3, d%nik), d%kweights(d%nik))
    d%kpoints     = M_ZERO
    d%kweights    = M_ONE
    kmax          = (M_ONE - coi*M_HALF)*sb%klat(1,1)
    total_weight  = M_ZERO

    do i = 1, d%nik
      d%kpoints(1, i) = (i-1)*kmax/real(d%nik-1)
      if (i /= 1 .and. i /= d%nik) then
        d%kweights(i) = d%kweights(i) + coi
      end if
      total_weight = total_weight + d%kweights(i)
    end do
    d%kweights = d%kweights/total_weight
    if (d%ispin == 2) then
      message(1) = 'Not implemented yet.'
      call write_fatal(1)
    end if
    return
  end if
    
  allocate(natom(geo%nspecies), coorat(geo%nspecies, geo%natoms, 3))
  
  natom  = 0
  coorat = M_ZERO
  do i = 1, geo%natoms
    is = geo%atom(i)%spec%index
    natom(is) = natom(is) + 1
    coorat(is, natom(is), :) = geo%atom(i)%x(:)
  end do
 
  do i = 1, 3
    coorat(:,:,i) = coorat(:,:,i) / sb%rlat(i, i) + M_HALF
  end do


  allocate(kp(3, nkmax), kw(nkmax))

  call crystal_init(sb%rlat, geo%nspecies, natom, geo%natoms, coorat, d%nik_axis, &
       kshifts, nk, kp, kw)

  ! double d%nik and copy points for spin polarized calc
  select case(d%ispin)
  case(1)
    d%nik = nk
    allocate(d%kpoints(3, d%nik), d%kweights(d%nik))  
    do i = 1, 3
      d%kpoints(i,:) = kp(i,:)*sb%klat(i,i)
    end do
    d%kweights = kw
  case(2)
    d%nik = 2 * nk
    allocate(d%kpoints(3, d%nik), d%kweights(d%nik))  
    do i = 1,3
      d%kpoints(i,::2)  = kp(i,:)*sb%klat(i,i)
      d%kpoints(i,2::2) = kp(i,:)*sb%klat(i,i)
    end do
    d%kweights(::2)  = kw(:)
    d%kweights(2::2) = kw(:)
  end select
  
  deallocate(natom, coorat)
  deallocate(kp, kw)

  call pop_sub()
end subroutine states_choose_kpoints

subroutine kpoints_write_info(d,iunit)
  type(states_dim_type), intent(IN) :: d
  integer, intent(IN) :: iunit

  integer :: ik
    
  call push_sub('kpoints_write_info')
  
  write(message(1),'(a,3(i3,1x))') 'Number of k points in each direction = ', d%nik_axis(:)
  call write_info(1, iunit)
  
  write(message(1), '(3x,a,7x,a,9x,a,9x,a,8x,a)') 'ik', 'K_x', 'K_y', 'K_z', 'Weight'
  message(2) = '       --------------------------------------------------'
  call write_info(2, iunit, verbose_limit=80)

  do ik = 1, d%nik
    write(message(1),'(i4,1x,4f12.4)') ik,d%kpoints(:,ik)*units_out%length%factor, d%kweights(ik)
    call write_info(1,iunit,verbose_limit=80)
  end do

  call pop_sub()
end subroutine kpoints_write_info
