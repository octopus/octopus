!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

! ---------------------------------------------------------
subroutine states_choose_kpoints(dd, sb, geo)
  type(states_dim_t), intent(inout) :: dd
  type(simul_box_t),  intent(in)    :: sb
  type(geometry_t),   intent(in)    :: geo

  integer :: ik, iq

  call push_sub('states_kpoints_inc.states_choose_kpoints')

  dd%nik = kpoints_number(sb%kpoints)

  if (dd%ispin == SPIN_POLARIZED) dd%nik = 2*dd%nik
  
  SAFE_ALLOCATE(dd%kpoints (1:MAX_DIM, 1:dd%nik))
  SAFE_ALLOCATE(dd%kweights(1:dd%nik))

  do iq = 1, dd%nik
    ik = states_dim_get_kpoint_index(dd, iq)
    dd%kpoints(1:3, iq) = kpoints_get_point(sb%kpoints, ik)
    dd%kweights(iq) = kpoints_get_weight(sb%kpoints, ik)
  end do
  
  call print_kpoints_debug
  call pop_sub('states_kpoints_inc.states_choose_kpoints')

contains
  subroutine print_kpoints_debug
    integer :: iunit

    call push_sub('states_kpoints_inc.states_choose_kpoints.print_kpoints_debug')

    if(in_debug_mode) then
      
      iunit = io_open('debug/kpoints', action = 'write')
      call kpoints_write_info(dd, sb, iunit)      
      call io_close(iunit)

    end if

    call pop_sub('states_kpoints_inc.states_choose_kpoints.print_kpoints_debug')
  end subroutine print_kpoints_debug

end subroutine states_choose_kpoints


! ---------------------------------------------------------
subroutine kpoints_write_info(dd, sb, iunit)
  type(states_dim_t), intent(in) :: dd
  type(simul_box_t),  intent(in) :: sb
  integer,            intent(in) :: iunit

  integer :: ik, idir

  call push_sub('states_kpoints_inc.kpoints_write_info')

  if(sb%kpoints%method == KPOINTS_MONKH_PACK) then
    write(message(1),'(a,9(i3,1x))') 'Number of k-points in each direction = ', sb%kpoints%nik_axis(1:sb%dim)
    call write_info(1, iunit)
  else
    ! a Monkhorst-Pack grid was not used
    write(message(1),'(a,9(i3,1x))') 'Number of k-points = ', kpoint_index(dd, dd%nik)
    call write_info(1, iunit)
  endif

  write(message(1), '(3x,a,7x,a,9x,a,9x,a,8x,a)') 'ik', 'K_x', 'K_y', 'K_z', 'Weight'
  message(2) = '   --------------------------------------------------'
  call write_info(2, iunit, verbose_limit=80)

  do ik = 1, dd%nik
     write(message(1),'(i4,1x,4f12.4)') &
       ik, (units_from_atomic(unit_one/units_out%length, dd%kpoints(idir, ik)), idir=1, sb%dim), dd%kweights(ik)
     call write_info(1, iunit, verbose_limit=80)
  end do

  call pop_sub('states_kpoints_inc.kpoints_write_info')
end subroutine kpoints_write_info


! ---------------------------------------------------------
logical pure function kpoint_is_gamma(this, ik)
  type(states_dim_t), intent(in) :: this
  integer,            intent(in) :: ik
  
  kpoint_is_gamma = (maxval(abs(this%kpoints(:, ik))) < M_EPSILON)

end function kpoint_is_gamma


! ---------------------------------------------------------
integer function kpoint_index(this, ik) result(index)
  type(states_dim_t), intent(in) :: this
  integer,            intent(in) :: ik

  call push_sub('states_kpoints_inc.kpoint_index')

  if (this%nspin == 2) then
    if (states_dim_get_spin_index(this, ik) == 1) then
      index = (ik + 1) / 2
    else
      index = ik / 2
    endif
  else
     index = ik
  endif

  call pop_sub('states_kpoints_inc.kpoint_index')
end function kpoint_index

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
