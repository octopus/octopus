!! Copyright (C) 2005-2006 M. Marques, X. Andrade
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$


subroutine X(restart_write_lr_density)(sys, lr, omega, tag)
  type(system_t),  intent(inout) :: sys
  type(lr_t),         intent(in) :: lr
  FLOAT,              intent(in) :: omega
  integer,            intent(in) :: tag

  character(len=80) :: fname
  integer :: is, ierr

  call block_signals()
  do is = 1, sys%st%d%nspin
    write(fname, '(a, f6.4, a, i1, a, i1)') 'density-', omega, '-', tag, '-', is
    call X(restart_write_function)(trim(tmpdir)//RESTART_DIR, fname, sys%gr,&
         lr%X(dl_rho)(:, is), ierr, size(lr%X(dl_rho),1))
  end do
  call unblock_signals()

end subroutine X(restart_write_lr_density)


subroutine X(restart_read_lr_density)(sys, lr, omega, tag, ierr)
  type(system_t),  intent(inout) :: sys
  type(lr_t),      intent(inout) :: lr
  FLOAT,              intent(in) :: omega
  integer,            intent(in) :: tag
  integer,         intent(inout) :: ierr

  character(len=80) :: fname
  integer :: is, s_ierr
  FLOAT :: closest_omega

  ierr = 0;
  do is = 1, sys%st%d%nspin
    write(fname, '(a, f6.4, a, i1, a, i1)') 'density-', omega, '-', tag, '-', is
    call X(restart_read_function)(trim(tmpdir)//RESTART_DIR, fname, sys%gr%m,&
         lr%X(dl_rho)(:, is), s_ierr)
    if( s_ierr /=0 ) ierr = s_ierr;
  end do


  if( ierr == 0 ) then 
    write(message(1),'(a, f6.4)') 'Loaded restart density for frequency ', omega/units_out%energy%factor
    call write_info(1)

  else

    write(message(1),'(a, f6.4)') 'Could not load restart density for frequency ', omega/units_out%energy%factor
    call write_info(1)
    
    !search for the density of the closest frequency
    closest_omega = omega
    call oct_search_file_lr(closest_omega, tag, ierr, trim(tmpdir)//RESTART_DIR)
    
    !atempt to read 
    if(ierr == 0 ) then 
      
      do is = 1, sys%st%d%nspin
        write(fname, '(a, f6.4, a, i1, a, i1)') 'density-', closest_omega, '-', tag, '-', is
        call X(restart_read_function)(trim(tmpdir)//RESTART_DIR, fname, sys%gr%m,&
             lr%X(dl_rho)(:, is), s_ierr)
        if( s_ierr /=0 ) ierr = s_ierr;
      end do
      
    end if

    if(ierr == 0 ) then 

      write(message(1),'(a, f6.4)') 'Using restart density from frequency ', closest_omega/units_out%energy%factor
      call write_info(1)

    end if

  endif

end subroutine X(restart_read_lr_density)


