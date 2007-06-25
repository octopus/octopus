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
!! $Id: opt_control.F90 2875 2007-04-30 16:54:15Z acastro $


  ! ---------------------------------------------------------
  subroutine output(oct, iterator)
    type(oct_t), intent(in)          :: oct
    type(oct_iterator_t), intent(in) :: iterator

    integer :: iunit, loop, ierr

    call push_sub('opt_control.output')
    
    iunit = io_open('opt-control/info', action='write')
    write(iunit, '(a,i4)')    'Total Iterations = ', iterator%ctr_iter
    write(iunit, '(a,f14.8)') 'Last Overlap    = ', iterator%overlap
    write(iunit, '(a,f14.8)') 'Last Functional = ', iterator%functional
    write(iunit, '(a)') 
    write(iunit, '(a)')       'Best value of functional'
    write(iunit, '(a,i4)')    'Iteration  = ', iterator%bestJ_ctr_iter
    write(iunit, '(a,f14.8)') 'Overlap    = ', iterator%bestJ_J1
    write(iunit, '(a,f14.8)') 'Functional = ', iterator%bestJ
    write(iunit, '(a,f14.8)') 'Fluence = ',    iterator%bestJ_fluence
    write(iunit, '(a)') 
    write(iunit, '(a)')       'Best value of target functional'
    write(iunit, '(a,i4)')    'Iteration  = ', iterator%bestJ1_ctr_iter
    write(iunit, '(a,f14.8)') 'Overlap    = ', iterator%bestJ1
    write(iunit, '(a,es18.8)') 'Functional = ', iterator%bestJ1_J
    write(iunit, '(a,f14.8)') 'Fluence = ',    iterator%bestJ1_fluence
    call io_close(iunit)
    message(1) = "Info: Output States"
    call write_info(1)
    ! should output wavefunctions ;)
    
    ! dump convergence: J,P,fluence,penalty
    iunit = io_open('opt-control/convergence', action='write')
    ! header
    write(iunit, '(4(a))') '# iteration ','functional ','overlap ','penalty '
    ! data
    do loop = 1, iterator%ctr_iter_max
       write(iunit, '(i6,3f18.8,es20.10)') loop, iterator%convergence(1,loop), iterator%convergence(2,loop), &
         iterator%convergence(3,loop), iterator%convergence(4,loop)
    end do
    call io_close(iunit)

    call pop_sub()
  end subroutine output



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
