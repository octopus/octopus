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

#include "global.h"

module opt_control_output_m
  use datasets_m
  use varinfo_m
  use global_m
  use loct_parser_m
  use loct_m
  use io_m
  use opt_control_constants_m
  use opt_control_parameters_m
  use opt_control_iter_m
  use grid_m
  use states_m
  use geometry_m
  use excited_states_m
  use states_output_m
  use output_m
  use messages_m

  implicit none

  private
  public :: oct_output


contains

  ! ---------------------------------------------------------
  subroutine oct_output(iterator, gr, outp, final_st)
    type(oct_iterator_t), intent(inout) :: iterator
    type(grid_t), intent(inout)         :: gr
    type(output_t), intent(in)          :: outp
    type(states_t), intent(inout)       :: final_st

    character(len=80)  :: filename
    integer :: iunit, loop

    call push_sub('output.oct_output')
    
    iunit = io_open('opt-control/info', action='write')
    write(iunit, '(a,i4)')     'Total Iterations = ', iterator%ctr_iter
    write(iunit, '(a)') 
    write(iunit, '(a)')        'Best value of target functional'
    write(iunit, '(a,i4)')     'Iteration  = ', iterator%bestJ1_ctr_iter
    write(iunit, '(a,f14.8)')  'Overlap    = ', iterator%bestJ1
    write(iunit, '(a,es18.8)') 'Functional = ', iterator%bestJ1_J
    write(iunit, '(a,f14.8)')  'Fluence = ',    iterator%bestJ1_fluence
    call io_close(iunit)
    message(1) = "Info: Output States"
    call write_info(1)
    ! should output wavefunctions ;)
    
    ! dump convergence: J,P,fluence,penalty
    iunit = io_open('opt-control/convergence', action='write')
    ! header
    write(iunit, '(4(a))') '# iteration', '  J[Psi,chi,epsilon]', &
                                          '            J_1[Psi]', &
                                          '        J_2[epsilon]'
    write(iunit, '(a)') &
     '#######################################################################'

    ! data
    do loop = 1, iterator%ctr_iter_max
       write(iunit, '(i11,3f20.8)') loop, iterator%convergence(1,loop), &
         iterator%convergence(2,loop), iterator%convergence(3,loop)
    end do
    call io_close(iunit)

    write(filename,'(a)') 'opt-control/laser.bestJ1'
    call parameters_write(filename, iterator%best_par, fourier = .true.)

    call states_output(final_st, gr, 'opt-control/final', outp)

    call pop_sub()
  end subroutine oct_output

end module opt_control_output_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
