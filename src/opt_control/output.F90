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
  use lib_oct_parser_m
  use lib_oct_m
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
  use restart_m
  use messages_m

  implicit none

  private
  public :: oct_output,            &
            oct_prop_output_t,     &
            oct_prop_output_init,  &
            oct_prop_output_check, &
            oct_prop_output_end,   &
            oct_prop_output

  type oct_prop_output_t
    integer, pointer :: iter(:)
    integer :: niter
    character(len=100) :: dirname
  end type oct_prop_output_t
  
  integer, parameter :: NUMBER_CHECKPOINTS = 10

contains

  ! ---------------------------------------------------------
  subroutine oct_prop_output_init(prop_output, dirname, niter)
    type(oct_prop_output_t), intent(inout) :: prop_output
    character(len=*),        intent(in)    :: dirname
    integer,                 intent(in)    :: niter

    integer :: j

    prop_output%dirname = 'opt-control/'//trim(dirname)
    call io_mkdir(trim(prop_output%dirname))
    prop_output%niter = niter

    ALLOCATE(prop_output%iter(NUMBER_CHECKPOINTS), NUMBER_CHECKPOINTS)
    do j = 1, NUMBER_CHECKPOINTS
      prop_output%iter(j) = nint( real(niter)/(NUMBER_CHECKPOINTS+1) * j)
    end do

  end subroutine oct_prop_output_init


  ! ---------------------------------------------------------
  subroutine oct_prop_output_end(prop_output)
    type(oct_prop_output_t), intent(inout) :: prop_output

    integer :: j
    character(len=100) :: filename
    
    deallocate(prop_output%iter)
    ! This routine should maybe delete the files?

  end subroutine oct_prop_output_end


  ! ---------------------------------------------------------
  subroutine oct_prop_output_check(prop_output, psi, gr, geo, iter)
    type(oct_prop_output_t), intent(in)    :: prop_output
    type(states_t),          intent(in)    :: psi
    type(grid_t),            intent(in)    :: gr
    type(geometry_t),        intent(in)    :: geo
    integer,                 intent(in)    :: iter

    type(states_t) :: stored_st
    character(len=100) :: filename
    integer :: j, ierr
    FLOAT :: overlap

    do j = 1, NUMBER_CHECKPOINTS
     if(prop_output%iter(j) .eq. iter) then
       call states_copy(stored_st, psi)
       write(filename,'(a,i4.4)') trim(prop_output%dirname)//'/', j
       call restart_read(trim(filename), stored_st, gr, geo, ierr)
       overlap = abs( zstates_mpdotp(gr%m, stored_st, psi) )**2
       if( abs(overlap - M_ONE) > CNST(1.0e-10) ) then
          write(message(1), '(a,es13.4)') "WARNING: forward-backward propagation produced an error of", abs(overlap-M_ONE)
          call write_warning(1)
       end if
       call states_end(stored_st)
     end if
    end do

  end subroutine oct_prop_output_check



  ! ---------------------------------------------------------
  subroutine oct_prop_output(prop_output, iter, psi, gr)
    type(oct_prop_output_t), intent(inout) :: prop_output
    integer, intent(in) :: iter
    type(states_t), intent(inout) :: psi
    type(grid_t), intent(inout) :: gr

    integer :: j, ierr
    character(len=100) :: filename

    do j = 1, NUMBER_CHECKPOINTS
      if(prop_output%iter(j) .eq. iter) then
        write(filename,'(a,i4.4)') trim(prop_output%dirname)//'/', j
        call restart_write(trim(filename), psi, gr, ierr, iter)
      end if
    end do

  end subroutine oct_prop_output
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_output(iterator, gr, outp, final_st)
    type(oct_iterator_t), intent(inout) :: iterator
    type(grid_t), intent(inout)         :: gr
    type(output_t), intent(in)          :: outp
    type(states_t), intent(inout)       :: final_st

    character(len=80)  :: filename
    integer :: iunit, loop, ierr

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
