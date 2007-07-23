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
!! $Id: td_transport.F90 3090 2007-07-20 09:39:28Z lorenzen $

! Entry point of the transport runmode. Work is delegated to the
! timedep_m module after a few initalization steps.

#include "global.h"

module td_transport_m
  use datasets_m
  use global_m
  use hamiltonian_m
  use io_m
  use messages_m
  use lib_oct_parser_m
  use profiling_m
  use simul_box_m
  use states_m
  use system_m
  use td_trans_rti_m
  use td_rti_m
  use timedep_m
  use units_m
  use varinfo_m

  implicit none

  private
  public ::           &
    td_transport_run

contains

  ! ---------------------------------------------------------
  ! Read and check input, allocate memory, and initialize
  ! substructures.
  subroutine td_transport_init(sys)
    type(system_t), intent(in) :: sys

    integer :: inp_int

    call push_sub('td_transport.td_transport_init')

    ! Restrictions.
    ! Some checks that the input file does not request too much unimplemented
    ! features:
    !   *) only non-interacting electrons work for now,
    !   *) the only simulation box permitted is the parallel epiped,
    !   *) the time evolution method cannot be chosen freely,
    !   *) only spin unpolarized calculations possible.
    if(.not.sys%ks%ip_app) then
      message(1) = 'Transport calculations for interacting electrons are'
      message(2) = 'not yet possible. Please include'
      message(3) = ''
      message(4) = '  NonInteractingElectrons = yes'
      message(5) = 'in your input file.'
      call write_fatal(5)
    end if
    
    if(sys%gr%sb%box_shape.ne.PARALLELEPIPED) then
      message(1) = 'Transport calculations are only possible with'
      message(2) = ''
      message(3) = '  BoxShape = parallelepided'
      message(4) = ''
      message(5) = 'in the input file.'
      call write_fatal(5)
    end if

    ! Check the propagator. If none is set, set it to PROP_CRANK_NICHOLSON_SRC_MEM,
    ! if another propagator is chosen, give an error message.
    call loct_parse_int(check_inp('TDEvolutionMethod'), -1, inp_int)
    if(varinfo_valid_option('TDEvolutionMethod', inp_int).and. &
      inp_int.ne.PROP_CRANK_NICHOLSON_SRC_MEM) then
      message(1) = 'The time evolution method for time dependent cannot'
      message(2) = 'be chosen freely. The Crank-Nicholson propagator'
      message(3) = 'with source and memory term has to be used. Either set'
      message(4) = ''
      message(5) = '  TDEvolutionmethod = crank_nicholson_src_mem'
      message(6) = ''
      message(7) = 'in your input or remove the TDEvolutionMethod variable.'
      call write_fatal(7)
    end if
    if(.not.varinfo_valid_option('TDEvolutionMethod', inp_int)) then
      call loct_parse_putsym('TDEvolutionMethod', PROP_CRANK_NICHOLSON_SRC_MEM)
    end if

    if(sys%st%d%ispin.ne.UNPOLARIZED) then
      message(1) = 'Only spin unpolarized transport calculations are possible.'
      message(2) = 'Set'
      message(3) = ''
      message(4) = '  SpinComponents = unpolarized'
      message(5) = ''
      message(6) = 'in your input file or remove the SpinComponents line entirely.'
      call write_fatal(6)
    end if

    call pop_sub()
  end subroutine td_transport_init


  ! ---------------------------------------------------------
  ! Write some status information to stdout.
!   subroutine td_transport_write_info(trans, sys, unit)
!     type(transport_t), intent(in) :: trans
!     type(system_t),    intent(in) :: sys
!     integer,           intent(in) :: unit

!     call push_sub('td_transport.td_transport_write_info')

!     call messages_print_stress(stdout, 'Transport')

!     write(message(1), '(a, i10, i10)') 'Points in interface regions:     ', &
!       trans%intface%np, trans%intface%np
!     write(message(2), '(a, i10)') 'MBytes required for memory term: ', &
!       mbytes_memory_term(trans%max_iter, trans%intface%np, NLEADS, sys)

!     call write_info(2, unit)

!     call messages_print_stress(stdout)

!     call pop_sub()
!   end subroutine td_transport_write_info


  ! ---------------------------------------------------------
  ! Initialize, run, and finalize the time propagation.
  subroutine td_transport_run(sys, h, from_scratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    logical,             intent(inout) :: from_scratch

    ! To avoid unused variable warnings. Restart will be implemented
    ! later.
    if(from_scratch) then
    end if

    call push_sub('td_transport.td_transport_run')

    call td_transport_init(sys)

    call td_run(sys, h, from_scratch)

    call td_transport_end()

    call pop_sub()
  end subroutine td_transport_run


  ! ---------------------------------------------------------
  ! Deallocate all memory.
  subroutine td_transport_end()
    call push_sub('td_transport.td_transport_end')

    call pop_sub()
  end subroutine td_transport_end
end module td_transport_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
