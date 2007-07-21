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

! Main module of the time dependent electronic transport implementation.

#include "global.h"

module td_transport_m
  use datasets_m
  use global_m
  use hamiltonian_m
  use io_m
  use messages_m
  use lib_oct_parser_m
  use simul_box_m
  use states_m
  use system_m
  use td_trans_mem_m
  use td_trans_intf_m
  use td_trans_rti_m
  use units_m
  use varinfo_m

  implicit none

  private
  public ::           &
    td_transport_run, &
    transport_t

  type transport_t
    FLOAT            :: dt
    FLOAT            :: delta    ! dt/2 (for convenience).
    integer          :: max_iter

    CMPLX, pointer   :: mem_coeff(:, :, :, :) ! i, j, t, il

    type(intface_t)  :: intface

    CMPLX, pointer   :: st_intface(:, :, :, :)
  end type transport_t

contains

  ! ---------------------------------------------------------
  ! Read and check input, allocate memory, and initialize
  ! substructures.
  subroutine td_transport_init(trans, h, sys)
    type(transport_t),   intent(out) :: trans
    type(hamiltonian_t), intent(in)  :: h
    type(system_t),      intent(in)  :: sys

    integer :: inp_int, allocsize
    logical :: file_exists

    call push_sub('td_transport.td_transport_init')

    ! Get the TD parameters, i. e. size of a timestep and the iteration
    ! number. Actually, this is copied from td_init. The problem here is
    ! that the td_transport run mode is actually a special case of the
    ! td runmode with a certain propagator, geometry, and special boundary
    ! conditions. For convenience during development, it is casted into its
    ! own runmode, perhaps later, the possibility arises to unify it with
    ! the td runmode. The same holds for the (deferred to a much later time)
    ! gs_transport runmode that calculates a Kohn-Sham groundstate for the
    ! transport geometry.
    ! For now, we stick to non-interacting electrons...
    call loct_parse_float(check_inp('TDTimeStep'), CNST(0.07)/units_inp%time%factor, trans%dt)
    trans%dt = trans%dt * units_inp%time%factor
    if(trans%dt <= M_ZERO) then
      write(message(1),'(a,f14.6,a)') "Input: '", trans%dt, "' is not a valid TDTimeStep"
      message(2) = '(0 < TDTimeStep)'
      call write_fatal(2)
    end if
    trans%delta = trans%dt/2

    call loct_parse_int(check_inp('TDMaximumIter'), 1500, trans%max_iter)
    if(trans%max_iter < 1) then
      write(message(1), '(a,i6,a)') "Input: '", trans%max_iter, "' is not a valid TDMaximumIter"
      message(2) = '(1 <= TDMaximumIter)'
      call write_fatal(2)
    end if

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
    
    call loct_parse_int(check_inp('TDEvolutionMethod'), -1, inp_int)
    if(varinfo_valid_option('TDEvolutionMethod', inp_int)) then
      message(1) = 'The time evolution method for time dependent cannot'
      message(2) = 'be chosen freely. The hard-wired Crank-Nicholson propagator'
      message(3) = 'with source and memory term is used in any case.'
      message(4) = 'Remove definition of the TDEvolutionMethod variable from'
      message(5) = 'your input file.'
      call write_fatal(5)
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

    ! Create output directory.
    inquire(file='transport', exist=file_exists)
    if(.not.file_exists) then
      call io_mkdir('transport')
    end if

    ! Initialization steps.
    call intface_init(sys%gr, trans%intface)
    call memory_init(trans%intface, trans%delta, trans%max_iter, &
      sys%gr%f_der%der_discr%lapl, trans%mem_coeff)
    allocsize = trans%intface%np*(sys%st%st_start-sys%st%st_end)*trans%max_iter*NLEADS
    ALLOCATE(trans%st_intface(trans%intface%np, sys%st%st_end-sys%st%st_start, trans%max_iter, NLEADS), allocsize)
#ifdef HAVE_SPARSKIT
    call cn_src_mem_init(sys%gr)
#endif

    call pop_sub()
  end subroutine td_transport_init


  ! ---------------------------------------------------------
  ! Write some status information to stdout.
  subroutine td_transport_write_info(trans, sys, unit)
    type(transport_t), intent(in) :: trans
    type(system_t),    intent(in) :: sys
    integer,           intent(in) :: unit

    call push_sub('td_transport.td_transport_write_info')

    call messages_print_stress(stdout, 'Transport')

    write(message(1), '(a, i10, i10)') 'Points in interface regions:     ', &
      trans%intface%np, trans%intface%np
    write(message(2), '(a, i10)') 'MBytes required for memory term: ', &
      mbytes_memory_term(trans%max_iter, trans%intface%np, NLEADS, sys)

    call write_info(2, unit)

    call messages_print_stress(stdout)

    call pop_sub()
  end subroutine td_transport_write_info


  ! ---------------------------------------------------------
  ! Initialize, run, and finalize the time propagation.
  subroutine td_transport_run(sys, h, from_scratch)
    type(system_t),      intent(in) :: sys
    type(hamiltonian_t), intent(in) :: h
    logical,             intent(in) :: from_scratch

    type(transport_t) :: trans

    ! To avoid unused variable warnings. Restart will be implemented
    ! later.
    if(from_scratch) then
    end if

    call push_sub('td_transport.td_transport_run')

    call td_transport_init(trans, h, sys)

    call td_transport_write_info(trans, sys, stdout)

    ! Do the propagation.


!     message(1) = 'Time dependent quantum transport not yet implemented.'
!     call write_fatal(1)

    call td_transport_end(trans)

    call pop_sub()
  end subroutine td_transport_run


  ! ---------------------------------------------------------
  ! Deallocate all memory.
  subroutine td_transport_end(trans)
    type(transport_t), intent(inout) :: trans

    call push_sub('td_transport.td_transport_end')

    call intface_end(trans%intface)
    call memory_end(trans%mem_coeff)
#ifdef HAVE_SPARSKIT
    call cn_src_mem_end()
#endif
    if(associated(trans%st_intface)) then
      deallocate(trans%st_intface)
      nullify(trans%st_intface)
    end if

    call pop_sub()
  end subroutine td_transport_end
end module td_transport_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
