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
  subroutine oct_finalcheck(oct, initial_st, target_st, sys, h, td, tdt)
    type(oct_t), intent(in)            :: oct
    type(states_t), intent(in)         :: initial_st
    type(states_t), intent(in)         :: target_st
    type(system_t), intent(inout)      :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(td_t), intent(inout)          :: td
    type(td_target_set_t), intent(inout) :: tdt


    type(states_t) :: psi
    FLOAT, allocatable :: td_fitness(:)
    FLOAT :: overlap

    call push_sub('opt_control_finalcheck.oct_finalcheck')

    if(.not.oct%oct_double_check) then
      call pop_sub(); return
    end if

    message(1) = "Info: Optimization finished...checking the field"
    call write_info(1)
    
    ALLOCATE(td_fitness(0:td%max_iter), td%max_iter+1)
    psi = initial_st
    
    call propagate_forward(oct, sys, h, td, tdt, psi, write_iter = .true.)

    if(oct%targetmode==oct_targetmode_td) then
      overlap = SUM(td_fitness) * abs(td%dt)
      write(message(1), '(a,f14.8)') " => overlap:", overlap
      call write_info(1)
    else
      overlap = abs(zstates_mpdotp(sys%gr%m, psi, target_st))
      write(message(1), '(6x,a,f14.8)') " => overlap:", overlap
      call write_info(1)
    end if

    call states_end(psi)
    deallocate(td_fitness)
    call pop_sub()
  end subroutine oct_finalcheck

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
