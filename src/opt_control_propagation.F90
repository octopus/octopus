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
!! $Id: opt_control.F90 2873 2007-04-29 22:05:29Z acastro $

#include "global.h"


  ! ---------------------------------------------------------
  subroutine propagate_forward(oct, sys, h, td, laser, v_old_i, v_old_f, &
                               td_tg, tdtarget, td_fitness, psi_n, write_iter)
    type(oct_t),         intent(in)    :: oct
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(td_t),          intent(inout) :: td
    FLOAT, pointer                     :: laser(:, :)
    FLOAT, pointer                     :: v_old_i(:, :, :), v_old_f(:, :, :)
    type(td_target_t), pointer         :: td_tg(:)
    CMPLX, intent(in)                  :: tdtarget(:)
    FLOAT, intent(out)                 :: td_fitness(:)
    type(states_t),      intent(inout) :: psi_n
    logical, optional, intent(in)      :: write_iter

    integer :: ierr, ii, i
    logical :: write_iter_ = .false.
    FLOAT   :: etime
    FLOAT, allocatable :: dens(:,:)
    type(grid_t),  pointer :: gr
    type(td_write_t)           :: write_handler

    call push_sub('opt_control.propagate_forward')

    message(1) = "Info: Forward propagating Psi"
    call write_info(1)

    write_iter_ = .false.
    if(present(write_iter)) write_iter_ = write_iter

    gr => sys%gr

    if(write_iter_) then
      call td_write_init(write_handler, gr, sys%st, sys%geo, (td%move_ions>0), h%ep%with_gauge_field, td%iter, td%dt)
    end if

    ALLOCATE(dens(NP_PART, psi_n%d%nspin), NP_PART*psi_n%d%nspin) 

    ! setup the hamiltonian
    call states_calc_dens(psi_n, NP_PART, dens)
    psi_n%rho = dens
    call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
    ! setup start of the propagation
    do i = 1, psi_n%d%nspin
      v_old_i(:, i, 2) = h%vhxc(:, i)
    end do
    v_old_f(:, :, 3) = v_old_i(:, :, 2)

    h%ep%lasers(1)%numerical => laser
    td%tr%v_old => v_old_i

    h%ep%lasers(1)%dt = td%dt      
    !call loct_progress_bar(-1, td%max_iter-1)

    ii = 1

    etime = loct_clock()
    do i = 1, td%max_iter 
      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), abs(td%dt))

      ! if td_target
      if(oct%targetmode==oct_targetmode_td) &
        call calc_tdfitness(td_tg, gr, psi_n, tdtarget, td_fitness(i))
      
      ! update
      call states_calc_dens(psi_n, NP_PART, dens)
      psi_n%rho = dens
      call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
      call hamiltonian_energy(h, sys%gr, sys%geo, psi_n, -1)

      ! only write in final run
      if(write_iter_) then

        write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
          i*td%dt       / units_out%time%factor, &
          (h%etot + sys%geo%kinetic_energy) / units_out%energy%factor, &
          loct_clock() - etime
        call write_info(1)
        etime = loct_clock()

        ii = ii + 1 
        if(ii==sys%outp%iter+1 .or. i == td%max_iter) then ! output 
          if(i == td%max_iter) sys%outp%iter = ii - 1 
          ii = 1 
          call td_write_iter(write_handler, gr, psi_n, h, sys%geo, td%kick, td%dt, i)
          call td_write_data(write_handler, gr, psi_n, h, sys%outp, sys%geo, td%dt, i) 
        end if

      end if

    end do
    
    call zoutput_function(sys%outp%how, 'opt-control', 'last_fwd', gr%m, gr%sb, psi_n%zpsi(:,1,1,1), M_ONE, ierr) 
    
    if(write_iter_) call td_write_end(write_handler)
    deallocate(dens)
    call pop_sub()
  end subroutine propagate_forward


  ! ---------------------------------------------------------
  subroutine propagate_backward(sys, h, td, laser, v_old_i, v_old_f, psi_n) ! give initial chi and laser
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(td_t),          intent(inout) :: td
    FLOAT, pointer                     :: laser(:, :)
    FLOAT, pointer                     :: v_old_i(:, :, :), v_old_f(:, :, :)
    type(states_t), intent(inout)      :: psi_n

    integer :: i
    FLOAT, allocatable :: dens(:,:)
    type(grid_t),  pointer :: gr

    call push_sub('opt_control.propagate_backward')
    
    message(1) = "Info: Backward propagating Chi"
    call write_info(1)

    gr => sys%gr

    ALLOCATE(dens(NP_PART, psi_n%d%nspin), NP_PART*psi_n%d%nspin) 

    ! setup the hamiltonian
    call states_calc_dens(psi_n, NP_PART, dens)
    psi_n%rho = dens
    call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)

    ! setup start of the propagation
    do i = 1, psi_n%d%nspin
      v_old_f(:, i, 2) = h%vhxc(:, i)
    end do
    v_old_f(:, :, 3) = v_old_f(:, :, 2)

    h%ep%lasers(1)%numerical => laser
    td%tr%v_old => v_old_f

    td%dt = -td%dt
    h%ep%lasers(1)%dt = td%dt

    do i = td%max_iter-1, 0, -1
      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), td%dt)
      ! update
      call states_calc_dens(psi_n, NP_PART, dens)
      psi_n%rho = dens
      call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
      
    end do
    td%dt = -td%dt
    
    deallocate(dens)
    call pop_sub()
  end subroutine propagate_backward
