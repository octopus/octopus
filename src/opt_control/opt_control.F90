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

#include "global.h"

module opt_control_m
  use excited_states_m
  use datasets_m
  use varinfo_m
  use global_m
  use filter_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lasers_m
  use lib_oct_parser_m
  use lib_oct_m
  use messages_m
  use mesh_m
  use mesh_function_m
  use output_m
  use geometry_m
  use states_m
  use states_output_m
  use string_m
  use system_m
  use td_rti_m
  use td_write_m
  use timedep_m
  use units_m
  use v_ks_m
  use external_pot_m
  use restart_m
  use tdf_m
  use mix_m
  use opt_control_constants_m
  use opt_control_propagation_m
  use opt_control_parameters_m
  use opt_control_target_m

  implicit none

  private
  public :: opt_control_run

  type oct_iterator_t
    FLOAT              :: eps
    integer            :: ctr_iter_max
    integer            :: ctr_iter
    FLOAT              :: old_functional
    FLOAT, pointer     :: convergence(:,:)
    FLOAT              :: bestJ1, bestJ1_fluence, bestJ1_J
    integer            :: bestJ1_ctr_iter
    type(oct_control_parameters_t) :: best_par
  end type oct_iterator_t

contains

  ! ---------------------------------------------------------
  subroutine opt_control_run(sys, h)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h

    type(oct_t)                :: oct
    type(oct_iterator_t)       :: iterator
    type(td_t)                 :: td
    type(grid_t),      pointer :: gr   ! some shortcuts
    type(states_t),    pointer :: st
    type(states_t)             :: chi, psi, initial_st
    type(target_t)             :: target
    type(filter_t)             :: filter
    type(oct_control_parameters_t) :: par, par_tmp
    integer :: i
    character(len=80)  :: filename


    call push_sub('opt_control.opt_control_run')

    ! Checks that the run is actually possible with the current settings
    call check_runmode_constrains(sys, h)

    call io_mkdir('opt-control')

    ! some shortcuts
    gr  => sys%gr
    st  => sys%st

    call td_init(sys, h, td)
    call states_allocate_wfns(st, gr%m, M_CMPLX)

    ! Initialize a bunch of states: initial, target, chi.
    call states_copy(chi, st)
    call states_copy(initial_st, st)

    call oct_read_inp(oct)

    ! Initial guess for the laser: read from the input file.
    call laser_init(h%ep%no_lasers, h%ep%lasers, gr%m)

    do i = 1, h%ep%no_lasers
      call laser_to_numerical(h%ep%lasers(i), td%dt, td%max_iter)
    end do

    call parameters_init(par, h%ep%no_lasers, td%dt, td%max_iter)
    call parameters_set(par, h%ep)
    call parameters_write('opt-control/initial_laser', par)
    call parameters_copy(par_tmp, par)

    call oct_iterator_init(iterator, par)

    if(oct%use_mixing) call mix_init(parameters_mix, td%max_iter + 1, par%no_parameters, 1)

    write(message(1),'(a,f14.8)') 'Input: Fluence of Initial laser ', laser_fluence(par)
    call write_info(1)

    call filter_init(td%max_iter, td%dt, filter)
    call filter_write(filter)

    call def_istate(gr, sys%geo, initial_st)
    call target_init(gr, sys%geo, sys%st, td, target)

    call check_faulty_runmodes(oct)

    call states_output(initial_st, gr, 'opt-control/initial', sys%outp)
    call target_output(gr, 'opt-control/target', sys%outp, target)

    ! psi is the "working state".
    psi = initial_st

    ! mode switcher
    select case(oct%algorithm_type)
      case(oct_algorithm_zbr98)
        message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
        call write_info(1)
        call scheme_zbr98
      case(oct_algorithm_wg05)
        message(1) = "Info: Starting OCT iteration using scheme: WG05"
        call write_info(1)
        call scheme_wg05
      case(oct_algorithm_zr98)
        message(1) = "Info: Starting OCT iteration using scheme: ZR98"
        call write_info(1)
        call scheme_mt03
      case(oct_algorithm_mt03)
        message(1) = "Info: Starting OCT iteration using scheme: MT03"
        call write_info(1)
        call scheme_mt03
      case(oct_algorithm_krotov)
        message(1) = "Info: Starting OCT iteration using scheme: KROTOV"
        call write_info(1)
        call scheme_mt03
    case default
      call input_error('OCTScheme')
    end select

    ! Some informative output.
    call output(iterator)
    call states_output(psi, gr, 'opt-control/final', sys%outp)

    ! do final test run: propagate initial state with optimal field
    call oct_finalcheck(oct, initial_st, target, iterator%best_par, sys, h, td)

    ! clean up
    call parameters_end(par)
    call parameters_end(par_tmp)
    call oct_iterator_end(iterator)
    if(oct%use_mixing) call mix_end(parameters_mix)
    call filter_end(filter)
    call td_end(td)

    nullify(gr)
    nullify(st)
    call states_end(psi)
    call states_end(chi)
    call states_end(initial_st)
    call target_end(target)
   
    call pop_sub()

  contains


    ! ---------------------------------------------------------
    subroutine scheme_zr98
      type(oct_control_parameters_t) :: par_new, par_prev
      logical :: stop_loop
      call push_sub('opt_control.scheme_zr98')

      call parameters_copy(par_new, par)
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_iter(oct, sys, h, td, iterator, psi, initial_st, target, par, stop_loop)
        call iterator_write(iterator, psi, par, gr, sys%outp)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call parameters_mixing(iterator%ctr_iter, par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
      end do ctr_loop

      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_zr98
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_mt03
      type(oct_control_parameters_t) :: par_new, par_prev
      logical :: stop_loop
      call push_sub('opt_control.scheme_mt03')

      call parameters_copy(par_new, par)
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_iter(oct, sys, h, td, iterator, psi, initial_st, target, par, stop_loop)
        call iterator_write(iterator, psi, par, gr, sys%outp)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call parameters_mixing(iterator%ctr_iter, par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
      end do ctr_loop

      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_mt03


    ! ---------------------------------------------------------
    subroutine scheme_wg05
      integer :: j
      FLOAT :: fluence, new_penalty, old_penalty

      call push_sub('opt_control.scheme_wg05')

      if (oct%mode_fixed_fluence) then
        par%alpha(1) = (M_ONE/sqrt(oct%targetfluence)) * sqrt ( laser_fluence(par) )
        par_tmp%alpha(1) = par%alpha(1)
      end if
      
      ctr_loop: do

        call states_end(psi)
        call states_copy(psi, initial_st)
        call parameters_to_h(par, h%ep)
        call propagate_forward(sys, h, td, target, psi) 

        write(filename,'(a,i3.3)') 'opt-control/PsiT.', iterator%ctr_iter
        call states_output(psi, gr, trim(filename), sys%outp)
        write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
        call parameters_write(filename, par)

        call calc_chi(oct, gr, target, psi, chi)
        
        if(iteration_manager(oct, gr, par, td, psi, target, iterator)) exit ctr_loop
        if(clean_stop()) exit ctr_loop
        
        call bwd_step(oct, sys, td, h, target, par, par_tmp, chi, psi)

        ! WARNING: Untested.
        do j = 1, par_tmp%no_parameters
          call filter_apply(par_tmp%f(j), filter)
        end do

        ! recalc field
        if (oct%mode_fixed_fluence) then
          fluence = laser_fluence(par_tmp) 
          old_penalty = par%alpha(1)
          new_penalty = sqrt( fluence * old_penalty**2 / oct%targetfluence )
          do j = 1, par_tmp%no_parameters
            par%alpha(1) = new_penalty
            par_tmp%alpha(1) = new_penalty
            call tdf_scalar_multiply( old_penalty / new_penalty , par_tmp%f(j) )
          end do
        end if

        do j = 1, par_tmp%no_parameters
          par%f(j) = par_tmp%f(j)
        end do
        
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_wg05


    ! ---------------------------------------------------------
    subroutine scheme_zbr98
      type(oct_control_parameters_t) :: par_new, par_prev

      call push_sub('opt_control.scheme_zbr98')
      
      if(oct%use_mixing) then      
        call parameters_copy(par_new, par_tmp)
        call parameters_copy(par_prev, par_tmp)
      else
        call calc_chi(oct, gr, target, psi, chi)
        call parameters_to_h(par, h%ep)
        call propagate_backward(sys, h, td, chi)
      end if

      ctr_loop: do
        if(oct%use_mixing) call parameters_copy(par_prev, par_tmp)

        if(oct%use_mixing) then
          call calc_chi(oct, gr, target, psi, chi)
          call parameters_to_h(par_tmp, h%ep)
          call propagate_backward(sys, h, td, chi)
        end if

        ! forward propagation
        call states_end(psi)
        call states_copy(psi, initial_st)
        call fwd_step(oct, sys, td, h, target, par, par_tmp, chi, psi)

        write(filename,'(a,i3.3)') 'opt-control/PsiT.', iterator%ctr_iter
        call states_output(psi, gr, filename, sys%outp)
        write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
        call parameters_write(filename, par)
        
        if(iteration_manager(oct, gr, par, td, psi, target, iterator)) exit ctr_loop
        if(clean_stop()) exit ctr_loop

        ! and now backward
        call calc_chi(oct, gr, target, psi, chi)
        call bwd_step(oct, sys, td, h, target, par, par_tmp, chi, psi)

        if(oct%use_mixing) then
          call parameters_mixing(iterator%ctr_iter, par_prev, par_tmp, par_new)
          call parameters_copy(par_tmp, par_new)
        end if

      end do ctr_loop

      if(oct%use_mixing) then
        call parameters_end(par_new)
        call parameters_end(par_prev)
      end if
      call pop_sub()
    end subroutine scheme_zbr98
    ! ---------------------------------------------------------
 
  end subroutine opt_control_run


  ! ---------------------------------------------------------
  subroutine f_iter(oct, sys, h, td, iterator, psi, initial_st, target, par, stop_loop)
    type(oct_t), intent(in)                       :: oct
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: h
    type(td_t), intent(inout)                     :: td
    type(oct_iterator_t), intent(inout)           :: iterator
    type(states_t), intent(inout)                 :: psi
    type(states_t), intent(in)                    :: initial_st
    type(target_t), intent(inout)                 :: target
    type(oct_control_parameters_t), intent(inout) :: par
    logical, intent(out)                          :: stop_loop

    type(states_t) :: chi
    type(oct_control_parameters_t) :: parp

    call push_sub('opt_control.f_zbr98')

    call parameters_copy(parp, par)

    if( (iterator%ctr_iter .eq. 0) .or. oct%use_mixing) then
      call states_end(psi)
      psi = initial_st
      call parameters_to_h(par, h%ep)
      call propagate_forward(sys, h, td, target, psi)
    end if

    stop_loop = iteration_manager(oct, sys%gr, par, td, psi, target, iterator)

    chi = psi
    call calc_chi(oct, sys%gr, target, psi, chi)
    call bwd_step(oct, sys, td, h, target, par, parp, chi, psi)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call fwd_step(oct, sys, td, h, target, par, parp, chi, psi)

    call states_end(chi)
    call parameters_end(parp)
    call pop_sub()
  end subroutine f_iter
  ! ---------------------------------------------------------



  ! ---------------------------------------------------------
  ! This subroutine just stops the run if some of the settings
  ! are not compatible with the run mode, either because it is
  ! meaningless, or because the run mode is still not fully
  ! implemented.
  !
  ! TODO: Right now just a couple of checks are made, but there are many other constrains.
  subroutine check_runmode_constrains(sys, h)
    type(system_t), target, intent(in) :: sys
    type(hamiltonian_t),    intent(in) :: h

    integer :: no_electrons, n_filled, n_partially_filled, n_half_filled

    ! Only dipole approximation in length gauge.
    if(h%gauge.ne.LENGTH) then
      write(message(1),'(a)') "So far only length gauge is supported in optimal control runs."
      call write_fatal(1)
    end if

    ! This should check that we really have occupation one for
    ! one of the spin-orbitals, and occupation zero for all the others.
    ! Otherwise the algorithms are bound to fail.
    select case(sys%st%d%ispin)
    case(UNPOLARIZED)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      no_electrons = 2*n_filled + n_half_filled
      if(n_partially_filled > 0 ) then
        write(message(1),'(a)') 'No partially filled orbitals are allowd in OCT calculations'
        call write_fatal(1)
      end if
    case(SPIN_POLARIZED)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowd in OCT calculations'
        call write_fatal(1)
      end if
      no_electrons = n_filled
      call occupied_states(sys%st, 2, n_filled, n_partially_filled, n_half_filled)
      no_electrons = n_filled + no_electrons
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowd in OCT calculations'
        call write_fatal(1)
      end if
    case(SPINORS)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      no_electrons = n_filled
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowd in OCT calculations'
        call write_fatal(1)
      end if
    end select

    if(abs(-sys%st%val_charge - real(no_electrons, REAL_PRECISION) ) > CNST(1.0e-8)) then
      write(message(1), '(a)') 'Error in check_runmode_constrains'
      call write_fatal(1)
    end if

  end subroutine check_runmode_constrains

#include "read.F90"
#include "aux.F90"
#include "defstates.F90"
#include "finalcheck.F90"
#include "iter.F90"
#include "output.F90"

end module opt_control_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
