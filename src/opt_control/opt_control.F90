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

  type oct_t
    FLOAT   :: targetfluence
    logical :: mode_fixed_fluence
    integer :: algorithm_type
    logical :: use_mixing
    logical :: oct_double_check
    logical :: dump_intermediate
  end type oct_t

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
    type(filter_t),    pointer :: f(:)
    type(oct_control_parameters_t) :: par, par_tmp
    integer :: i, ierr
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

    call oct_iterator_init(iterator, oct, par)

    if(oct%use_mixing) call mix_init(parameters_mix, td%max_iter + 1, par%no_parameters, 1)

    write(message(1),'(a,f14.8)') 'Input: Fluence of Initial laser ', laser_fluence(par)
    call write_info(1)

    nullify(f)
    call def_filter(td%max_iter, td%dt, f)
    call filter_write(f, NDIM, td%dt, td%max_iter)

    call def_istate(gr, sys%geo, initial_st)
    call target_init(gr, sys%geo, sys%st, td, target)

    call check_faulty_runmodes(oct)

    call states_output(initial_st, gr, 'opt-control/initial', sys%outp)
    call target_output(gr, 'opt-control/target', sys%outp, target)

    ! psi is the "working state".
    psi = initial_st

    ! mode switcher
    select case(oct%algorithm_type)
      case(oct_algorithm_zbr98) ;  call scheme_zbr98
      case(oct_algorithm_zr98)  ;  call scheme_zr98
      case(oct_algorithm_wg05)  ;  call scheme_wg05
    case default
      write(message(1),'(a)') "The OCT algorithm has not been implemented yet."
      write(message(2),'(a)') "Choose: ZR98, ZBR98, WG05."
      call write_fatal(2)
    end select

    ! Some informative output.
    call output(oct, iterator)
    call states_output(psi, gr, 'opt-control/final', sys%outp)

    ! do final test run: propagate initial state with optimal field
    call oct_finalcheck(oct, initial_st, target, iterator%best_par, sys, h, td)

    ! clean up
    call parameters_end(par)
    call parameters_end(par_tmp)
    call oct_iterator_end(iterator)
    if(oct%use_mixing) call mix_end(parameters_mix)
    call filter_end(f)
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
      integer :: ierr
      type(oct_control_parameters_t) :: par_new, par_prev
      call push_sub('opt_control.scheme_ZR98')

      message(1) = "Info: Starting Optimal Control Run  using Scheme: ZR98"
      call write_info(1)
      
      if(oct%use_mixing) then
        call parameters_copy(par_new, par)
        call parameters_copy(par_prev, par)
      else
        psi = initial_st
        call parameters_to_h(par, h%ep)
        call propagate_forward(sys, h, td, target, psi) 
        write(filename,'(a,i3.3)') 'opt-control/PsiT.', iterator%ctr_iter
        call states_output(psi, gr, filename, sys%outp)
      end if
      
      ctr_loop: do
        if(oct%use_mixing) call parameters_copy(par_prev, par)

        if(oct%use_mixing) then
          call states_end(psi)
          psi = initial_st
          call parameters_to_h(par, h%ep)
          call propagate_forward(sys, h, td, target, psi)
        end if

        if(iteration_manager(oct, gr, par, td, psi, target, iterator)) exit ctr_loop
        if(clean_stop()) exit ctr_loop
        
        call calc_chi(oct, gr, target, psi, chi)
        call bwd_step(oct_algorithm_zr98, sys, td, h, target, par, par_tmp, chi, psi)
        
        ! forward propagation
        call states_end(psi)
        call states_copy(psi, initial_st)
        call fwd_step(oct_algorithm_zr98, sys, td, h, target, par, par_tmp, chi, psi)

        if(oct%use_mixing) then
          call parameters_mixing(iterator%ctr_iter, par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
 
        ! Print state and laser after the forward propagation.
        write(filename,'(a,i3.3)') 'opt-control/PsiT.', iterator%ctr_iter
        call states_output(psi, gr, filename, sys%outp)
        write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
        call parameters_write(filename, par)
        
      end do ctr_loop

      if(oct%use_mixing) then
        call parameters_end(par_new)
        call parameters_end(par_prev)
      end if
      call pop_sub()
    end subroutine scheme_zr98
    

    !---------------------------------------
    subroutine scheme_wg05
      integer :: ierr, j
      FLOAT :: fluence, new_penalty, old_penalty
      call push_sub('opt_control.scheme_WG05')
      
      message(1) = "Info: Starting OCT iteration using scheme: WG05"
      call write_info(1)
      
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
        
        call bwd_step(oct_algorithm_wg05, sys, td, h, target, par, par_tmp, chi, psi)

        ! WARNING: Untested.
        do j = 1, par_tmp%no_parameters
          call apply_filter(td%max_iter, f, par_tmp%f(j))
        end do

        ! recalc field
        if (oct%mode_fixed_fluence) then
          fluence = laser_fluence(par_tmp)
          old_penalty = tdf(par%td_penalty(1), 1)
          new_penalty = sqrt( fluence * old_penalty**2 / oct%targetfluence )

          write(message(1), '(a,e15.6)') 'current penalty =', old_penalty
          write(message(2), '(a,e15.6)') 'new penalty     =', new_penalty
          call write_info(2)

          do j = 1, par_tmp%no_parameters
            call tdf_set_numerical(par%td_penalty(j), &
                                   spread(new_penalty, 1, td%max_iter+1) )
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
      integer :: ierr
      type(oct_control_parameters_t) :: par_new, par_prev

!!!!DEBUG
      FLOAT :: etime
!!!!ENDOFDEBUG

      call push_sub('opt_control.scheme_zbr98')
      
      message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
      call write_info(1)

      if(oct%use_mixing) then      
        call parameters_copy(par_new, par_tmp)
        call parameters_copy(par_prev, par_tmp)
      else
        call calc_chi(oct, gr, target, psi, chi)
        call parameters_to_h(par, h%ep)
        etime =  loct_clock()
        call propagate_backward(sys, h, td, chi)
        etime = loct_clock() - etime
        write(0, *) 'PROPAGATE_BACKWARD OUTPUT', etime
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
        etime = loct_clock()
        call fwd_step(oct_algorithm_zbr98, sys, td, h, target, par, par_tmp, chi, psi)
        etime = loct_clock() - etime
        write(0, *) 'FWD_STEP OUTPUT', etime

        write(filename,'(a,i3.3)') 'opt-control/PsiT.', iterator%ctr_iter
        call states_output(psi, gr, filename, sys%outp)
        write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
        call parameters_write(filename, par)
        
        if(iteration_manager(oct, gr, par, td, psi, target, iterator)) exit ctr_loop
        if(clean_stop()) exit ctr_loop

        ! and now backward
        call calc_chi(oct, gr, target, psi, chi)
        etime = loct_clock()
        call bwd_step(oct_algorithm_zbr98, sys, td, h, target, par, par_tmp, chi, psi)
        etime = loct_clock() - etime
        write(0, *) 'BWD_STEP OUTPUT', etime

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

 
  end subroutine opt_control_run


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

    integer :: no_electrons

    ! Only dipole approximation in length gauge.
    if(h%gauge.ne.LENGTH) then
      write(message(1),'(a)') "So far only length gauge is supported in optimal control runs."
      call write_fatal(1)
    end if

    ! Only single-electron calculations.
    ! WARNING: maybe this should check that we really have occupation one for 
    ! one of the spin-orbitals, and occupation zero for all the others. 
    ! Otherwise the algorithms are bound to fail.
    no_electrons = nint(sum(sys%st%occ(:, :)))
    if(no_electrons .ne. 1) then
      write(message(1),'(a,i4,a)') 'Number of electrons (currently ', &
        no_electrons,') should be just one.'
      write(message(2),'(a)') 'Optimal control theory for many electron systems not yet developed!'
      call write_fatal(2)
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
