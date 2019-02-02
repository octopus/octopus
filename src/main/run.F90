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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module run_oct_m
  use fft_oct_m
  use global_oct_m
  use ground_state_oct_m
  use hamiltonian_oct_m
  use parser_oct_m
  use messages_oct_m
  use mpi_debug_oct_m
  use multicomm_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use restart_oct_m
  use system_oct_m
  use td_oct_m
  use test_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                      &
    run

  integer :: calc_mode_id

  integer, parameter :: LR = 1, FD = 2

  integer, public, parameter ::   &
    CM_NONE               =   0,  &
    CM_GS                 =   1,  &
    CM_TD                 =   3,  &
    CM_TEST               =  19

contains

  ! ---------------------------------------------------------
  integer function get_resp_method()

    PUSH_SUB(get_resp_method)
    
    !%Variable ResponseMethod
    !%Type integer
    !%Default sternheimer
    !%Section Linear Response
    !%Description
    !% Some response properties can be calculated either via
    !% Sternheimer linear response or by using finite
    !% differences. You can use this variable to select how you want
    !% them to be calculated, it applies to <tt>em_resp</tt> and <tt>vib_modes</tt>
    !% calculation modes. By default, the Sternheimer linear-response
    !% technique is used.
    !%Option sternheimer 1
    !% The linear response is obtained by solving a self-consistent
    !% Sternheimer equation for the variation of the orbitals. This
    !% is the recommended method.
    !%Option finite_differences 2
    !% Properties are calculated as a finite-differences derivative of
    !% the energy obtained by several ground-state calculations. This
    !% method, slow and limited only to static response, is kept
    !% mainly because it is simple and useful for testing purposes.
    !%End
    
    call parse_variable('ResponseMethod', LR, get_resp_method)

    if(.not.varinfo_valid_option('ResponseMethod', get_resp_method)) then
      call messages_input_error('ResponseMethod')
    end if

    POP_SUB(get_resp_method)
  end function get_resp_method
  
  ! ---------------------------------------------------------
  subroutine run(cm)
    integer, intent(in) :: cm

    type(system_t)      :: sys
    type(hamiltonian_t) :: hm
    type(profile_t), save :: calc_mode_prof
    logical :: fromScratch

    PUSH_SUB(run)

    calc_mode_id = cm

    call messages_print_stress(stdout, "Calculation Mode")
    call messages_print_var_option(stdout, "CalculationMode", calc_mode_id)
    call messages_print_stress(stdout)

    call calc_mode_init()

    call restart_module_init()

    ! initialize FFTs
    call fft_all_init()

    call unit_system_init()

    if(calc_mode_id == CM_TEST) then
      call test_run()
      call fft_all_end()
#ifdef HAVE_MPI
      call mpi_debug_statistics()
#endif
      POP_SUB(run)
      return
    end if

    call system_init(sys)

    call hamiltonian_init(hm, sys%gr, sys%geo, sys%st, sys%ks%theory_level, &
      sys%ks%xc_family, sys%ks%xc_flags, &
      family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
    
    message(1) = "Info: Generating external potential"
    call messages_info(1)
    call hamiltonian_epot_generate(hm, sys%gr, sys%geo, sys%st)
    message(1) = "      done."
    call messages_info(1)
    
    if(sys%ks%theory_level /= INDEPENDENT_PARTICLES) then
      call poisson_async_init(sys%ks%hartree_solver, sys%mc)
      ! slave nodes do not call the calculation routine
      if(multicomm_is_slave(sys%mc))then
        !for the moment we only have one type of slave
        call poisson_slave_work(sys%ks%hartree_solver)
      end if
    end if

    if(.not. multicomm_is_slave(sys%mc)) then
      call messages_write('Info: Octopus initialization completed.', new_line = .true.)
      call messages_write('Info: Starting calculation mode.')
      call messages_info()

      !%Variable FromScratch
      !%Type logical
      !%Default false
      !%Section Execution
      !%Description
      !% When this variable is set to true, <tt>Octopus</tt> will perform a
      !% calculation from the beginning, without looking for restart
      !% information.
      !%End

      call parse_variable('FromScratch', .false., fromScratch)

      call profiling_in(calc_mode_prof, "CALC_MODE")

      select case(calc_mode_id)
      case(CM_GS)
        call ground_state_run(sys, hm, fromScratch)
      case(CM_TD)
        if(sys%gr%sb%kpoints%use_symmetries) &
          call messages_experimental("KPoints symmetries with CalculationMode = td")
        call td_run(sys, hm, fromScratch)
      end select

      call profiling_out(calc_mode_prof)
    end if

    if(sys%ks%theory_level /= INDEPENDENT_PARTICLES) call poisson_async_end(sys%ks%hartree_solver, sys%mc)
    
    call hamiltonian_end(hm)
    call system_end(sys)

    call fft_all_end()

#ifdef HAVE_MPI
    call mpi_debug_statistics()
#endif

    POP_SUB(run)

  contains

    subroutine calc_mode_init()

      PUSH_SUB(calc_mode_init)

      select case(calc_mode_id)
      case(CM_GS)
        call ground_state_run_init()
      case(CM_TD)
        call td_run_init()
      end select

      POP_SUB(calc_mode_init)
    end subroutine calc_mode_init

  end subroutine run

end module run_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
