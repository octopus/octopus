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

module run_m
  use calc_mode_m
  use casida_m
  use datasets_m
  use em_resp_m
  use epot_m
  use fft_m
  use geom_opt_m
  use global_m
  use ground_state_m
  use hamiltonian_m
  use invert_ks_m
  use parser_m
  use messages_m
  use mpi_debug_m
  use memory_m
  use multicomm_m
  use one_shot_m
  use opt_control_m
  use phonons_fd_m
  use phonons_lr_m
  use kdotp_m
  use gcm_m
  use pulpo_m
  use raman_m
  use restart_m
  use static_pol_m
  use system_m
  use td_m
  use unit_m
  use unit_system_m
  use unocc_m
  use varinfo_m
  use vdw_m

  implicit none

  private
  public ::                      &
    run_init,                    &
    run,                         &
    run_end

  type(system_t)      :: sys
  type(hamiltonian_t) :: hm

  integer, parameter :: LR = 1, FD = 2


contains

  ! ---------------------------------------------------------
  subroutine run()
    logical :: fromScratch

    ! slave nodes do not call the calculation routine
    if(multicomm_is_slave(sys%mc)) return

    PUSH_SUB(run)

    !%Variable FromScratch
    !%Type logical
    !%Default false
    !%Section Execution
    !%Description
    !% When this variable is set to true, <tt>Octopus</tt> will perform a
    !% calculation from the beginning, without looking for restart
    !% information.
    !%End

    call parse_logical(datasets_check('fromScratch'), .false., fromScratch)

    select case(calc_mode())
    case(CM_GS)
      call ground_state_run(sys, hm, fromScratch)
    case(CM_UNOCC)
      call unocc_run(sys, hm, fromScratch)
    case(CM_TD)
      call td_run(sys, hm, fromScratch)
    case(CM_LR_POL)
      select case(get_resp_method())
      case(FD)
        call static_pol_run(sys, hm, fromScratch)
      case(LR)
        call em_resp_run(sys, hm, fromScratch)
      end select
    case(CM_VDW)
      call vdW_run(sys, hm, fromScratch)
    case(CM_GEOM_OPT)
      call geom_opt_run(sys, hm, fromScratch)
    case(CM_PHONONS_LR)
      select case(get_resp_method())
      case(FD)
        call phonons_run(sys, hm)
      case(LR)
        call phonons_lr_run(sys, hm, fromscratch)
      end select
    case(CM_OPT_CONTROL)
      call opt_control_run(sys, hm)
    case(CM_CASIDA)
      call casida_run(sys, hm, fromScratch)
    case(CM_RAMAN)
      call raman_run(sys, hm, fromscratch)
    case(CM_ONE_SHOT)
      call one_shot_run(sys, hm)
    case(CM_KDOTP)
      call kdotp_lr_run(sys, hm, fromScratch)
    case(CM_MEMORY)
      call memory_run(sys, hm)
    case(CM_GCM)
      call gcm_run(sys, hm)
    case(CM_INVERTKDS)
      call invert_ks_run(sys, hm)
    case(CM_PULPO_A_FEIRA)
      call pulpo_print()
    end select

    POP_SUB(run)
    
  end subroutine run
  

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
    
    call parse_integer(datasets_check('ResponseMethod'), LR, get_resp_method)

    if(.not.varinfo_valid_option('ResponseMethod', get_resp_method)) then
      call input_error('ResponseMethod')
    end if

    POP_SUB(get_resp_method)
  end function get_resp_method
  
  ! ---------------------------------------------------------
  subroutine run_init()

    PUSH_SUB(run_init)

    call messages_print_stress(stdout, "Calculation Mode")
    call messages_print_var_option(stdout, "CalculationMode", calc_mode())
    call messages_print_stress(stdout)

    if(.not. calc_mode_is(CM_PULPO_A_FEIRA)) then
      ! initialize FFTs
      call fft_all_init()

      call unit_system_init()
      call system_init(sys)
      call hamiltonian_init(hm, sys%gr, sys%geo, sys%st, sys%ks%theory_level, sys%ks%xc_family)
      if(.not. calc_mode_is(CM_MEMORY)) then
        message(1) = "Info: Generating external potential"
        call write_info(1)
        call hamiltonian_epot_generate(hm, sys%gr, sys%geo, sys%st)
      endif
      call restart_init()
    end if

    POP_SUB(run_init)
  end subroutine run_init


  ! ---------------------------------------------------------
  subroutine run_end()
    
    PUSH_SUB(run_end)

    if(.not. calc_mode_is(CM_PULPO_A_FEIRA)) then
      call hamiltonian_end(hm, sys%gr, sys%geo)
      call system_end(sys)
      call fft_all_end()
    end if

#ifdef HAVE_MPI
    call mpi_debug_statistics()
#endif

    POP_SUB(run_end)
  end subroutine run_end

end module run_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
