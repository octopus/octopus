
!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module run_prog
  use global
  use messages
  use lib_oct_parser
  use units
  use fft
  use restart
  use external_pot
  use hamiltonian
  use system
  use scf
  use timedep
  use unocc
  use static_pol
  use static_pol_lr
  use casida
  use wave_matching
  use geom_opt
  use phonons
  use opt_control
  use ground_state
  use pulpo
  use io
  use multicomm_mod
  use mpi_debug_mod

  implicit none

  private
  public ::                      &
    run_init,                    &
    run,                         &
    run_end

  integer, public, parameter ::  &
    M_GS                 =   1,  &
    M_UNOCC              =   2,  &
    M_TD                 =   3,  &
    M_STATIC_POL         =   4,  &
    M_GEOM_OPT           =   5,  &
    M_PHONONS            =   6,  &
    M_OPT_CONTROL        =   7,  &
    M_LR_STATIC_POL      =   8,  &
    M_CASIDA             =   9,  &
    M_WAVE_MATCHING      =  10,  &
    M_BO_MD              =  98,  &
    M_PULPO_A_FEIRA      =  99

  type(system_type)      :: sys
  type(hamiltonian_type) :: h


contains

  ! ---------------------------------------------------------
  subroutine run()
    logical :: fromScratch

    call push_sub('run.run')

    call loct_parse_logical(check_inp('fromScratch'), .false., fromScratch)

    select case(calc_mode)
    case(M_GS)
      call ground_state_run(sys, h, fromScratch)
    case(M_UNOCC)
      call unocc_run(sys, h, fromScratch)
    case(M_TD)
      call td_run(sys, h, fromScratch)
    case(M_STATIC_POL)
      call static_pol_run(sys, h, fromScratch)
    case(M_LR_STATIC_POL)
      call static_pol_lr_run(sys, h, fromScratch)
    case(M_GEOM_OPT)
      call geom_opt_run(sys, h)
    case(M_PHONONS)
      call phonons_run(sys, h)
    case(M_OPT_CONTROL)
      call opt_control_run(sys, h)
    case(M_CASIDA)
      call casida_run(sys, h, fromScratch)
    case(M_WAVE_MATCHING)
      call wave_matching_run()
    case(M_PULPO_A_FEIRA)
      call pulpo_print()
    end select

    call pop_sub()

  end subroutine run


  ! ---------------------------------------------------------
  subroutine run_init()
    integer :: mode

    call messages_print_stress(stdout, "Calculation Mode")
    call messages_print_var_option(stdout, "CalculationMode", calc_mode)
    call messages_print_stress(stdout)

    ! initialize ffts
#ifdef HAVE_FFT
    call fft_all_init()
#endif

    if(calc_mode .ne. M_PULPO_A_FEIRA) then
      call units_init()
      call get_mode_parallelized(mode)

      call system_init(sys, mode)
      call hamiltonian_init(h, sys%gr, sys%st%d, sys%ks%ip_app)
      call epot_generate(h%ep, sys%gr%m, sys%gr%sb, sys%gr%geo, sys%st, h%reltype)

      call restart_init()
    end if


  contains

    subroutine get_mode_parallelized(mode)
      integer, intent(out) :: mode

      mode = ibset(mode, P_STRATEGY_DOMAINS - 1) ! all modes are parallel in domains
      select case(calc_mode)
      case(M_TD)
        mode = ibset(mode, P_STRATEGY_STATES - 1)
      case(M_CASIDA)
        mode = ibset(mode, P_STRATEGY_OTHER - 1)
      end select

    end subroutine get_mode_parallelized

  end subroutine run_init


  ! ---------------------------------------------------------
  subroutine run_end()
    if(calc_mode .ne. M_PULPO_A_FEIRA) then
       call hamiltonian_end(h, sys%gr)
       call system_end(sys)
    end if

#ifdef HAVE_FFT
    call fft_all_end()
#endif

#ifdef HAVE_MPI
    call MPI_Debug_Statistics()
#endif

  end subroutine run_end

end module run_prog
