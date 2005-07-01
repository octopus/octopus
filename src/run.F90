
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

implicit none

private
public :: run_init, &
          run,      &
          run_end

type(system_type) :: sys
type(hamiltonian_type) :: h


! run stack
integer, private :: i_stack(100), instr
integer, private, parameter ::   &
     M_GS                 =  1, &
     M_UNOCC              =  2, &
     M_TD                 =  3, &
     M_STATIC_POL         =  4, &
     M_GEOM_OPT           =  5, &
     M_PHONONS            =  6, &
     M_OPT_CONTROL        =  7, &
     M_LR_STATIC_POL      =  8, &
     M_CASIDA             =  9, &
     M_WAVE_MATCHING      = 10, &
     M_BO_MD              = 98, &
     M_PULPO_A_FEIRA      = 99, &
     M_MULTI_SUBSYSTEM    = multi_subsys_mode

integer, private, parameter :: &
     I_GS_INIT            = 100, &
     I_GS                 = 101, &
     I_UNOCC              = 102, &
     I_TD                 = 103, &
     I_STATIC_POL         = 104, &
     I_GEOM_OPT           = 105, &
     I_PHONONS            = 106, &
     I_OPT_CONTROL        = 107, &
     I_LR_STATIC_POL      = 108, &
     I_CASIDA             = 109, &
     I_WAVE_MATCHING      = 110, &
     I_PULPO              = 999

contains

subroutine run()
  integer :: ierr

  logical :: fromScratch(M_GS:M_WAVE_MATCHING)

  call push_sub('run')

  instr = 0
  call define_run_modes()

  program: do
    if(instr <= 0) exit program

    select case(i_stack(instr))
    case(I_GS_INIT)
      call ground_state_init(sys%gr, sys%st, sys%ks, h)

    case(I_GS)
      if(ground_state_run(sys, h).ne.0) then ! could not load wfs
        i_stack(instr) = I_GS;      instr = instr + 1
        i_stack(instr) = I_GS_INIT
        cycle program
      end if

    case(I_UNOCC)
      if(unocc_run(sys, h, fromScratch(M_UNOCC)).ne.0) then
        i_stack(instr) = I_UNOCC;    instr = instr + 1
        i_stack(instr) = I_GS
        cycle program
      end if

    case(I_TD)
      if(td_run(sys, h, fromScratch(M_TD)).ne.0) then
        i_stack(instr) = I_TD;    instr = instr + 1
        i_stack(instr) = I_GS
        cycle program
      end if

    case(I_STATIC_POL)
      if(static_pol_run(sys, h, fromScratch(M_STATIC_POL)).ne.0) then ! could not load wfs
        i_stack(instr) = I_STATIC_POL;      instr = instr + 1
        i_stack(instr) = I_GS_INIT
        cycle program
      end if

    case(I_GEOM_OPT)
      if(geom_opt_run(sys, h).ne.0) then ! could not load wfs
        i_stack(instr) = I_GEOM_OPT;      instr = instr + 1
        i_stack(instr) = I_GS_INIT
        cycle program
      end if

    case(I_PHONONS)
      if(phonons_run(sys, h).ne.0) then ! could not load wfs
        i_stack(instr) = I_PHONONS;      instr = instr + 1
        i_stack(instr) = I_GS_INIT
        cycle program
      end if

    case(I_OPT_CONTROL)
      message(1) = 'Info: Optimum control.'
      call write_info(1)

      ierr = opt_control_run(sys, h)

    case(I_LR_STATIC_POL)
      if(static_pol_lr_run(sys, h, fromScratch(M_LR_STATIC_POL)).ne.0) then ! could not load wfs
        i_stack(instr) = I_LR_STATIC_POL;      instr = instr + 1
        i_stack(instr) = I_GS
        cycle program
      end if

    case(I_CASIDA)
      if(casida_run(sys, h, fromScratch(M_CASIDA)).ne.0) then ! could not load wfs
        i_stack(instr) = I_CASIDA;      instr = instr + 1
        i_stack(instr) = I_UNOCC
        cycle program
      end if

    case(I_WAVE_MATCHING)
      ierr = wave_matching_run(sys, h)

    case(I_PULPO)
      call pulpo_print()

    case default
      write(message(1), '(a,i3,a)') "Instruction ", i_stack(instr), " not defined!"
      call write_warning(1)
    end select

    instr = instr - 1
  end do program
      
  call pop_sub()

contains

  subroutine define_run_modes()
    logical :: fS

    call loct_parse_logical(check_inp('fromScratch'), .false., fS)
    fromScratch(:) = .false.

#if defined(HAVE_MPI)
    if((calc_mode.ne.M_TD).and.(calc_mode.ne.M_CASIDA)) then
      message(1) = "Code is only parallelized for run modes 'td' and 'casida'"
      message(2) = "Please use the serial version for the current run mode"
      call write_fatal(2)
    end if
#endif

    select case(calc_mode)
    case(M_GS)
      fromScratch(M_GS) = fS
      instr = instr + 1; i_stack(instr) = I_GS
      if(fS) then
        instr = instr + 1; i_stack(instr) = I_GS_INIT
      end if

    case(M_UNOCC)
      fromScratch(M_UNOCC) = fS
      instr = instr + 1; i_stack(instr) = I_UNOCC

    case(M_TD)
      fromScratch(M_TD) = fS
      instr = instr + 1; i_stack(instr) = I_TD

    case(M_STATIC_POL)
      fromScratch(M_STATIC_POL) = fS
      instr = instr + 1; i_stack(instr) = I_STATIC_POL

    case(M_LR_STATIC_POL)
      fromScratch(M_LR_STATIC_POL) = fS
      instr = instr + 1; i_stack(instr) = I_LR_STATIC_POL

    case(M_GEOM_OPT)
      fromScratch(M_GEOM_OPT) = fS
      instr = instr + 1; i_stack(instr) = I_GEOM_OPT

    case(M_PHONONS)
      fromScratch(M_PHONONS) = fS
      instr = instr + 1; i_stack(instr) = I_PHONONS

    case(M_OPT_CONTROL)
      fromScratch(M_OPT_CONTROL) = fS
      instr = instr + 1; i_stack(instr) = I_OPT_CONTROL

    case(M_CASIDA)
      fromScratch(M_CASIDA) = fS
      instr = instr + 1; i_stack(instr) = I_CASIDA

    case(M_WAVE_MATCHING)
      fromScratch(M_WAVE_MATCHING) = fS
      instr = instr + 1; i_stack(instr) = I_WAVE_MATCHING

    case(M_PULPO_A_FEIRA)
      instr = instr + 1; i_stack(instr) = I_PULPO
  end select

end subroutine define_run_modes

end subroutine run

subroutine run_init()
  character(len=70) :: mode_string

  !%Variable Dimensions
  !%Type integer
  !%Section 1 Generalities
  !%Description
  !% octopus can run in 1, 2 or 3 dimensions, depending on the value of this
  !% variable. Note that not all input variables may be available in all cases.
  !%Option 1
  !% The system is 1-dimensional
  !%Option 2
  !% The system is 2-dimensional
  !%Option 3
  !% The system is 3-dimensional (default)
  !%End
  call loct_parse_int(check_inp('Dimensions'), 3, conf%dim)
  if(conf%dim<1 .or. conf%dim>3) call input_error('Dimensions')

  ! handle periodic directions
  call loct_parse_int(check_inp('PeriodicDimensions'), 0, conf%periodic_dim)
  if ((conf%periodic_dim < 0) .or. (conf%periodic_dim > 3)) then
    message(1) = 'Periodic dimensions must be either 0, 1, 2, or 3'
    call write_fatal(1)
  endif
  if(conf%periodic_dim > conf%dim) then
    message(1) = 'PeriodicDimensions must be <= Dimensions'
    call write_fatal(1)
  end if

  call loct_parse_logical(check_inp('BoundaryZeroDerivative'), .false., conf%boundary_zero_derivative)

  ! do we treat only userdefined species
  call loct_parse_logical(check_inp('OnlyUserDef'), .false., conf%only_user_def)

  select case(calc_mode)
    case(M_GS);              mode_string = "gs [Calculation of the ground state]"
    case(M_UNOCC);           mode_string = "unocc [Calculation of unoccupied/virtual KS states]"
    case(M_TD);              mode_string = "td [Time-dependent calculation]"
    case(M_STATIC_POL);      mode_string = "pol [Calculation of the static polarizability]"
    case(M_GEOM_OPT);        mode_string = "geom [Optimization of the geometry]"
    case(M_PHONONS);         mode_string = "phonons [Calculation of the vibrational modes]"
    case(M_OPT_CONTROL);     mode_string = "opt_control [Optimal control]"
    case(M_LR_STATIC_POL);   mode_string = "pol_lr [Linear-response calculation of the polarizability]"
    case(M_CASIDA);          mode_string = "casida [Excitations via linear-response TDDFT]"
    case(M_WAVE_MATCHING);   mode_string = "wave_matching [Wave-matching a la Heiko]"
    case(M_BO_MD);           mode_string = "bo [Born-Oppenheimer-like Molecular Dynamics]"
    case(M_PULPO_A_FEIRA);   mode_string = "recipe [Prints out a tasty recipe]"
    case(M_MULTI_SUBSYSTEM); mode_string = "multi_subsystem [Multi subsystem mode]"
    case default; call input_error('CalculationMode')
  end select
  write(message(1), '(a,a)')    'Calculation Mode = ', trim(mode_string)
  write(message(2), '(a,i1,a)') 'The octopus will run in ', conf%dim, ' dimension(s).'
  write(message(3), '(a,i1,a)') 'The octopus will treat system as periodic in ', &
                                 conf%periodic_dim, ' dimension(s)'

  message(5) = "Boundary conditions:"
  if(conf%boundary_zero_derivative) then
    write(message(4), '(2a)') trim(message(5)), " zero derivatives"
  else
    write(message(4), '(2a)') trim(message(5)), " zero wave-functions"
  end if

  call write_info(4, stress = .true.)

  ! initialize ffts
#ifdef HAVE_FFT
  call fft_all_init()
#endif

  if(calc_mode .ne. M_PULPO_A_FEIRA) then
    call units_init()
    call system_init(sys)
    call hamiltonian_init(h, sys%gr%m, sys%gr%geo, sys%st%d, sys%ks%ip_app)
    call epot_generate(h%ep, sys%gr%m, sys%gr%sb, sys%gr%geo, sys%st, h%reltype)
  endif

  call restart_init()

end subroutine run_init

subroutine run_end()
  if(calc_mode .ne. M_PULPO_A_FEIRA) then
    call hamiltonian_end(h, sys%gr%geo)
    call system_end(sys)
  endif

#ifdef HAVE_FFT
  call fft_all_end()
#endif

end subroutine run_end

end module run_prog
