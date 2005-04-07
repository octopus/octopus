
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

#include "global.h"

module run_prog
use global
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
use geom_opt
use phonons
use opt_control
use ground_state
use pulpo

implicit none

private
public :: run_init, &
          run,      &
          run_end

type(system_type) :: sys
type(hamiltonian_type) :: h
integer :: calc_mode

! run stack
integer, private :: i_stack(100), instr
integer, private, parameter ::   &
     M_GS                 = 1, &
     M_UNOCC              = 2, &
     M_TD                 = 3, &
     M_STATIC_POL         = 4, &
     M_GEOM_OPT           = 5, &
     M_PHONONS            = 6, &
     M_OPT_CONTROL        = 7, &
     M_LR_STATIC_POL      = 8, &
     M_CASIDA             = 9, &
     M_BO_MD              = 98,&
     M_PULPO_A_FEIRA      = 99

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
     I_PULPO              = 999

contains

subroutine run()
  integer :: i, ierr
  FLOAT :: x
  logical :: log
  character(len=100) :: filename

  logical :: fromScratch(M_GS:M_CASIDA)

  call push_sub('run')

  instr = 0
  call define_run_modes()

  program: do
    if(instr <= 0) exit program

    select case(i_stack(instr))
    case(I_GS_INIT)
      call ground_state_init(sys, h)

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

    call loct_parse_logical("fromScratch", .false., fS)
    fromScratch(:) = .false.

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

    case(M_PULPO_A_FEIRA)
      instr = instr + 1; i_stack(instr) = I_PULPO
  end select

end subroutine define_run_modes

end subroutine run

subroutine run_init()
  ! initialize some stuff

  call loct_parse_int('CalculationMode', 1, calc_mode)
  if( (calc_mode < 1 .or. calc_mode > 10) .and. (calc_mode .ne. M_PULPO_A_FEIRA)) then
    write(message(1), '(a,i2,a)') "Input: '", calc_mode, "' is not a valid CalculationMode"
    message(2) = '  Calculation Mode = '
    message(3) = '    gs          <= ground-state calculation'
    message(4) = '    unocc       <= calculate unocuppied states'
    message(5) = '    td          <= time-dependent simulation'
    message(6) = '    pol         <= calculate static polarizability'
    message(6) = '    pol_lr      <= calculate static polarizability from LR theory'
    message(7) = '    bo          <= perform Born-Oppenheimer MD'
    message(8) = '    geom        <= geometry optimization'
    message(9) = '    phonon      <= calculate phonon frequencies'
    message(10)= '    opt_control <= optimum control'
    message(11)= '    casida      <= calculate excitation a la Marc Casida'
    message(12)= '    recipe      <= prints out the "Pulpo a Feira" recipe'
    call write_fatal(12)
  end if

  write(message(1), '(a,i2)')   'Info: Calculation Mode = ', calc_mode
  write(message(2), '(a,i1,a)') 'Info: The octopus will run in ', conf%dim, ' dimension(s).'
  write(message(3), '(a,i1,a)') '      The octopus will treat system as periodic in ', &
                                 conf%periodic_dim, ' dimension(s)'

  message(5) = "Info: Boundary conditions:"
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
    call hamiltonian_init(h, sys%m, sys%geo, sys%st%d, sys%ks%ip_app)
    call epot_generate(h%ep, sys%m, sys%st, sys%geo, h%reltype)
  endif

  call restart_init()

end subroutine run_init

subroutine run_end()
  if(calc_mode .ne. M_PULPO_A_FEIRA) then
    call hamiltonian_end(h, sys%geo)
    call system_end(sys)
  endif

#ifdef HAVE_FFT
  call fft_all_end()
#endif

end subroutine run_end

end module run_prog
