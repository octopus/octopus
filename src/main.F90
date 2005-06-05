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

program octopus
  use string
  use global
  use messages
  use syslabels
  use lib_oct
  use lib_oct_parser
  use run_prog
  use io

  implicit none

  integer :: ns
  character(len=256) :: sys_name
#if defined(HAVE_MPI)
  integer :: ierr
#endif

  call global_init()
  call parser_init()
  ! initialize input/output system
  call io_init()

  ! need to find out calc_mode already here since some of the variables here (e.g.
  ! periodic dimensions) can be different for the subsystems
  !%Variable CalculationMode
  !%Type integer
  !%Section 1 Generalities
  !%Description
  !% Decides what kind of calculation is to be performed
  !%Option gs 1
  !% Calculation of the ground state
  !%Option unocc 2
  !% Calculation of unoccupied/virtual KS states
  !%Option td 3
  !% Time-dependent calculation
  !%Option pol 4
  !% Calculation of the static polarizability
  !%Option geom 5
  !% Optimization of the geometry
  !%Option phonons 6
  !% Calculation of the vibrational modes
  !%Option opt_control 7
  !% Optimal control.
  !%Option pol_lr 8
  !% Linear-response calculation of the polarizability
  !%Option casida 9
  !% Excitations via linear-response TDDFT
  !%Option wave_matching 10
  !% Wave-matching a la Heiko
  !%Option bo 98
  !% Born-Oppenheimer-like Molecular Dynamics
  !%Option recipe 99
  !% Prints out a tasty recipe
  !%Option multi_subsystem 1000
  !% Multi subsystem mode
  !%End
  call loct_parse_int('CalculationMode', 1, calc_mode)
  if(calc_mode == multi_subsys_mode) then
     call read_system_labels()
  else
     call syslabels_init(calc_mode)
  endif


  ! loop over all subsystems
  subsystems: do ns = 1,no_syslabels

     ! set system label
     current_label = trim(subsys_label(subsys_run_order(ns)))
     current_subsystem = subsys_run_order(ns)
     calc_mode = subsys_runmode(subsys_run_order(ns))

     ! Let us print our logo
     if(mpiv%node == 0) then
        ! Let us print our logo
        if(conf%verbose >= VERBOSE_NORMAL) call io_dump_file(stdout, trim(trim(conf%share) // '/logo'))
#ifdef DEBUG
        if(conf%verbose >= VERBOSE_DEBUG) write(stderr, '(5a)') "# ", " A ", "Time", "Mem", "Call"
#endif
     end if

     ! Let us print the version
     message(1) = ""
     message(2) = str_center("Running octopus, version " // trim(conf%version), 70)
     message(3) = str_center("(build time - " // trim(conf%build_time) // ")", 70)
     message(4) = str_center("(latest cvs changes: " // trim(conf%latest_cvs) // ")", 70)
     message(5) = ""
     call write_info(5)

     ! Let us print where we are running
     call loct_sysname(sys_name)
     write(message(1), '(a)') str_center("The octopus is swimming in " // trim(sys_name), 70)
     message(2) = ""
     call write_info(2)

#if defined(HAVE_MPI)
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

     call print_date("Calculation started on ")

     if(no_syslabels > 1) then
        message(1) = 'Info: Multi-Subsystem Mode'
        message(2) = 'Info: Running '//current_label
        call write_info(2, stress = .true.)
     endif

     ! create temporary dir (we will need it)
     tmpdir = 'tmp/'
     call io_mkdir(tmpdir)

     ! create debug directory if in debugging mode
     if(conf%verbose>=VERBOSE_DEBUG) call io_mkdir('debug')

     ! now we really start
     call run_init()
     call run()
     call run_end()

     call print_date("Calculation ended on ")

#if defined(HAVE_MPI)
     ! wait for all processors to finish
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

  enddo subsystems


  call io_end()
  call syslabels_end()
  call parser_end()
  call global_end()

end program octopus
