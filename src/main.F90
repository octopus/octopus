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
  use profiling_mod
  use varinfo
  use mpi_mod

  implicit none

  integer :: ns
  integer(POINTER_SIZE) :: blk
  character(len=256) :: sys_name
#if defined(HAVE_MPI)
  integer :: ierr
#endif

  call global_init()
  call parser_init()

  ! need to find out calc_mode already here since some of the variables here (e.g.
  ! periodic dimensions) can be different for the subsystems
  !%Variable CalculationMode
  !%Type integer
  !%Default gs
  !%Section Generalities
  !%Description
  !% Decides what kind of calculation is to be performed
  !%Option gs 01
  !% Calculation of the ground state
  !%Option unocc 02
  !% Calculation of unoccupied/virtual KS states
  !%Option td 03
  !% Time-dependent calculation
  !%Option pol 04
  !% Calculation of the static polarizability
  !%Option geom 05
  !% Optimization of the geometry
  !%Option phonons 06
  !% Calculation of the vibrational modes
  !%Option opt_control 07
  !% Optimal control.
  !%Option pol_lr 08
  !% Linear-response calculation of the polarizability
  !%Option casida 09
  !% Excitations via linear-response TDDFT
  !%Option wave_matching 10
  !% Wave-matching a la Heiko
  !%Option bo 98
  !% Born-Oppenheimer-like Molecular Dynamics
  !%Option recipe 99
  !% Prints out a tasty recipe
  !%End
  if(loct_parse_block('CalculationMode', blk) == 0) then
    call syslabels_init(calc_mode, blk)
  else
    call loct_parse_int('CalculationMode', M_GS, calc_mode)
    if(.not.varinfo_valid_option('CalculationMode', calc_mode)) call input_error('CalculationMode')
    call syslabels_init(calc_mode)
  end if

  !%Variable Dimensions
  !%Type integer
  !%Section Generalities
  !%Default 3
  !%Description
  !% octopus can run in 1, 2 or 3 dimensions, depending on the value of this
  !% variable. Note that not all input variables may be available in all cases.
  !%Option 1
  !% The system is 1-dimensional
  !%Option 2
  !% The system is 2-dimensional
  !%Option 3
  !% The system is 3-dimensional
  !%End
  call loct_parse_int(check_inp('Dimensions'), 3, calc_dim)
  if( calc_dim > 3 .or. calc_dim < 1) call input_error('Dimensions')

  ! loop over all subsystems
  subsystems: do ns = 1, no_subsystems

    ! set system label
    current_subsystem = subsys_run_order(ns)
    current_label = trim(subsys_label(current_subsystem))
    calc_mode = subsys_runmode(current_subsystem)

    ! syslabels have to be available before calling the _init() functions below
    call io_init()
    call profiling_init()

    call profiling_in(C_PROFILING_COMPLETE_SUBSYS)

    ! Let us print our logo
    if(mpi_grp_is_root(mpi_world)) then
      call io_dump_file(stdout, trim(trim(conf%share) // '/logo'))
    end if

    ! Let us print the version
    message(1) = ""
    message(2) = str_center("Running octopus, version " // trim(conf%version), 70)
    message(3) = str_center("(build time - " // trim(conf%build_time) // ")", 70)
    message(4) = str_center("(latest cvs changes: " // trim(conf%latest_cvs) // ")", 70)
    message(5) = ""
    message(6) = str_center("Compiler: "//trim(conf%compiler), 70)
    message(7) = str_center("Compiler flags: "//trim(conf%fcflags), max(70, len(trim(conf%fcflags))+17))
    message(8) = ""
    call write_info(8)

    ! Let us print where we are running
    call loct_sysname(sys_name)
    write(message(1), '(a)') str_center("The octopus is swimming in " // trim(sys_name), 70)
    message(2) = ""
    call write_info(2)

#if defined(HAVE_MPI)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

    call print_date("Calculation started on ")

    if(no_subsystems > 1) then
      message(1) = 'Info: Multi-Subsystem Mode'
      message(2) = 'Info: Running '//current_label
      call write_info(2, stress = .true.)
    end if

    ! now we really start
    call run_init()
    call run()
    call run_end()

#if defined(HAVE_MPI)
    ! wait for all processors to finish
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

    call io_end()

    call profiling_out(C_PROFILING_COMPLETE_SUBSYS)
    call profiling_output()

    call print_date("Calculation ended on ")

  end do subsystems

  call syslabels_end()
  call parser_end()
  call global_end()

end program octopus
