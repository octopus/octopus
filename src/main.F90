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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

program octopus
  use string_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use run_prog_m
  use io_m
  use profiling_m
  use varinfo_m
  use mpi_m

  implicit none

  integer :: ns
  integer(POINTER_SIZE) :: blk
  character(len=256) :: sys_name

  call global_init()
  call parser_init()

  !%Variable DevelVersion
  !%Type logical
  !%Default no
  !%Section Generalities::Debug
  !%Description
  !% If true, allows the use of certains parts of the code that are
  !% still under development. This should not be used for production runs.
  !%End
  call loct_parse_logical('DevelVersion', .false., conf%devel_version)

  !%Variable DebugLevel
  !%Type integer
  !%Default 1
  !%Section Generalities::Debug
  !%Description
  !% This variable decides wether or not to enter debug-mode. In debugging mode,
  !% the program prints to standard error when it enters and exits the subroutines,
  !% what is the memory it is using (only, for the moment being, in Linux systems),
  !% and some other information. Useful for developers, and mandatory if you want
  !% to send a bug report to the developers and being considered.
  !% You have two options: (i) setting it to zero -- or less than zero, in which
  !% case you do not run in debugging mode (this is the default), or (ii) setting
  !% it to a positive number. In this case the entries and exits to nested subroutines
  !% are only printed down to the level that is given in this variable.
  !%End
  call loct_parse_int('DebugLevel', 0, conf%debug_level)
  if(conf%debug_level>0) then
    in_debug_mode = .true.
  else
    in_debug_mode = .false.
  end if

  ! need to find out calc_mode already here since some of the variables here (e.g.
  ! periodic dimensions) can be different for the subsystems
  !%Variable CalculationMode
  !%Type integer
  !%Default gs
  !%Section Generalities
  !%Description
  !% Decides what kind of calculation is to be performed
  !%Option gs 1
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
    call datasets_init(calc_mode, blk)
  else
    call loct_parse_int('CalculationMode', M_GS, calc_mode)
    if(.not.varinfo_valid_option('CalculationMode', calc_mode)) call input_error('CalculationMode')
    call datasets_init(calc_mode)
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

  ! loop over all datasets
  datasets: do ns = 1, no_datasets

    ! set system label
    current_dataset = dataset_run_order(ns)
    current_label = trim(dataset_label(current_dataset))
    calc_mode = dataset_runmode(current_dataset)

    ! datasets have to be available before calling the _init() functions below
    call io_init()
    ! now we declare octopus as running
    call io_switch_status('running')

    call profiling_init()

    call profiling_in(C_PROFILING_COMPLETE_DATASET)

    ! Let us print our logo
    if(mpi_grp_is_root(mpi_world)) then
      call io_dump_file(stdout, trim(trim(conf%share) // '/logo'))
    end if

    ! Let us print the version
    message(1) = ""
    message(2) = str_center("Running octopus, version " // trim(conf%version), 70)
    message(3) = str_center("build time - " // trim(conf%build_time) , 70)
    message(4) = str_center("svn revision: " // trim(conf%latest_svn) , 70)
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
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
#endif

    call print_date("Calculation started on ")

    if(no_datasets > 1) then
      message(1) = 'Info: Multi-Dataset Mode'
      message(2) = 'Info: Running dataset "'//trim(current_label)//'"'
      call write_info(2, stress = .true.)
    end if

    ! now we really start
    call run_init()
    call run()
    call run_end()

#if defined(HAVE_MPI)
    ! wait for all processors to finish
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
#endif

    ! run finished successfully
    call io_switch_status('finished')
    call io_end()

    call profiling_out(C_PROFILING_COMPLETE_DATASET)
    call profiling_output()

    call print_date("Calculation ended on ")

  end do datasets

  call datasets_end()
  call parser_end()
  call global_end()

end program octopus
