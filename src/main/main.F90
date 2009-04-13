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

program octopus
  use calc_mode_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use loct_parser_m
  use messages_m
  use mpi_m
  use profiling_m
  use run_m
  use string_m
  use varinfo_m

  implicit none

  integer :: ns, inp_calc_mode, lenFCFLAGS
  type(block_t) :: blk
  character(len=256) :: sys_name

  call global_init()
  call parser_init()

  !%Variable DevelVersion
  !%Type logical
  !%Default no
  !%Section Execution::Debug
  !%Description
  !% If true, allows the use of certain parts of the code that are
  !% still under development. This should not be used for production runs.
  !%End
  call loct_parse_logical('DevelVersion', .false., conf%devel_version)

  !%Variable DebugLevel
  !%Type integer
  !%Default 0
  !%Section Execution::Debug
  !%Description
  !% This variable decides whether or not to enter debug mode.
  !% If its value is 0 (the default) Octopus does not enter debug mode.
  !% If it is greater than 0 different amounts of additional information
  !% are written to standard output and additional assertion checks are performed.
  !% Level 1: moderate amount of debug output but assertion checks enabled.
  !% Level 2: the code prints a stack trace as it enters end exits subroutines.
  !% This is useful for developers and you should include this output when
  !% submitting a bug report.
  !% Level 100: the debug output is additionally written to files in the debug
  !% directory. For each node (when running in parallel) there is a file called
  !% debug_trace.<rank>. Writing these files slows down the code by a huge factor and
  !% it is usually only necessary for parallel runs. In the serial case all
  !% the information can be obtained from standard out.
  !%End
  call loct_parse_int('DebugLevel', 0, conf%debug_level)
  if(conf%debug_level>0) then
    in_debug_mode = .true.
  else
    in_debug_mode = .false.
  end if

  ! Now we can initialize the io
  call io_init()

  !%Variable ReportMemory
  !%Type logical
  !%Default no
  !%Section Execution::Debug
  !%Description
  !% If true, octopus will print as part of the screen output
  !% information about the memory the code is using. The quantity
  !% reported is an approximation to the size of the heap and
  !% generally it is a lower bound to the actual memory octopus is
  !% using. By default this variable is set to false.
  !%End
  call loct_parse_logical('ReportMemory', .false., conf%report_memory)

  ! need to find out calc_mode already here since some of the variables here (e.g.
  ! periodic dimensions) can be different for the subsystems

  !%Variable CalculationMode
  !%Type integer
  !%Default gs
  !%Section Calculation Modes
  !%Description
  !% Decides what kind of calculation is to be performed.
  !%Option gs 01
  !% Calculation of the ground state
  !%Option unocc 02
  !% Calculation of unoccupied/virtual KS states
  !%Option td 03
  !% Time-dependent calculation
  !%Option go 05
  !% Optimization of the geometry
  !%Option opt_control 07
  !% Optimal control.
  !%Option em_resp 08
  !% Calculation of the electromagnetic response: electric
  !% polarizabilities and hyperpolarizabilities and magnetic
  !% susceptibilities.
  !% For periodic system, currently available only in development version.
  !%Option casida 09
  !% Excitations via Casida linear-response TDDFT
  !%Option td_transport 10
  !% Time-dependent quantum transport
  !%Option vdw 11
  !% Calculate van der Waals coefficients
  !%Option vib_modes 12
  !% Calculation of the vibrational modes.
  !%Option raman 13
  !% Calculation of Raman response properties.
  !% Currently available only in development version.
  !%Option one_shot 14
  !% Use the self-consistent wave-functions in the restart directory to
  !% evaluate the total energy using a different xc functional.
  !% This is effectively a first-order perturbative calculation of the total energy, 
  !% the perturbation being the difference between the two xc used.
  !%Option kdotp 15
  !% Calculation of effective masses by k.p perturbation theory.
  !% Currently available only in development version.
  !%Option gcm 16
  !% Generator-Coordinates calculation (experimental).
  !%Option memory 17
  !% It tells you the approximate amount of memory octopus will need to run.
  !%Option recipe 99
  !% Prints out a tasty recipe
  !%End
  if(loct_parse_block('CalculationMode', blk) == 0) then
    call datasets_init(inp_calc_mode, blk)
  else
    call loct_parse_int('CalculationMode', CM_GS, inp_calc_mode)
    if(.not.varinfo_valid_option('CalculationMode', inp_calc_mode)) call input_error('CalculationMode')
    call datasets_init(inp_calc_mode)
  end if

  call calc_mode_set(inp_calc_mode)

  ! loop over all datasets
  datasets: do ns = 1, no_datasets

    ! set system label
    current_dataset = dataset_run_order(ns)
    current_label = trim(dataset_label(current_dataset))
    call calc_mode_set(dataset_runmode(current_dataset))

    ! datasets have to be available before calling the _init() functions below
    call io_init_datasets()

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

    lenFCFLAGS = len(trim(conf%fcflags))
    if(lenFCFLAGS + 17 .le. 80) then
       message(7) = str_center("Compiler flags: "//trim(conf%fcflags), max(70, len(trim(conf%fcflags))+17))
    else
       message(7) = "Compiler flags: "//trim(conf%fcflags)
       ! allow full string to be written out
    endif

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

    !%Variable Dimensions
    !%Type integer
    !%Section System
    !%Default 3
    !%Description
    !% octopus can run in 1, 2 or 3 dimensions, depending on the value of this
    !% variable. Note that not all input variables may be available in all cases.
    !%End
    call loct_parse_int(datasets_check('Dimensions'), 3, calc_dim)
    if( calc_dim > MAX_DIM .or. calc_dim < 1) call input_error('Dimensions')

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
    call profiling_end()

    call print_date("Calculation ended on ")

  end do datasets

  call datasets_end()
  call parser_end()
  call global_end()

end program octopus

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
