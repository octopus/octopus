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
  use parser_m
  use messages_m
  use mpi_m
  use profiling_m
  use run_m
  use string_m
  use varinfo_m

  implicit none

  integer :: ns, inp_calc_mode
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
  call parse_logical('DevelVersion', .false., conf%devel_version)

  !%Variable DebugLevel
  !%Type integer
  !%Default 0
  !%Section Execution::Debug
  !%Description
  !% This variable decides whether or not to enter debug mode.
  !% If it is greater than 0, different amounts of additional information
  !% are written to standard output and additional assertion checks are performed.
  !%Option 0
  !% (default) <tt>Octopus</tt> does not enter debug mode.
  !%Option 1
  !% Moderate amount of debug output; assertion checks enabled.
  !%Option 2
  !% The code prints a stack trace as it enters end exits subroutines.
  !% This is useful for developers and you should include this output when
  !% submitting a bug report.
  !%Option 100
  !% The debug output is additionally written to files in the debug
  !% directory. For each node (when running in parallel) there is a file called
  !% <tt>debug_trace.&lt;rank&gt;</tt>. Writing these files slows down the code by a huge factor and
  !% it is usually only necessary for parallel runs. In the serial case all
  !% the information can be obtained from standard out.
  !%End
  call parse_integer('DebugLevel', 0, conf%debug_level)
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
  !% If true, <tt>Octopus</tt> will print as part of the screen output
  !% information about the memory the code is using. The quantity
  !% reported is an approximation to the size of the heap and
  !% generally it is a lower bound to the actual memory <tt>Octopus</tt> is
  !% using. By default this variable is set to false.
  !%End
  call parse_logical('ReportMemory', .false., conf%report_memory)

  ! need to find out calc_mode already here since some of the variables here (e.g.
  ! periodic dimensions) can be different for the subsystems

  !%Variable CalculationMode
  !%Type integer
  !%Default gs
  !%Section Calculation Modes
  !%Description
  !% Decides what kind of calculation is to be performed.
  !%Option gs 01
  !% Calculation of the ground state.
  !%Option unocc 02
  !% Calculation of unoccupied/virtual KS states.
  !%Option td 03
  !% Time-dependent calculation.
  !%Option go 05
  !% Optimization of the geometry.
  !%Option opt_control 07
  !% Optimal control.
  !%Option em_resp 08
  !% Calculation of the electromagnetic response: electric
  !% polarizabilities and hyperpolarizabilities and magnetic
  !% susceptibilities.
  !% For periodic system, currently available only in development version.
  !%Option casida 09
  !% Excitations via Casida linear-response TDDFT; for finite systems only.
  !%Option td_transport 10
  !% Time-dependent quantum transport.
  !%Option vdw 11
  !% Calculate van der Waals coefficients.
  !%Option vib_modes 12
  !% Calculation of the vibrational modes.
  !%Option raman 13
  !% Calculation of Raman response properties.
  !% Currently available only in development version.
  !%Option one_shot 14
  !% Use the self-consistent wavefunctions in the <tt>restart</tt> directory to
  !% evaluate the total energy using a different XC functional.
  !% This is effectively a first-order perturbative calculation of the total energy, 
  !% the perturbation being the difference between the two XC potentials used.
  !%Option kdotp 15
  !% Calculation of effective masses by <i>k.p</i> perturbation theory.
  !% Currently available only in development version.
  !%Option gcm 16
  !% Generator-Coordinates calculation (experimental).
  !%Option memory 17
  !% It tells you the approximate amount of memory <tt>Octopus</tt> will need to run.
  !%Option invert_ks 18
  !% Run mode used to invert the Kohn-Sham equations.
  !%Option recipe 99
  !% Prints out a tasty recipe.
  !%End
  if(parse_block('CalculationMode', blk) == 0) then
    call datasets_init(inp_calc_mode, blk)
  else
    call parse_integer('CalculationMode', CM_GS, inp_calc_mode)
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
    message(2) = str_center("Running octopus", 70)
    message(3) = ""
    call write_info(3)

    message(1) = &
      "Version                : " // trim(conf%version)
    message(2) = &
      "Revision               : "// trim(conf%latest_svn)
    message(3) = &
      "Build time             : "// trim(conf%build_time)
    call write_info(3)

    write(message(1), '(a, i1)') &
      'Configuration options  : max-dim=', MAX_DIM

#ifdef USE_OMP
    message(1) = trim(message(1))//' openmp'
#endif
#ifdef HAVE_MPI
    message(1) = trim(message(1))//' mpi'
#endif

    message(2) = &
      'Optional libraries     :'
#ifdef HAVE_MPI2
    message(2) = trim(message(2))//' mpi2'
#endif
#ifdef HAVE_NETCDF
    message(2) = trim(message(2))//' netcdf'
#endif
#ifdef HAVE_METIS
    message(2) = trim(message(2))//' metis'
#endif
#ifdef HAVE_GDLIB
    message(2) = trim(message(2))//' gdlib'
#endif
#ifdef HAVE_LIBNBC
    message(2) = trim(message(2))//' libnbc'
#endif
#ifdef HAVE_PAPI
    message(2) = trim(message(2))//' papi'
#endif
#ifdef HAVE_SPARSKIT
    message(2) = trim(message(2))//' sparskit'
#endif

    message(3) = &
      'Architecture           : '// TOSTRING(OCT_ARCH)
    call write_info(3)

    message(1) = &
      "C compiler             : "//trim(conf%cc)
    message(2) = &
      "C compiler flags       : "//trim(conf%cflags)
    message(3) = &
      "Fortran compiler       : "//trim(conf%fc)
    message(4) = &
      "Fortran compiler flags : "//trim(conf%fcflags)
    call write_info(4)

    message(1) = ""
    call write_info(1)

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
    !% <tt>Octopus</tt> can run in 1, 2 or 3 dimensions, depending on the value of this
    !% variable. Note that not all input variables may be available in all cases.
    !%End
    call parse_integer(datasets_check('Dimensions'), 3, calc_dim)
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
