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
  use command_line_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use run_m
  use string_m
  use varinfo_m

  implicit none

  integer :: ns, inp_calc_mode, ierr
  type(block_t) :: blk
  character(len=256) :: sys_name

  call getopt_init(ierr)
  if(ierr .eq. 0) call getopt_octopus()
  call getopt_end()

  call global_init()
  call parser_init()
  call messages_init()

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
  !% Time-dependent calculation (experimental for periodic systems).
  !%Option go 05
  !% Optimization of the geometry.
  !%Option opt_control 07
  !% Optimal control.
  !%Option em_resp 08
  !% Calculation of the electromagnetic response: electric
  !% polarizabilities and hyperpolarizabilities and magnetic
  !% susceptibilities (experimental for periodic systems).
  !%Option casida 09
  !% Excitations via Casida linear-response TDDFT; for finite systems only.
  !%Option td_transport 10
  !% Time-dependent quantum transport (experimental).
  !%Option vdw 11
  !% Calculate van der Waals coefficients.
  !%Option vib_modes 12
  !% Calculation of the vibrational modes.
  !%Option one_shot 14
  !% Use the self-consistent wavefunctions in the <tt>restart</tt> directory to
  !% evaluate the total energy using a different XC functional.
  !% This is effectively a first-order perturbative calculation of the total energy, 
  !% the perturbation being the difference between the two XC potentials used.
  !%Option kdotp 15
  !% Calculation of effective masses by <i>k.p</i> perturbation theory (experimental).
  !%Option gcm 16
  !% Generator-Coordinates Method calculation (experimental).
  !% Ref. K. Capelle, <i>J. Chem. Phys.</i> <b>119</b>, 1285 (2003).
  !%Option dummy 17
  !% This calculation mode does nothing. Useful for debugging, testing and benchmarking.  
  !%Option invert_ks 18
  !% Invert the Kohn-Sham equations (experimental).
  !%Option recipe 99
  !% Prints out a tasty recipe.
  !%
  !% May also be used as a block for multi-dataset mode. The first line is a list of calculation modes,
  !% the second is labels (optional), and the third is the order for the runs (optional). Example:
  !%
  !% <pre>%CalculationMode
  !%   gs     | unocc  | td
  !%   "run1" | "run2" | "run3"
  !%   1      | 2      | 3
  !% %</pre>
  !%End
  if(parse_block('CalculationMode', blk) == 0) then
    call datasets_init(inp_calc_mode, blk)
  else
    call parse_integer('CalculationMode', CM_GS, inp_calc_mode)
    if(.not.varinfo_valid_option('CalculationMode', inp_calc_mode)) call input_error('CalculationMode')
    call datasets_init(inp_calc_mode)
  end if

  ! Now we can initialize the I/O
  call io_init()

  ! loop over all datasets
  datasets: do ns = 1, no_datasets

    ! set system label
    current_dataset = dataset_run_order(ns)
    current_label = trim(dataset_label(current_dataset))
    call calc_mode_init()

    ! datasets have to be available before calling the _init() functions below
    call io_init_datasets()

    ! now we declare octopus as running
    call io_switch_status('running')

    call profiling_init()

    ! Let us print our logo
    if(mpi_grp_is_root(mpi_world)) then
      call io_dump_file(stdout, trim(trim(conf%share) // '/logo'))
    end if

    ! Let us print the version
    message(1) = ""
    message(2) = str_center("Running octopus", 70)
    message(3) = ""
    call messages_info(3)

    message(1) = &
      "Version                : " // trim(conf%version)
    message(2) = &
      "Revision               : "// trim(conf%latest_svn)
    message(3) = &
      "Build time             : "// trim(conf%build_time)
    call messages_info(3)

    write(message(1), '(a, i1)') &
      'Configuration options  : max-dim=', MAX_DIM

#ifdef HAVE_OPENMP
    message(1) = trim(message(1))//' openmp'
#endif
#ifdef HAVE_MPI
    message(1) = trim(message(1))//' mpi'
#endif
#ifdef HAVE_OPENCL
    message(1) = trim(message(1))//' opencl'
#endif
#ifdef HAVE_M128D
    message(1) = trim(message(1))//' sse2'
#endif
#ifdef HAVE_M256D
    message(1) = trim(message(1))//' avx'
#endif
#ifdef HAVE_BLUE_GENE
    message(1) = trim(message(1))//' bluegene'
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
#ifdef HAVE_PAPI
    message(2) = trim(message(2))//' papi'
#endif
#ifdef HAVE_SPARSKIT
    message(2) = trim(message(2))//' sparskit'
#endif

    message(3) = &
      'Architecture           : '// TOSTRING(OCT_ARCH)
    call messages_info(3)

    message(1) = &
      "C compiler             : "//trim(conf%cc)
    message(2) = &
      "C compiler flags       : "//trim(conf%cflags)
    message(3) = &
      "Fortran compiler       : "//trim(conf%fc)
    message(4) = &
      "Fortran compiler flags : "//trim(conf%fcflags)
    call messages_info(4)

    message(1) = ""
    call messages_info(1)

    ! Let us print where we are running
    call loct_sysname(sys_name)
    write(message(1), '(a)') str_center("The octopus is swimming in " // trim(sys_name), 70)
    message(2) = ""
    call messages_info(2)

#if defined(HAVE_MPI)
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

    call print_date("Calculation started on ")

    if(no_datasets > 1) then
      message(1) = 'Info: Multi-Dataset Mode'
      message(2) = 'Info: Running dataset "'//trim(current_label)//'"'
      call messages_info(2, stress = .true.)
    end if

    ! now we really start
    call run_init(dataset_runmode(current_dataset))
    call run()
    call run_end()
    
#if defined(HAVE_MPI)
    ! wait for all processors to finish
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

    ! run finished successfully
    call io_switch_status('finished')
    call io_end()

    call profiling_end()

    call calc_mode_end()
    
    call print_date("Calculation ended on ")
    call print_walltime()
  end do datasets

  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()

contains

  subroutine print_walltime()
    integer :: days, hours, min, sec, usec

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)
    
    days  = sec / 86400
    hours = (sec / 3600) - (days * 24)
    min   = (sec / 60) - (days * 1440) - (hours * 60)
    sec   = modulo(sec, 60)

    message(2) = ''
    if(days  > 0) write(message(2), '(i3,a)') days, ' days,'
    if(hours > 0.or.message(2).ne.'') &
      write(message(2), '(a,1x,i2.2,a)') trim(message(2)), hours, 'h'
    if(min   > 0.or.message(1).ne.'') &
      write(message(2), '(a,1x,i2.2,a)') trim(message(2)), min, 'm'
    write(message(2), '(a,1x,i2.2,a,i3,a)') trim(message(2)), sec, '.', usec/1000, 's'
    message(1) = str_center('Walltime: ' // trim(message(2)), 70)
    call messages_info(1)

  end subroutine print_walltime

end program octopus

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
