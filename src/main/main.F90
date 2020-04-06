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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

program octopus
  use calc_mode_par_oct_m
  use command_line_oct_m
  use io_oct_m
  use global_oct_m
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use run_oct_m
  use string_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use walltimer_oct_m

  implicit none

  character(len=256) :: config_str
  integer :: inp_calc_mode, ierr
  type(block_t) :: blk

  call getopt_init(ierr)
  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(trim(config_str))
  call getopt_end()

  call global_init()

  call parser_init()
  
  call messages_init(global_namespace)

  call walltimer_init(global_namespace)
  
  !%Variable ReportMemory
  !%Type logical
  !%Default no
  !%Section Execution::Debug
  !%Description
  !% If true, after each SCF iteration <tt>Octopus</tt> will print
  !% information about the memory the code is using. The quantity
  !% reported is an approximation to the size of the heap and
  !% generally it is a lower bound to the actual memory <tt>Octopus</tt> is
  !% using.
  !%End
  call parse_variable(global_namespace, 'ReportMemory', .false., conf%report_memory)

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
  !% Calculation of unoccupied/virtual KS states. Can also be used for a non-self-consistent
  !% calculation of states at arbitrary k-points, if <tt>density.obf</tt> from <tt>gs</tt>
  !% is provided in the <tt>restart/gs</tt> directory.
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
  !%Option vdw 11
  !% Calculate van der Waals coefficients.
  !%Option vib_modes 12
  !% Calculation of the vibrational modes.
  !%Option one_shot 14
  !% Obsolete. Use <tt>gs</tt> with <tt>MaximumIter = 0</tt> instead.
  !%Option kdotp 15
  !% Calculation of effective masses by <math>\vec{k} \cdot \vec{p}</math> perturbation theory (experimental).
  !%Option dummy 17
  !% This calculation mode does nothing. Useful for debugging, testing and benchmarking.  
  !%Option invert_ks 18
  !% Invert the Kohn-Sham equations (experimental).
  !%Option test 19
  !%Option recipe 99
  !% Prints out a tasty recipe.
  !%End
  if(parse_block(global_namespace, 'CalculationMode', blk) == 0) then
    call messages_write('The datasets mode has been deprecated,', new_line = .true.)
    call messages_write('please use several Octopus runs.')
    call messages_fatal()
  end if

  call parse_variable(global_namespace, 'CalculationMode', CM_GS, inp_calc_mode)
  if(.not.varinfo_valid_option('CalculationMode', inp_calc_mode)) call messages_input_error('CalculationMode')

  ! Now we can initialize the I/O
  call io_init(global_namespace)

  call calc_mode_par_init()
  
  ! now we declare octopus as running
  call messages_switch_status('running')
  
  call profiling_init(global_namespace)
  
  call print_header()

#if !defined(HAVE_LIBXC4)
  call messages_write('You have compiled Octopus with version 3 of Libxc.', new_line = .true.)
  call messages_write('Support for this version of Libxc has been deprecated and', new_line = .true.)
  call messages_write('will be removed in one of the next major releases of Octopus.', new_line = .true.)
  call messages_warning()
#endif
  
  ! now we really start
  call run(global_namespace, inp_calc_mode)
  
#if defined(HAVE_MPI)
  ! wait for all processors to finish
  call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
  
  ! run finished successfully
  call messages_switch_status('finished')
  call io_end()
  
  call profiling_end(global_namespace)
  
  call calc_mode_par_end()

  call walltimer_end()
  
  call print_date("Calculation ended on ")
  call print_walltime()

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
    if(hours > 0.or.message(2) /= '') &
      write(message(2), '(a,1x,i2.2,a)') trim(message(2)), hours, 'h'
    if(min   > 0.or.message(1) /= '') &
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
