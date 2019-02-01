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
  use global_oct_m
  use calc_mode_par_oct_m
  use command_line_oct_m
  use io_oct_m
  use loct_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use run_oct_m
  use string_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use messages_oct_m

  implicit none

  character(len=256) :: config_str
  integer :: inp_calc_mode, ierr
  type(block_t) :: blk

  call getopt_init(ierr)
  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(trim(config_str))
  call getopt_end()

  call global_init()
  call messages_init()

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
  call parse_variable('ReportMemory', .false., conf%report_memory)

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
  !%Option td 03
  !% Time-dependent calculation (experimental for periodic systems).
  !%Option test 19
  !%End
  if(parse_block('CalculationMode', blk) == 0) then
    call messages_write('The datasets mode has been deprecated,', new_line = .true.)
    call messages_write('please use several Octopus runs.')
    call messages_fatal()
  end if

  call parse_variable('CalculationMode', CM_GS, inp_calc_mode)
  if(.not.varinfo_valid_option('CalculationMode', inp_calc_mode)) call messages_input_error('CalculationMode')

  ! Now we can initialize the I/O
  call io_init()

  call calc_mode_par_init()
  
  ! now we declare octopus as running
  call messages_switch_status('running')
  
  call profiling_init()
  
  call print_header()
  
  ! now we really start
  call run(inp_calc_mode)
  
#if defined(HAVE_MPI)
  ! wait for all processors to finish
  call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
  
  ! run finished successfully
  call messages_switch_status('finished')
  call io_end()
  
  call profiling_end()
  
  call calc_mode_par_end()
  
  call print_date("Calculation ended on ")
  call print_walltime()

  call messages_end()
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
