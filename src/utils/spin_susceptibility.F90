!! Copyright (C) 2008 X. Andrade
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

program spin_susceptibility
  use batch_oct_m
  use command_line_oct_m
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use spectrum_oct_m
  use simul_box_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  integer :: in_file, out_file, ref_file, ii, jj, kk, idir, jdir, ierr
  integer :: time_steps, time_steps_ref, energy_steps, istart, iend, ntiter
  FLOAT   :: dt, dt_ref, tt, ww
  FLOAT, allocatable :: ftreal(:,:), ftimag(:,:)
  FLOAT, allocatable :: magnetization(:,:)
  type(spectrum_t)  :: spectrum
  type(block_t)     :: blk
  type(batch_t)     :: vecpotb, ftrealb, ftimagb
  character(len=120) :: header
  FLOAT :: delta_strength
  character(len=MAX_PATH_LEN) :: ref_filename

  ! Initialize stuff
  call global_init(is_serial = .true.)

  call getopt_init(ierr)
  if(ierr == 0) call getopt_dielectric_function()
  call getopt_end()

  call messages_init()

  call io_init()

  call unit_system_init()

  call spectrum_init(spectrum)

  call parse_variable('TDDeltaStrength', M_ZERO, delta_strength )

  in_file = io_open('td.general/transverse_magnetization', action='read', status='old', die=.false.)
  if(in_file < 0) then 
    message(1) = "Cannot open file '"//trim(io_workpath('td.general/transverse_magnetization'))//"'"
    call messages_fatal(1)
  end if
  call io_skip_header(in_file)
  call spectrum_count_time_steps(in_file, time_steps, dt)

  time_steps = time_steps + 1

  SAFE_ALLOCATE(magnetization(1:time_steps, 1:4))
  
  call io_skip_header(in_file)
  
  do ii = 1, time_steps
    read(in_file, *) jj, tt, magnetization(ii,1), magnetization(ii,2), magnetization(ii,3), magnetization(ii,4)
  end do

  call io_close(in_file)

  write(message(1), '(a, i7, a)') "Info: Read ", time_steps, " steps from file '"// &
    trim(io_workpath('td.general/magnetization'))//"'"
  call messages_info(1)


  ! Find out the iteration numbers corresponding to the time limits.
  call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

  istart = max(1, istart)

  energy_steps = int(spectrum%max_energy / spectrum%energy_step)

  SAFE_ALLOCATE(ftreal(0:energy_steps, 4))
  SAFE_ALLOCATE(ftimag(0:energy_steps, 4))

  call batch_init(vecpotb, 4)
  call batch_init(ftrealb, 4)
  call batch_init(ftimagb, 4)

  call batch_add_state(vecpotb, magnetization(:, 1))
  call batch_add_state(vecpotb, magnetization(:, 2))
  call batch_add_state(vecpotb, magnetization(:, 3))
  call batch_add_state(vecpotb, magnetization(:, 4))
  call batch_add_state(ftrealb, ftreal(0:energy_steps,1))
  call batch_add_state(ftimagb, ftimag(0:energy_steps,1))
  call batch_add_state(ftrealb, ftreal(0:energy_steps,2))
  call batch_add_state(ftimagb, ftimag(0:energy_steps,2))
  call batch_add_state(ftrealb, ftreal(0:energy_steps,3))
  call batch_add_state(ftimagb, ftimag(0:energy_steps,3))
  call batch_add_state(ftrealb, ftreal(0:energy_steps,4))
  call batch_add_state(ftimagb, ftimag(0:energy_steps,4))


  call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart, iend, spectrum%start_time, dt, vecpotb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, vecpotb, 1, energy_steps + 1, spectrum%energy_step, ftrealb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, vecpotb, 1, energy_steps + 1, spectrum%energy_step, ftimagb)


  call batch_end(vecpotb)
  call batch_end(ftrealb)
  call batch_end(ftimagb)

  do kk = 0, energy_steps
    ftreal(kk, :) = ftreal(kk, :)/delta_strength
    ftimag(kk, :) = ftimag(kk, :)/delta_strength
  end do

  write(header, '(5a15)') '#        energy', 'Re[\chi_{+-}]', 'Im[\chi_{+-}]', 'Re[\chi_{-+}]', 'Im[\chi_{-+}]'

  out_file = io_open('td.general/spin_susceptibility', action='write')
  write(out_file,'(a)') trim(header)
  do kk = 0, energy_steps
    ww = kk*spectrum%energy_step
    write(out_file, '(5e15.6)') ww,                                   &
         (ftreal(kk, 1)-ftimag(kk, 2)), (ftimag(kk, 1)+ftreal(kk,2)), &
         (ftreal(kk, 3)-ftimag(kk, 4)), (ftimag(kk, 3)+ftreal(kk,4))
  end do
  call io_close(out_file)
 
  SAFE_DEALLOCATE_A(magnetization)
  SAFE_DEALLOCATE_A(ftreal)
  SAFE_DEALLOCATE_A(ftimag)
    
  call io_end()
  call messages_end()
  call global_end()

end program spin_susceptibility

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
