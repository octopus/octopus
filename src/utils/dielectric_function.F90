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

program dielectric_function
  use batch_oct_m
  use command_line_oct_m
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use lalg_adv_oct_m
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
  FLOAT   :: dt, dt_ref, tt, ww, n0
  FLOAT, allocatable :: vecpot(:, :), vecpot0(:), ftreal(:, :), ftimag(:, :)
  CMPLX, allocatable :: dielectric(:, :), chi(:, :), invdielectric(:, :), fullmat(:, :)
  FLOAT, allocatable :: vecpot_ref(:, :)
  type(spectrum_t)  :: spectrum
  type(block_t)     :: blk
  type(space_t)     :: space
  type(geometry_t)  :: geo
  type(simul_box_t) :: sb
  type(batch_t)     :: vecpotb, ftrealb, ftimagb
  character(len=120) :: header
  FLOAT :: start_time
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

  call space_init(space)
  call geometry_init(geo, space)
  call simul_box_init(sb, geo, space)
    
  SAFE_ALLOCATE(vecpot0(1:space%dim))

  if(parse_block('GaugeVectorField', blk) == 0) then
    
    do ii = 1, space%dim
      call parse_block_float(blk, 0, ii - 1, vecpot0(ii))
    end do
    
    call parse_block_end(blk)

  else
    
    if(in_file < 0) then 
      message(1) = "Cannot find the GaugeVectorField in the input file"
      call messages_fatal(1)
    end if

  end if

  message(1) = "This program assumes that the gauge field is in the 'x'"
  message(2) = "direction, and that the 'y' and 'z' directions are equivalent."
  message(3) = "If this is not the case the dielectric function and the"
  message(4) = "susceptibility will be wrong."
  call messages_warning(4)

  start_time = spectrum%start_time
  call parse_variable('GaugeFieldDelay', start_time, spectrum%start_time )

  in_file = io_open('td.general/gauge_field', action='read', status='old', die=.false.)
  if(in_file < 0) then 
    message(1) = "Cannot open file '"//trim(io_workpath('td.general/gauge_field'))//"'"
    call messages_fatal(1)
  end if
  call io_skip_header(in_file)
  call spectrum_count_time_steps(in_file, time_steps, dt)

  if(parse_is_defined('TransientAbsorptionReference')) then
    !%Variable TransientAbsorptionReference
    !%Type string
    !%Default "."
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% In case of delayed kick, the calculation of the transient absorption requires 
    !% to substract a reference calculation, containing the gauge-field without the kick
    !% This reference must be computed using GaugeFieldPropagate=yes and to have
    !% TDOutput = gauge_field.
    !% This variables defined the directory in which the reference gauge_field field is,
    !% relative to the current folder
    !%End

    call parse_variable('TransientAbsorptionReference', '.', ref_filename)
    ref_file = io_open(trim(ref_filename)//'/gauge_field', action='read', status='old', die=.false.)
    if(ref_file < 0) then
      message(1) = "Cannot open reference file '"//trim(io_workpath(trim(ref_filename)//'/gauge_field'))//"'"
      call messages_fatal(1)
    end if
    call io_skip_header(ref_file)
    call spectrum_count_time_steps(ref_file, time_steps_ref, dt_ref)
    if(time_steps_ref < time_steps) then
      message(1) = "The reference calculation does not contain enought time steps"
      call messages_fatal(1)
    end if
 
    if(dt_ref /= dt) then
      message(1) = "The time step of the reference calculation is different from the current calculation"
      call messages_fatal(1)
    end if

  end if
  
  time_steps = time_steps + 1

  SAFE_ALLOCATE(vecpot(1:time_steps, space%dim*3))
  
  call io_skip_header(in_file)
  
  do ii = 1, time_steps
    read(in_file, *) jj, tt, (vecpot(ii, kk), kk = 1, space%dim*3)
  end do

  call io_close(in_file)

  !We remove the reference
  if(parse_is_defined('TransientAbsorptionReference')) then
    time_steps_ref = time_steps_ref + 1
    SAFE_ALLOCATE(vecpot_ref(1:time_steps_ref, space%dim*3))
    call io_skip_header(ref_file)
    do ii = 1, time_steps_ref
      read(ref_file, *) jj, tt, (vecpot_ref(ii, kk), kk = 1, space%dim*3)
    end do
    call io_close(ref_file)
    do ii = 1, time_steps
      do kk = 1, space%dim*3
        vecpot(ii, kk) = vecpot(ii, kk) - vecpot_ref(ii, kk)
      end do
    end do
  end if

  write(message(1), '(a, i7, a)') "Info: Read ", time_steps, " steps from file '"// &
    trim(io_workpath('td.general/gauge_field'))//"'"
  call messages_info(1)


  ! Find out the iteration numbers corresponding to the time limits.
  call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

  istart = max(1, istart)

  energy_steps = int(spectrum%max_energy / spectrum%energy_step)

  n0 = sqrt(sum(vecpot0(1:space%dim))**2)

  SAFE_ALLOCATE(ftreal(0:energy_steps, 1:space%dim))
  SAFE_ALLOCATE(ftimag(0:energy_steps, 1:space%dim))

  call batch_init(vecpotb, space%dim)
  call batch_init(ftrealb, space%dim)
  call batch_init(ftimagb, space%dim)

  do ii = 1, space%dim
    call batch_add_state(vecpotb, vecpot(:,  space%dim + ii))
    call batch_add_state(ftrealb, ftreal(0:,  ii))
    call batch_add_state(ftimagb, ftimag(0:,  ii))
  end do

  call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart, iend, spectrum%start_time, dt, vecpotb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, vecpotb, 1, energy_steps + 1, spectrum%energy_step, ftrealb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, vecpotb, 1, energy_steps + 1, spectrum%energy_step, ftimagb)


  call batch_end(vecpotb)
  call batch_end(ftrealb)
  call batch_end(ftimagb)

  SAFE_ALLOCATE(invdielectric(1:space%dim, 0:energy_steps))
  SAFE_ALLOCATE(dielectric(1:space%dim, 0:energy_steps))
  SAFE_ALLOCATE(chi(1:space%dim, 0:energy_steps))
  SAFE_ALLOCATE(fullmat(1:space%dim, 1:space%dim))

  do kk = 0, energy_steps
    ww = kk*spectrum%energy_step

    invdielectric(1:space%dim, kk) = (vecpot0(1:space%dim) + TOCMPLX(ftreal(kk, 1:space%dim), ftimag(kk, 1:space%dim)))/n0

    ! calculate the full inverse dielectric matrix
    do idir = 1, space%dim
      do jdir = 1, space%dim
        fullmat(idir, mod(jdir + idir - 2, space%dim) + 1) = invdielectric(jdir, kk)
      end do
    end do

    ! and invert it
    call lalg_sym_inverter('u', space%dim, fullmat)

    ! symmetrize the rest dielectric matrix
    do idir = 2, space%dim
      do jdir = 1, idir - 1
        fullmat(jdir, idir) = fullmat(idir, jdir)
      end do
    end do
    
    dielectric(1:space%dim, kk) = fullmat(1:space%dim, 1)

    chi(1:space%dim, kk) = (dielectric(1:space%dim, kk) - vecpot0(1:space%dim)/n0)*sb%rcell_volume/(M_FOUR*M_PI)
  end do

  SAFE_DEALLOCATE_A(fullmat)

  write(header, '(7a15)') '#        energy', 'Re x', 'Im x', 'Re y', 'Im y', 'Re z', 'Im z'

  out_file = io_open('td.general/inverse_dielectric_function', action='write')
  write(out_file,'(a)') trim(header)
  do kk = 0, energy_steps
    ww = kk*spectrum%energy_step
    write(out_file, '(7e15.6)') ww,                                         &
         real(invdielectric(1, kk), REAL_PRECISION), aimag(invdielectric(1, kk)), &
         real(invdielectric(2, kk), REAL_PRECISION), aimag(invdielectric(2, kk)), &
         real(invdielectric(3, kk), REAL_PRECISION), aimag(invdielectric(3, kk))
  end do
  call io_close(out_file)
 
  out_file = io_open('td.general/dielectric_function', action='write')
  write(out_file,'(a)') trim(header)
  do kk = 0, energy_steps
    ww = kk*spectrum%energy_step
    write(out_file, '(7e15.6)') ww,                                         &
         real(dielectric(1, kk), REAL_PRECISION), aimag(dielectric(1, kk)), &
         real(dielectric(2, kk), REAL_PRECISION), aimag(dielectric(2, kk)), &
         real(dielectric(3, kk), REAL_PRECISION), aimag(dielectric(3, kk))
  end do
  call io_close(out_file)

  out_file = io_open('td.general/chi', action='write')
  write(out_file,'(a)') trim(header)
  do kk = 0, energy_steps
    dielectric(1:3, kk) = (dielectric(1:3, kk) - M_ONE)/(CNST(4.0)*M_PI)
    ww = kk*spectrum%energy_step
    write(out_file, '(7e15.6)') ww, &
      real(chi(1, kk), REAL_PRECISION), aimag(chi(1, kk)), &
      real(chi(2, kk), REAL_PRECISION), aimag(chi(2, kk)), &
      real(chi(3, kk), REAL_PRECISION), aimag(chi(3, kk))
  end do
  call io_close(out_file)

  SAFE_DEALLOCATE_A(dielectric)
  SAFE_DEALLOCATE_A(invdielectric)
  SAFE_DEALLOCATE_A(chi)
  SAFE_DEALLOCATE_A(vecpot)
  SAFE_DEALLOCATE_A(vecpot_ref)
  SAFE_DEALLOCATE_A(vecpot0)
  SAFE_DEALLOCATE_A(ftreal)
  SAFE_DEALLOCATE_A(ftimag)
    
  call simul_box_end(sb)
  call geometry_end(geo)
  call space_end(space)
  call io_end()
  call messages_end()
  call global_end()

end program dielectric_function

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
