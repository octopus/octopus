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
  use namespace_oct_m
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

  call parser_init()

  call messages_init()

  call io_init()

  call unit_system_init(global_namespace)

  call spectrum_init(spectrum, global_namespace)

  call space_init(space, global_namespace)
  call geometry_init(geo, global_namespace, space)
  call simul_box_init(sb, global_namespace, geo, space)
    
  SAFE_ALLOCATE(vecpot0(1:space%dim))

  if(parse_block(global_namespace, 'GaugeVectorField', blk) == 0) then
    
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
  call parse_variable(global_namespace, 'GaugeFieldDelay', start_time, spectrum%start_time )

  in_file = io_open('td.general/gauge_field', global_namespace, action='read', status='old', die=.false.)
  if(in_file < 0) then 
    message(1) = "Cannot open file '"//trim(io_workpath('td.general/gauge_field', global_namespace))//"'"
    call messages_fatal(1)
  end if
  call io_skip_header(in_file)
  call spectrum_count_time_steps(global_namespace, in_file, time_steps, dt)

  if(parse_is_defined(global_namespace, 'TransientAbsorptionReference')) then
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

    call parse_variable(global_namespace, 'TransientAbsorptionReference', '.', ref_filename)
    ref_file = io_open(trim(ref_filename)//'/gauge_field', global_namespace, action='read', status='old', die=.false.)
    if(ref_file < 0) then
      message(1) = "Cannot open reference file '"// &
        trim(io_workpath(trim(ref_filename)//'/gauge_field', global_namespace))//"'"
      call messages_fatal(1)
    end if
    call io_skip_header(ref_file)
    call spectrum_count_time_steps(global_namespace, ref_file, time_steps_ref, dt_ref)
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
  if(parse_is_defined(global_namespace, 'TransientAbsorptionReference')) then
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
    trim(io_workpath('td.general/gauge_field', global_namespace))//"'"
  call messages_info(1)


  ! Find out the iteration numbers corresponding to the time limits.
  call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

  istart = max(1, istart)

  energy_steps = spectrum_nenergy_steps(spectrum)

  n0 = sqrt(sum(vecpot0(1:space%dim))**2)

  SAFE_ALLOCATE(ftreal(1:energy_steps, 1:space%dim))
  SAFE_ALLOCATE(ftimag(1:energy_steps, 1:space%dim))

  call batch_init(vecpotb, 1, 1, space%dim, vecpot(:, space%dim+1:space%dim*2))
  call batch_init(ftrealb, 1, 1, space%dim, ftreal)
  call batch_init(ftimagb, 1, 1, space%dim, ftimag)

  call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart, iend, spectrum%start_time, dt, vecpotb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, vecpotb, spectrum%min_energy, &
    spectrum%max_energy, spectrum%energy_step, ftrealb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, vecpotb, spectrum%min_energy, &
    spectrum%max_energy, spectrum%energy_step, ftimagb)


  call vecpotb%end()
  call ftrealb%end()
  call ftimagb%end()

  SAFE_ALLOCATE(invdielectric(1:space%dim, 1:energy_steps))
  SAFE_ALLOCATE(dielectric(1:space%dim, 1:energy_steps))
  SAFE_ALLOCATE(chi(1:space%dim, 1:energy_steps))
  SAFE_ALLOCATE(fullmat(1:space%dim, 1:space%dim))

  do kk = 1, energy_steps
    ww = (kk-1)*spectrum%energy_step + spectrum%min_energy

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

  out_file = io_open('td.general/inverse_dielectric_function', global_namespace, action='write')
  select case(space%dim)
  case(1)
    write(header, '(7a15)') '#        energy', 'Re x', 'Im x'
  case(2)
    write(header, '(7a15)') '#        energy', 'Re x', 'Im x', 'Re y', 'Im y'
  case(3)
    write(header, '(7a15)') '#        energy', 'Re x', 'Im x', 'Re y', 'Im y', 'Re z', 'Im z'
  end select


  write(out_file,'(a)') trim(header)
  do kk = 1, energy_steps
    ww = (kk-1)*spectrum%energy_step + spectrum%min_energy
    write(out_file, '(e15.6)', advance='no') ww
    do idir = 1, space%dim
      write(out_file, '(2e15.6)', advance='no') TOFLOAT(invdielectric(idir, kk)), aimag(invdielectric(idir, kk))
    end do
    write(out_file, '()')
  end do
  call io_close(out_file)
 
  out_file = io_open('td.general/dielectric_function', global_namespace, action='write')
  write(out_file,'(a)') trim(header)
  do kk = 1, energy_steps
    ww = (kk-1)*spectrum%energy_step + spectrum%min_energy
    write(out_file, '(e15.6)', advance='no') ww
    do idir = 1, space%dim
      write(out_file, '(2e15.6)', advance='no') TOFLOAT(dielectric(idir, kk)), aimag(dielectric(idir, kk))
    end do
    write(out_file, '()')
  end do
  call io_close(out_file)

  out_file = io_open('td.general/chi', global_namespace, action='write')
  write(out_file,'(a)') trim(header)
  do kk = 1, energy_steps
    dielectric(1:3, kk) = (dielectric(1:3, kk) - M_ONE)/(CNST(4.0)*M_PI)
    ww = (kk-1)*spectrum%energy_step + spectrum%min_energy
    write(out_file, '(e15.6)', advance='no') ww
    do idir = 1, space%dim
      write(out_file, '(2e15.6)', advance='no') TOFLOAT(chi(idir, kk)), aimag(chi(idir, kk))
    end do
    write(out_file, '()')
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

  call parser_end()
  call global_end()

end program dielectric_function

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
