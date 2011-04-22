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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: dielectric_function.F90 3613 2007-11-29 16:47:41Z xavier $

#include "global.h"

program dielectric_function
  use datasets_m
  use command_line_m
  use geometry_m
  use global_m
  use io_m
  use lalg_adv_m
  use loct_m
  use messages_m
  use parser_m
  use profiling_m
  use space_m
  use spectrum_m
  use simul_box_m
  use unit_m
  use unit_system_m

  implicit none

  integer :: in_file, out_file, ii, jj, kk, idir, jdir, ierr
  integer :: time_steps, energy_steps, istart, iend, ntiter
  FLOAT   :: dt, tt, ww, n0
  FLOAT, allocatable :: vecpot(:, :), dumpa(:), vecpot0(:)
  CMPLX, allocatable :: dielectric(:, :), chi(:, :), invdielectric(:, :), fullmat(:, :)
  type(spec_t)      :: spectrum
  type(block_t)     :: blk
  type(space_t)     :: space
  type(geometry_t)  :: geo
  type(simul_box_t) :: sb

  ! Initialize stuff
  call global_init()
  call getopt_init(ierr)
  if(ierr.eq.0) call getopt_dielectric_function
  call parser_init()
  call messages_init()

  call parse_logical('ExperimentalFeatures', .false., conf%devel_version)

  call datasets_init(1)
  call io_init()

  call unit_system_init()

  call spectrum_init(spectrum)

  call space_init(space)
  call geometry_init(geo, space)
  call simul_box_init(sb, geo, space)
    
  SAFE_ALLOCATE(vecpot0(1:space%dim))

  if(parse_block(datasets_check('GaugeVectorField'), blk) == 0) then
    
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

  in_file = io_open('td.general/gauge_field', action='read', status='old', die=.false.)
  call io_skip_header(in_file)
  call count_time_steps(in_file, time_steps, dt)
  if(in_file < 0) then 
    message(1) = "Cannot open file '"//trim(io_workpath('td.general/gauge_field'))//"'"
    call messages_fatal(1)
  end if

  SAFE_ALLOCATE(vecpot(1:space%dim*3, 1:time_steps))

  call io_skip_header(in_file)
  
  do ii = 1, time_steps
    read(in_file, *) jj, tt, vecpot(1:space%dim*3, ii)
  end do

  call io_close(in_file)

  write(message(1), '(a, i7, a)') "Info: Read ", time_steps, " steps from file '"// &
    trim(io_workpath('td.general/gauge_field'))//"'"
  call messages_info(1)


  ! Find out the iteration numbers corresponding to the time limits.
  call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

  istart = max(1, istart)

  SAFE_ALLOCATE(dumpa(istart:iend))

  do ii = istart, iend
    jj = ii - istart
    select case(spectrum%damp)
    case(SPECTRUM_DAMP_NONE)
      dumpa(ii) = M_ONE
    case(SPECTRUM_DAMP_LORENTZIAN)
      dumpa(ii)= exp(-jj*dt*spectrum%damp_factor)
    case(SPECTRUM_DAMP_POLYNOMIAL)
      dumpa(ii) = M_ONE - M_THREE*(real(jj)/ntiter)**2+ M_TWO*(real(jj)/ntiter)**3
    case(SPECTRUM_DAMP_GAUSSIAN)
      dumpa(ii)= exp(-(jj*dt)**2*spectrum%damp_factor**2)
    end select
  end do

  energy_steps = spectrum%max_energy / spectrum%energy_step
  SAFE_ALLOCATE(invdielectric(1:space%dim, 0:energy_steps))
  SAFE_ALLOCATE(dielectric(1:space%dim, 0:energy_steps))
  SAFE_ALLOCATE(chi(1:space%dim, 0:energy_steps))
  SAFE_ALLOCATE(fullmat(1:space%dim, 1:space%dim))

  n0 = sqrt(sum(vecpot0(1:space%dim))**2)

  do kk = 0, energy_steps
    ww = kk*spectrum%energy_step

    invdielectric(1:space%dim, kk) = M_ZERO

    do ii = istart, iend
      tt = ii*dt
      invdielectric(1:space%dim, kk) = &
        invdielectric(1:space%dim, kk) + vecpot(space%dim + 1:2*space%dim, ii)*exp(M_zI*ww*tt)*dumpa(ii)*dt
    end do
    
    invdielectric(1:space%dim, kk) = (vecpot0(1:space%dim) + invdielectric(1:space%dim, kk))/n0

    ! calculate the full inverse dielectric matrix
    do idir = 1, space%dim
      do jdir = 1, space%dim
        fullmat(idir, mod(jdir + idir - 2, space%dim) + 1) = invdielectric(jdir, kk)
      end do
    end do

    ! and invert it
    call lalg_sym_inverter('u', space%dim, fullmat)

    ! symmetrize the rest dielectric matrix
    do idir = 1, space%dim
      do jdir = 1, jdir - 1
        fullmat(jdir, idir) = fullmat(idir, jdir)
      end do
    end do

    dielectric(1:space%dim, kk) = fullmat(1:space%dim, 1)

    chi(1:space%dim, kk) = (dielectric(1:space%dim, kk) - vecpot0(1:space%dim)/n0)*sb%rcell_volume/(CNST(4.0)*M_PI)
  end do

  SAFE_DEALLOCATE_A(fullmat)

  out_file = io_open('td.general/inverse_dielectric_function', action='write')
  do kk = 0, energy_steps
    ww = kk*spectrum%energy_step
    write(out_file, '(7e15.6)') ww,                                         &
         real(invdielectric(1, kk), REAL_PRECISION), aimag(invdielectric(1, kk)), &
         real(invdielectric(2, kk), REAL_PRECISION), aimag(invdielectric(2, kk)), &
         real(invdielectric(3, kk), REAL_PRECISION), aimag(invdielectric(3, kk))
  end do
  call io_close(out_file)
 
  out_file = io_open('td.general/dielectric_function', action='write')
  do kk = 0, energy_steps
    ww = kk*spectrum%energy_step
    write(out_file, '(7e15.6)') ww,                                         &
         real(dielectric(1, kk), REAL_PRECISION), aimag(dielectric(1, kk)), &
         real(dielectric(2, kk), REAL_PRECISION), aimag(dielectric(2, kk)), &
         real(dielectric(3, kk), REAL_PRECISION), aimag(dielectric(3, kk))
  end do
  call io_close(out_file)

  out_file = io_open('td.general/chi', action='write')
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
  SAFE_DEALLOCATE_A(chi)
  SAFE_DEALLOCATE_A(vecpot)
  SAFE_DEALLOCATE_A(dumpa)

  call simul_box_end(sb)
  call geometry_end(geo)
  call space_end(space)
  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()

end program dielectric_function

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
