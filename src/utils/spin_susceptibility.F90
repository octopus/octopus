!! Copyright (C) 2018 N. Tancogne-Dejean
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
  use global_oct_m
  use io_oct_m
  use ions_oct_m
  use kick_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use spectrum_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  integer :: in_file, out_file, ref_file, ii, jj, kk, ierr, ib, num_col, num_col_cart
  integer :: time_steps, time_steps_ref, energy_steps, istart, iend, ntiter, iq
  FLOAT   :: dt, tt, ww, dt_ref
  FLOAT, allocatable :: ftreal(:,:), ftimag(:,:)
  FLOAT, allocatable :: m_cart(:,:), m_cart_ref(:,:), magnetization(:,:,:)
  type(spectrum_t)  :: spectrum
  type(batch_t)     :: vecpotb, ftrealb, ftimagb
  type(kick_t)      :: kick
  character(len=256) :: header
  character(len=MAX_PATH_LEN) :: fname, ref_filename
  CMPLX, allocatable :: chi(:,:)

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

  in_file = io_open('td.general/total_magnetization', global_namespace, action='read', status='old', die=.false.)
  if(in_file < 0) then
    message(1) = "Cannot open file '"//trim(io_workpath('td.general/total_magnetization', global_namespace))//"'"
    call messages_fatal(1)
  end if
  rewind(in_file)
  read(in_file,*)
  read(in_file,*)
  call kick_read(kick, in_file, global_namespace)
  call io_skip_header(in_file)
  call spectrum_count_time_steps(global_namespace, in_file, time_steps, dt)

  time_steps = time_steps + 1

  num_col_cart = 12*kick%nqvec
  SAFE_ALLOCATE(m_cart(1:time_steps, 1:num_col_cart))

  call io_skip_header(in_file)

  do ii = 1, time_steps
    read(in_file, *) jj, tt, (m_cart(ii,ib),ib = 1, num_col_cart)
  end do

  call io_close(in_file)

  if(parse_is_defined(global_namespace, 'TransientMagnetizationReference')) then
    !%Variable TransientMagnetizationReference
    !%Type string
    !%Default "."
    !%Section Utilities::oct-spin_susceptibility
    !%Description
    !% In case of delayed kick, the calculation of the transient spin susceptibility requires 
    !% to substract a reference calculation, containing dynamics of the magnetization without the kick
    !% This reference must be computed having
    !% TDOutput = total_magnetization.
    !% This variables defined the directory in which the reference total_magnetization file is,
    !% relative to the current folder
    !%End

    call parse_variable(global_namespace, 'TransientMagnetizationReference', '.', ref_filename)
    ref_file = io_open(trim(ref_filename)//'/total_magnetization', global_namespace, action='read', status='old', die=.false.)
    if(ref_file < 0) then
      message(1) = "Cannot open reference file '"// &
        trim(io_workpath(trim(ref_filename)//'/total_magnetization', global_namespace))//"'"
      call messages_fatal(1)
    end if
    call io_skip_header(ref_file)
    call spectrum_count_time_steps(global_namespace, ref_file, time_steps_ref, dt_ref)
    time_steps_ref = time_steps_ref + 1

    if(time_steps_ref < time_steps) then
      message(1) = "The reference calculation does not contain enought time steps"
      call messages_fatal(1)
    end if

    if(dt_ref /= dt) then
      message(1) = "The time step of the reference calculation is different from the current calculation"
      call messages_fatal(1)
    end if

    !We remove the reference
    SAFE_ALLOCATE(m_cart_ref(1:time_steps_ref, num_col_cart))
    call io_skip_header(ref_file)
    do ii = 1, time_steps_ref
      read(ref_file, *) jj, tt, (m_cart_ref(ii, kk), kk = 1, num_col_cart)
    end do
    call io_close(ref_file)
    do ii = 1, time_steps
      do kk = 1, num_col_cart
        m_cart(ii, kk) = m_cart(ii, kk) - m_cart_ref(ii, kk)
      end do
    end do
    SAFE_DEALLOCATE_A(m_cart_ref)
  end if

  !We now perform the change of basis to the rotating basis
  !In this basis we have only m_+(q), m_-(q), and m_z(+/-q)
  !where z means here along the easy axis
  num_col = 8
  SAFE_ALLOCATE(magnetization(1:time_steps, 1:num_col, 1:kick%nqvec))

  do iq = 1, kick%nqvec
    !Real part of m_x
    magnetization(:,1,iq) = m_cart(:,(iq-1)*12+1)*kick%trans_vec(1,1)  &
                          + m_cart(:,(iq-1)*12+3)*kick%trans_vec(2,1) &
                          + m_cart(:,(iq-1)*12+5)*kick%trans_vec(3,1)
    !We add -Im(m_y)
    magnetization(:,1,iq) = magnetization(:,1,iq) -(m_cart(:,(iq-1)*12+2)*kick%trans_vec(1,2) &
                           + m_cart(:,(iq-1)*12+4)*kick%trans_vec(2,2) &
                           + m_cart(:,(iq-1)*12+6)*kick%trans_vec(3,2))

    !Im part of m_x
    magnetization(:,2,iq) = m_cart(:,(iq-1)*12+2)*kick%trans_vec(1,1) &
                          + m_cart(:,(iq-1)*12+4)*kick%trans_vec(2,1) &
                          + m_cart(:,(iq-1)*12+6)*kick%trans_vec(3,1)
    !We add +Re(m_y)
    magnetization(:,2,iq) = magnetization(:,2,iq) +(m_cart(:,(iq-1)*12+1)*kick%trans_vec(1,2) &
                         + m_cart(:,(iq-1)*12+3)*kick%trans_vec(2,2) &
                         + m_cart(:,(iq-1)*12+5)*kick%trans_vec(3,2))

    !Real part of m_x
    magnetization(:,3,iq) = m_cart(:,(iq-1)*12+7)*kick%trans_vec(1,1)  &
                           + m_cart(:,(iq-1)*12+9)*kick%trans_vec(2,1) &
                           + m_cart(:,(iq-1)*12+11)*kick%trans_vec(3,1)
    !We add +Im(m_y)
    magnetization(:,3,iq) = magnetization(:,3,iq) +(m_cart(:,(iq-1)*12+8)*kick%trans_vec(1,2) &
                           + m_cart(:,(iq-1)*12+10)*kick%trans_vec(2,2) &
                           + m_cart(:,(iq-1)*12+12)*kick%trans_vec(3,2))

    !Im part of m_x
    magnetization(:,4,iq) = m_cart(:,(iq-1)*12+8)*kick%trans_vec(1,1)  &
                          + m_cart(:,(iq-1)*12+10)*kick%trans_vec(2,1) &
                          + m_cart(:,(iq-1)*12+12)*kick%trans_vec(3,1)
    !We add -Re(m_y)
    magnetization(:,4,iq) = magnetization(:,4,iq) -(m_cart(:,(iq-1)*12+7)*kick%trans_vec(1,2) &
                         + m_cart(:,(iq-1)*12+9)*kick%trans_vec(2,2) &
                         + m_cart(:,(iq-1)*12+11)*kick%trans_vec(3,2))

    !Real and Im part of m_z
    magnetization(:,5,iq) = m_cart(:,(iq-1)*12+1)*kick%easy_axis(1)  &
                          + m_cart(:,(iq-1)*12+3)*kick%easy_axis(2)  &
                          + m_cart(:,(iq-1)*12+5)*kick%easy_axis(3)
    magnetization(:,6,iq) = m_cart(:,(iq-1)*12+2)*kick%easy_axis(1)  &
                          + m_cart(:,(iq-1)*12+4)*kick%easy_axis(2)  &
                          + m_cart(:,(iq-1)*12+6)*kick%easy_axis(3)
    magnetization(:,7,iq) = m_cart(:,(iq-1)*12+7)*kick%easy_axis(1)  &
                          + m_cart(:,(iq-1)*12+9)*kick%easy_axis(2)  &
                          + m_cart(:,(iq-1)*12+11)*kick%easy_axis(3)
    magnetization(:,8,iq) = m_cart(:,(iq-1)*12+8)*kick%easy_axis(1)  &
                          + m_cart(:,(iq-1)*12+10)*kick%easy_axis(2) &
                          + m_cart(:,(iq-1)*12+12)*kick%easy_axis(3)
  end do

  SAFE_DEALLOCATE_A(m_cart)

  write(message(1), '(a, i7, a)') "Info: Read ", time_steps, " steps from file '"// &
    trim(io_workpath('td.general/total_magnetization', global_namespace))//"'"
  call messages_info(1)

  write(header, '(9a15)') '#  time', 'Re[m_+(q,t)]', 'Im[m_+(q,t)]', &
            'Re[m_-(-q, t)]', 'Im[m_-(-q,t)]', &
            'Re[m_z(q,t)]', 'Im[m_z(q,t)]', 'Re[m_z(-q,t)]', 'Im[m_z(-q,t)]'

  do iq = 1, kick%nqvec

    write(fname, '(a,i3.3)') 'td.general/transverse_magnetization_q', iq
    out_file = io_open(trim(fname), global_namespace, action='write')
    write(out_file,'(a)') trim(header)
    do kk = 1, time_steps
      write(out_file, '(9e15.6)') (kk - 1)*dt, (magnetization(kk,ii,iq), ii = 1,8)
    end do
    call io_close(out_file)


    ! Find out the iteration numbers corresponding to the time limits.
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    istart = max(1, istart)

    energy_steps = spectrum_nenergy_steps(spectrum)

    SAFE_ALLOCATE(ftreal(1:energy_steps, num_col))
    SAFE_ALLOCATE(ftimag(1:energy_steps, num_col))

    call batch_init(vecpotb, 1, 1, num_col, magnetization(:, :, iq))
    call batch_init(ftrealb, 1, 1, num_col, ftreal)
    call batch_init(ftimagb, 1, 1, num_col, ftimag)


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

    ASSERT(abs(anint(num_col*M_HALF) - num_col*M_HALF) <= M_EPSILON)
    SAFE_ALLOCATE(chi(1:energy_steps, 1:num_col/2))
    do ii = 1, num_col/2
      do kk = 1, energy_steps
        chi(kk,ii) = (ftreal(kk,(ii-1)*2+1) + M_zI*ftimag(kk, (ii-1)*2+1)&
                    -ftimag(kk, (ii-1)*2+2) + M_zI*ftreal(kk, (ii-1)*2+2))/kick%delta_strength
      end do
    end do

    SAFE_DEALLOCATE_A(ftreal)
    SAFE_DEALLOCATE_A(ftimag)


    write(header, '(9a15)') '#        energy', 'Re[\chi_{+-}(q)]', 'Im[\chi_{+-}(q)]', &
              'Re[\chi_{-+}(-q)]', 'Im[\chi_{-+}(-q)]', &
              'Re[\chi_{zz}(q)]', 'Im[\chi_{zz}(q)]', 'Re[\chi_{zz}(-q)]', 'Im[\chi_{zz}(-q)]'

    write(fname, '(a,i3.3)') 'td.general/spin_susceptibility_q', iq
    out_file = io_open(trim(fname), global_namespace, action='write')
    write(out_file,'(a)') trim(header)
    do kk = 1, energy_steps
      ww = (kk-1)*spectrum%energy_step + spectrum%min_energy
      write(out_file, '(13e15.6)') ww,                                   &
               (TOFLOAT(chi(kk,ii)), aimag(chi(kk,ii)), ii = 1, num_col/2)
    end do
    call io_close(out_file)

    SAFE_DEALLOCATE_A(chi)
  end do !iq

  SAFE_DEALLOCATE_A(magnetization)

  call io_end()
  call messages_end()
  call parser_end()
  call global_end()

end program spin_susceptibility

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
