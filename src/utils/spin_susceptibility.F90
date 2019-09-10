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
  use kick_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
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

  integer :: in_file, out_file, ref_file, ii, jj, kk, idir, jdir, ierr, ib, num_col, num_col_cart
  integer :: time_steps, time_steps_ref, energy_steps, istart, iend, ntiter, n_rows
  FLOAT   :: dt, dt_ref, tt, ww, norm, dot
  FLOAT, allocatable :: ftreal(:,:), ftimag(:,:)
  FLOAT, allocatable :: m_cart(:,:), magnetization(:,:)
  type(spectrum_t)  :: spectrum
  type(block_t)     :: blk
  type(batch_t)     :: vecpotb, ftrealb, ftimagb
  type(namespace_t) :: namespace
  type(kick_t)      :: kick
  character(len=256) :: header
  FLOAT :: delta_strength
  character(len=MAX_PATH_LEN) :: ref_filename
  CMPLX, allocatable :: chi(:,:)

  ! Initialize stuff
  call global_init(is_serial = .true.)

  call getopt_init(ierr)
  if(ierr == 0) call getopt_dielectric_function()
  call getopt_end()

  call parser_init()
  namespace = namespace_t("")


  call messages_init(namespace)

  call io_init(namespace)

  call unit_system_init(namespace)

  call spectrum_init(spectrum, namespace)

  call parse_variable(namespace, 'TDDeltaStrength', M_ZERO, delta_strength )

  if(parse_block(namespace, 'TDEasyAxis', blk)==0) then
    n_rows = parse_block_n(blk)

    do idir = 1, 3
      call parse_block_float(blk, 0, idir - 1, kick%easy_axis(idir))
    end do
    norm = sqrt(sum(kick%easy_axis(1:3)**2))
    if(norm < CNST(1e-9)) then
      message(1) = "TDEasyAxis norm is too small."
      call messages_fatal(1)
    end if
    kick%easy_axis(1:3) = kick%easy_axis(1:3)/norm
    call parse_block_end(blk)
  else
    message(1) = "For magnons, the variable TDEasyAxis must be defined."
    call messages_fatal(1)
  end if

  !We first two vectors defining a basis in the transverse plane
  !For this we take two vectors not align with the first one
  !and we perform a Gram-Schmidt orthogonalization
  kick%trans_vec(1,1) = -kick%easy_axis(2)
  kick%trans_vec(2,1) = M_TWO*kick%easy_axis(3)
  kick%trans_vec(3,1) = M_THREE*kick%easy_axis(1)

  dot = sum(kick%easy_axis(1:3)*kick%trans_vec(1:3,1))
  kick%trans_vec(1:3,1) = kick%trans_vec(1:3,1) - dot*kick%easy_axis(1:3)
  norm = sum(kick%trans_vec(1:3,1)**2)
  kick%trans_vec(1:3,1) = kick%trans_vec(1:3,1)/sqrt(norm)

  kick%trans_vec(1,2) = kick%easy_axis(2) * kick%trans_vec(3,1) - kick%easy_axis(3) * kick%trans_vec(2,1)
  kick%trans_vec(2,2) = kick%easy_axis(3) * kick%trans_vec(1,1) - kick%easy_axis(1) * kick%trans_vec(3,1)
  kick%trans_vec(3,2) = kick%easy_axis(1) * kick%trans_vec(2,1) - kick%easy_axis(2) * kick%trans_vec(1,1)


  in_file = io_open('td.general/total_magnetization', namespace, action='read', status='old', die=.false.)
  if(in_file < 0) then 
    message(1) = "Cannot open file '"//trim(io_workpath('td.general/total_magnetization', namespace))//"'"
    call messages_fatal(1)
  end if
  call io_skip_header(in_file)
  call spectrum_count_time_steps(in_file, time_steps, dt)

  time_steps = time_steps + 1

  num_col_cart = 12
  SAFE_ALLOCATE(m_cart(1:time_steps, 1:num_col_cart))
  
  call io_skip_header(in_file)
  
  do ii = 1, time_steps
    read(in_file, *) jj, tt, (m_cart(ii,ib),ib = 1, num_col_cart)
  end do

  call io_close(in_file)

  !We now perform the change of basis to the rotating basis
  !In this basis we have only m_+(q), m_-(q), and m_z(+/-q)
  !where z means here along the easy axis
  num_col = 8
  SAFE_ALLOCATE(magnetization(1:time_steps, 1:num_col))
  !Real part of m_x
  magnetization(:,1) = m_cart(:,1)*kick%trans_vec(1,1) + m_cart(:,3)*kick%trans_vec(2,1) &
                         + m_cart(:,5)*kick%trans_vec(3,1)
  !We add -Im(m_y)
  magnetization(:,1) = magnetization(:,1) -(m_cart(:,2)*kick%trans_vec(1,2) &
                         + m_cart(:,4)*kick%trans_vec(2,2) + m_cart(:,6)*kick%trans_vec(3,2))  
  
  !Im part of m_x
  magnetization(:,2) = m_cart(:,2)*kick%trans_vec(1,1) + m_cart(:,4)*kick%trans_vec(2,1) &
                         + m_cart(:,6)*kick%trans_vec(3,1)
  !We add +Re(m_y)
  magnetization(:,2) = magnetization(:,2) +(m_cart(:,1)*kick%trans_vec(1,2) &
                         + m_cart(:,3)*kick%trans_vec(2,2) + m_cart(:,5)*kick%trans_vec(3,2))

  !Real part of m_x
  magnetization(:,3) = m_cart(:,7)*kick%trans_vec(1,1) + m_cart(:,9)*kick%trans_vec(2,1) &
                         + m_cart(:,11)*kick%trans_vec(3,1)
  !We add +Im(m_y)
  magnetization(:,3) = magnetization(:,3) +(m_cart(:,8)*kick%trans_vec(1,2) &
                         + m_cart(:,10)*kick%trans_vec(2,2) + m_cart(:,12)*kick%trans_vec(3,2))

  !Im part of m_x
  magnetization(:,4) = m_cart(:,8)*kick%trans_vec(1,1) + m_cart(:,10)*kick%trans_vec(2,1) &
                         + m_cart(:,12)*kick%trans_vec(3,1)
  !We add -Re(m_y)
  magnetization(:,4) = magnetization(:,4) -(m_cart(:,7)*kick%trans_vec(1,2) &
                         + m_cart(:,9)*kick%trans_vec(2,2) + m_cart(:,11)*kick%trans_vec(3,2))
  
  !Real and Im part of m_z
  magnetization(:,5) = m_cart(:,1)*kick%easy_axis(1) + m_cart(:,3)*kick%easy_axis(2) &
                              + m_cart(:,5)*kick%easy_axis(3)
  magnetization(:,6) = m_cart(:,2)*kick%easy_axis(1) + m_cart(:,4)*kick%easy_axis(2) &
                              + m_cart(:,6)*kick%easy_axis(3)
  magnetization(:,7) = m_cart(:,7)*kick%easy_axis(1) + m_cart(:,9)*kick%easy_axis(2) &
                              + m_cart(:,11)*kick%easy_axis(3)
  magnetization(:,8) = m_cart(:,8)*kick%easy_axis(1) + m_cart(:,10)*kick%easy_axis(2) &
                              + m_cart(:,12)*kick%easy_axis(3)

  SAFE_DEALLOCATE_A(m_cart)

  write(message(1), '(a, i7, a)') "Info: Read ", time_steps, " steps from file '"// &
    trim(io_workpath('td.general/total_magnetization', namespace))//"'"
  call messages_info(1)

  write(header, '(9a15)') '#  time', 'Re[m_+(q,t)]', 'Im[m_+(q,t)]', &
            'Re[m_-(-q, t)]', 'Im[m_-(-q,t)]', &
            'Re[m_z(q,t)]', 'Im[m_z(q,t)]', 'Re[m_z(-q,t)]', 'Im[m_z(-q,t)]'

  out_file = io_open('td.general/tranverse_magnetization', namespace, action='write')
  write(out_file,'(a)') trim(header)
  do kk = 1, time_steps
    write(out_file, '(9e15.6)') (kk - 1)*dt, (magnetization(kk,ii), ii = 1, 8)
  end do
  call io_close(out_file)


  ! Find out the iteration numbers corresponding to the time limits.
  call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

  istart = max(1, istart)

  energy_steps = spectrum_nenergy_steps(spectrum)

  SAFE_ALLOCATE(ftreal(1:energy_steps, num_col))
  SAFE_ALLOCATE(ftimag(1:energy_steps, num_col))

  call batch_init(vecpotb, num_col)
  call batch_init(ftrealb, num_col)
  call batch_init(ftimagb, num_col)

  do ib = 1, num_col
    call batch_add_state(vecpotb, magnetization(:, ib))
    call batch_add_state(ftrealb, ftreal(1:energy_steps,ib))
    call batch_add_state(ftimagb, ftimag(1:energy_steps,ib))
  end do


  call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart, iend, spectrum%start_time, dt, vecpotb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, vecpotb, spectrum%min_energy, &
    spectrum%max_energy, spectrum%energy_step, ftrealb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, vecpotb, spectrum%min_energy, &
    spectrum%max_energy, spectrum%energy_step, ftimagb)


  call batch_end(vecpotb)
  call batch_end(ftrealb)
  call batch_end(ftimagb)
  SAFE_DEALLOCATE_A(magnetization)

  ASSERT(abs(anint(num_col*M_HALF) - num_col*M_HALF) <= M_EPSILON)
  SAFE_ALLOCATE(chi(1:energy_steps, 1:num_col/2))
  do ii = 1, num_col/2
    do kk = 1, energy_steps
      chi(kk,ii) = (ftreal(kk,(ii-1)*2+1) + M_zI*ftimag(kk, (ii-1)*2+1)&
                  -ftimag(kk, (ii-1)*2+2) + M_zI*ftreal(kk, (ii-1)*2+2))/delta_strength
    end do
  end do

  SAFE_DEALLOCATE_A(ftreal)
  SAFE_DEALLOCATE_A(ftimag)
  

  write(header, '(9a15)') '#        energy', 'Re[\chi_{+-}(q)]', 'Im[\chi_{+-}(q)]', &
            'Re[\chi_{-+}(-q)]', 'Im[\chi_{-+}(-q)]', &
            'Re[\chi_{zz}(q)]', 'Im[\chi_{zz}(q)]', 'Re[\chi_{zz}(-q)]', 'Im[\chi_{zz}(-q)]'

  out_file = io_open('td.general/spin_susceptibility', namespace, action='write')
  write(out_file,'(a)') trim(header)
  do kk = 1, energy_steps
    ww = (kk-1)*spectrum%energy_step + spectrum%min_energy
    write(out_file, '(13e15.6)') ww,                                   &
             (real(chi(kk,ii), REAL_PRECISION), aimag(chi(kk,ii)), ii = 1, num_col/2)
  end do
  call io_close(out_file)
 
  SAFE_DEALLOCATE_A(chi)

  call io_end()
  call messages_end()
  call parser_end()
  call global_end()

end program spin_susceptibility

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
