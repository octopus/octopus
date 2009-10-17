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

program broad
  use datasets_m
  use global_m
  use io_m
  use parser_m
  use messages_m
  use profiling_m
  use units_m

  implicit none

  type broad_t
    FLOAT :: b, energy_step, min_energy, max_energy
  end type broad_t

  type(broad_t) :: b

  ! Initialize stuff
  call global_init()
  call parser_init()
  call datasets_init(1)
  call io_init()
  if(in_debug_mode) then
     call io_mkdir('debug')
  end if
  call units_init()

  ! broadening to use
  call parse_float(datasets_check('LinBroadening'),  units_from_atomic(units_inp%energy, CNST(0.005)), b%b)
  call parse_float(datasets_check('LinEnergyStep'),  units_from_atomic(units_inp%energy, CNST(0.001)), b%energy_step)
  call parse_float(datasets_check('LinMinEnergy'),   M_ZERO, b%min_energy)
  call parse_float(datasets_check('LinMaxEnergy'),   units_from_atomic(units_inp%energy, M_ONE), b%max_energy)
  b%b = units_to_atomic(units_inp%energy, b%b)
  b%energy_step = units_to_atomic(units_inp%energy, b%energy_step)
  b%min_energy  = units_to_atomic(units_inp%energy, b%min_energy)
  b%max_energy  = units_to_atomic(units_inp%energy, b%max_energy)

  call calc_broad(b, CASIDA_DIR, 'eps-diff', .true.)
  call calc_broad(b, CASIDA_DIR, 'petersilka', .true.)
  call calc_broad(b, CASIDA_DIR, 'casida', .false.)

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()

contains

  ! ---------------------------------------------------------
  subroutine calc_broad(b, dir, fname, extracols)
    type(broad_t), intent(in) :: b
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: fname
    logical, intent(in) :: extracols

    FLOAT, allocatable :: s(:,:)
    FLOAT :: w, e, f(4)
    integer :: n, iunit, j1, j2

    n = (b%max_energy - b%min_energy) / b%energy_step
    SAFE_ALLOCATE(s(1:4, 1:n))
    s = M_ZERO

    iunit = io_open(trim(dir)//"/"// fname, action='read', status='old', die = .false.)
    if(iunit < 0) return

    read(iunit, *) ! skip header
    do
      if(extracols) then
        read(iunit, *, end=100) j1, j2, e, f
      else
        read(iunit, *, end=100) e, f
      end if

      e = units_to_atomic(units_inp%energy, e)

      do j1 = 1, n
        w = b%min_energy + real(j1-1, REAL_PRECISION)*b%energy_step
        s(:, j1) = s(:, j1) + f(:)*b%b/((w-e)**2 + b%b**2)/M_PI ! lorentzian
      end do
    end do
100   continue
    call io_close(iunit)

    ! print spectra
    iunit = io_open(trim(dir)//"/spectrum."//fname, action='write')
    do j1 = 1, n
      w = units_from_atomic(units_inp%energy, b%min_energy + real(j1 - 1, REAL_PRECISION)*b%energy_step)
      s(:, j1) = units_from_atomic(units_inp%energy, s(:, j1))
      write(iunit, '(5es14.6)') w, s(:, j1)
    end do

    ! IMPORTANT: s has changed units here

    call io_close(iunit)

    SAFE_DEALLOCATE_A(s)
  end subroutine calc_broad

end program broad

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
