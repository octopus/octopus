!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "config_F90.h"

module hartree
use global
use oct_parser
use mesh
use math
#ifdef HAVE_FFTW
use fft
#endif
implicit none

private

type hartree_type
  integer :: solver

  ! used by conjugated gradients (method 1)
  integer :: ncgiter
  type(mesh_type), pointer :: m_aux

#if defined(HAVE_FFTW) || defined(HAVE_FFTW3)
  ! used by the fft methods
  type(fft_type) :: fft
  real(r8), pointer :: ff(:, :, :)   
#endif
end type hartree_type

public :: hartree_type, hartree_init, hartree_solve, hartree_end

contains

subroutine hartree_init(h, m)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(inout) :: m
  
  call push_sub('hartree_init')
  
  if(conf%dim==1.or.conf%dim==2) then
    h%solver = -conf%dim ! internal type
    message(1) = 'Info: Using direct integration method to solve poisson equation'
    call write_info(1)
  else
    call oct_parse_int('PoissonSolver', 3, h%solver)
    if(h%solver<1 .or. h%solver>3 ) then
      write(message(1), '(a,i2,a)') "Input: '", h%solver, &
           "' is not a valid PoissonSolver"
      message(2) = 'PoissonSolver = 1(cg) | 2(fft) | 3(fft spherical cutoff)'
      call write_fatal(2)
    end if

    if(conf%periodic_dim>0 .and. h%solver.ne.2) then
      message(1) = "For periodic systems we cannot use"
      message(2) = "a spherical cutoff for the Poisson solver."
      message(3) = "PoissonSolver = 2 will be used."
      call write_warning(3)
      h%solver = 2
    end if

    call hartree3D_init(h, m)
  end if
end subroutine hartree_init

subroutine hartree_end(h)
  type(hartree_type), intent(inout) :: h

  call push_sub('hartree_end')

  select case(h%solver)
  case(1)
    nullify(h%m_aux%d) ! this is a copy from sys%m, so it should not be dealloated here
    call mesh_end(h%m_aux)
    deallocate(h%m_aux); nullify(h%m_aux)
#if defined(HAVE_FFTW) || defined(HAVE_FFTW3)
  case(2,3)
    call fft_end(h%fft)
    deallocate(h%ff); nullify(h%ff)
#endif
  end select

  call pop_sub()
  return
end subroutine hartree_end

subroutine hartree_solve(h, m, pot, dist)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(inout) :: pot
  real(r8), dimension(:, :), intent(IN) :: dist

  call push_sub('hartree_solve')

  select case(h%solver)
  case(-1)
    call hartree1d_solve(h, m, pot, dist)
  case(-2)
    call hartree2d_solve(h, m, pot, dist)
  case(1)
    call hartree_cg(h, m, pot, dist)
#if defined(HAVE_FFTW) || defined(HAVE_FFTW3)
  case(2,3)
    call hartree_fft(h, m, pot, dist)
#endif

  case default
    message(1) = "Hartree structure not initialized"
    call write_fatal(1)
  end select

  call pop_sub()
  return
end subroutine hartree_solve

#include "hartree1D.F90"
#include "hartree2D.F90"
#include "hartree3D.F90"

end module hartree
