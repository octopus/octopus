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
use liboct
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

#ifdef HAVE_FFTW
  ! used in method 3
  real(r8), pointer :: ff(:, :, :)   
#endif
end type hartree_type

public :: hartree_type, hartree_init, hartree_solve, hartree_end

contains

subroutine hartree_init(h, m)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(inout) :: m
  
  sub_name = 'hartree_init'; call push_sub()
  
  if(conf%dim==1.or.conf%dim==2) then
    h%solver = -conf%dim ! internal type
    message(1) = 'Info: Using direct integration method to solve poisson equation'
    call write_info(1)
  else
    call oct_parse_int(C_string('PoissonSolver'), 3, h%solver)
    if(h%solver<1 .or. h%solver>3 ) then
      write(message(1), '(a,i2,a)') "Input: '", h%solver, &
           "' is not a valid PoissonSolver"
      message(2) = 'PoissonSolver = 1(cg) | 2(fft) | 3(fft spherical cutoff)'
      call write_fatal(2)
    end if

    if(conf%periodic_dim>0 .and. h%solver.ne.2) then
      message(1) = "Sorry, but for periodic systemps only the fft method"
      message(2) = "is available for solving the poisson equation"
      call write_fatal(2)
    end if

    call hartree3D_init(h, m)
  end if
end subroutine hartree_init

subroutine hartree_end(h)
  type(hartree_type), intent(inout) :: h

  sub_name = 'hartree_end'; call push_sub()

  select case(h%solver)
  case(1)
#ifdef HAVE_FFTW
  case(2,3)
    if(associated(h%ff)) then ! has been allocated => destroy
      deallocate(h%ff); nullify(h%ff)
    end if
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

  sub_name = 'hartree_solve'; call push_sub()

  select case(h%solver)
  case(-1)
    call hartree1d_solve(h, m, pot, dist)
  case(-2)
    call hartree2d_solve(h, m, pot, dist)
  case(1)
    call hartree_cg(h, m, pot, dist)
#ifdef HAVE_FFTW
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
