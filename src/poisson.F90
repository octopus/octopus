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

#include "global.h"

module poisson
  use global
  use oct_parser
  use mesh
#ifdef HAVE_FFT
  use fft
  use cube_function
#endif
  use functions

  implicit none
  private

  integer :: poisson_solver = -99

  ! used by conjugated gradients (method 1)
  type(mesh_type), pointer :: cg_m_aux

#ifdef HAVE_FFT
  type(dcf) :: fft_cf
  FLOAT, pointer :: fft_coulb_FS(:,:,:)
#endif

public :: poisson_init, poisson_solve, poisson_end

contains

subroutine poisson_init(m)
  type(mesh_type), intent(inout) :: m
  
  if(poisson_solver.ne.-99) return ! already initialized
  
  call push_sub('poisson_init')
  
  if(conf%dim==1.or.conf%dim==2) then
    poisson_solver = -conf%dim ! internal type
    message(1) = 'Info: Using direct integration method to solve poisson equation'
    call write_info(1)
  else
#ifdef HAVE_FFT
    call oct_parse_int('PoissonSolver', 3, poisson_solver)
    if(poisson_solver<1 .or. poisson_solver>3 ) then
      write(message(1), '(a,i2,a)') "Input: '", poisson_solver, &
           "' is not a valid PoissonSolver"
      message(2) = 'PoissonSolver = 1(cg) | 2(fft) | 3(fft spherical cutoff)'
      call write_fatal(2)
    end if

    if(conf%periodic_dim>0 .and. poisson_solver.ne.2) then
      message(1) = "For periodic systems we cannot use"
      message(2) = "a spherical cutoff for the Poisson solver."
      message(3) = "PoissonSolver = 2 will be used."
      call write_warning(3)
      poisson_solver = 2
    end if
#else
    poisson_solver = 1
#endif

    call poisson3D_init(m)
  end if
end subroutine poisson_init

subroutine poisson_end()
  integer :: j
  call push_sub('poisson_end')
  select case(poisson_solver)
  case(1)
    call derivatives_end(cg_m_aux%laplacian)
    do j = 1, conf%dim
       call derivatives_end(cg_m_aux%grad(j))
    enddo
    deallocate(cg_m_aux%grad); nullify(cg_m_aux%grad)

    call mesh_end(cg_m_aux)
    deallocate(cg_m_aux); nullify(cg_m_aux)
#ifdef HAVE_FFT
  case(2,3)
    call dcf_free(fft_cf)
    deallocate(fft_coulb_FS); nullify(fft_coulb_FS)
#endif
  end select

  call pop_sub()
  return
end subroutine poisson_end

subroutine poisson_solve(m, pot, dist)
  type(mesh_type), intent(IN) :: m
  FLOAT, intent(inout) :: pot(m%np)
  FLOAT, intent(IN)    :: dist(m%np)

  call push_sub('poisson_solve')

  ASSERT(poisson_solver.ne.-99)

  select case(poisson_solver)
  case(-1)
    call poisson1d_solve(m, pot, dist)
  case(-2)
    call poisson2d_solve(m, pot, dist)
  case(1)
    call poisson_cg(m, pot, dist)
#ifdef HAVE_FFT
  case(2,3)
    call poisson_fft(m, pot, dist)
#endif
  end select

  call pop_sub()
end subroutine poisson_solve

#include "poisson1D.F90"
#include "poisson2D.F90"
#include "poisson3D.F90"

end module poisson
