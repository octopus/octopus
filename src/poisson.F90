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
  use lib_oct_parser
  use geometry
  use mesh
#ifdef HAVE_FFT
  use fft
  use cube_function
#endif
  use functions
  use math

  implicit none
  private

  integer :: poisson_solver = -99

#ifdef HAVE_FFT
  type(dcf) :: fft_cf
  FLOAT, pointer :: fft_coulb_FS(:,:,:)

  integer, parameter :: FFT_SPH       = 0, &
                        FFT_CYL       = 1, &
                        FFT_PLA       = 2, &
                        FFT_NOCUT     = 3
#endif

  integer, parameter :: CG            = 4

public :: poisson_init, poisson_solve, poisson_end

contains

subroutine poisson_init(m)
  type(mesh_type),     intent(inout) :: m
  
  if(poisson_solver.ne.-99) return ! already initialized
  
  call push_sub('poisson_init')
  
  if(conf%dim == 1 .or. conf%dim == 2) then
    poisson_solver = -conf%dim ! internal type
    message(1) = 'Info: Using direct integration method to solve poisson equation'
    call write_info(1)
  else
#ifdef HAVE_FFT
    call loct_parse_int('PoissonSolver', conf%periodic_dim, poisson_solver)
    if(poisson_solver < 0 .or. poisson_solver > 4 ) then
      write(message(1), '(a,i2,a)') "Input: '", poisson_solver, &
           "' is not a valid PoissonSolver"
      message(2) = 'PoissonSolver = 0 (fft spherical cutoff)   | '
      message(3) = '                1 (fft cylindrical cutoff) | '
      message(4) = '                2 (fft planar cutoff)      | '
      message(5) = '                3 (fft no cutoff)          | '    
      message(6) = '                4 (conj grad)   '
      call write_fatal(6)
    end if
    if (poisson_solver /= conf%periodic_dim .and. poisson_solver /= CG) then
      write(message(1), '(a,i1,a)')'The System is periodic in ',conf%periodic_dim ,' dimension(s),'
      write(message(2), '(a,i1,a)')'but Poisson Solver is set for ',poisson_solver,' dimensions.'
      message(3) =                 'You know what you are doing, right?'
      call write_warning(3)
    end if
#else
    poisson_solver = CG
#endif

    call poisson3D_init(m)
  end if
end subroutine poisson_init

subroutine poisson_end()
  call push_sub('poisson_end')

  select case(poisson_solver)
#ifdef HAVE_FFT
  case(FFT_SPH,FFT_CYL,FFT_PLA,FFT_NOCUT)
    call dcf_free(fft_cf)
    deallocate(fft_coulb_FS); nullify(fft_coulb_FS)
#endif
  end select

  call pop_sub()
  return
end subroutine poisson_end

subroutine poisson_solve(m, f_der, pot, rho)
  type(mesh_type),  intent(in)    :: m
  type(f_der_type), intent(in)    :: f_der
  FLOAT,            intent(inout) :: pot(:)  ! pot(m%np)
  FLOAT,            intent(in)    :: rho(:)  ! rho(m%np)

  call push_sub('poisson_solve')

  ASSERT(poisson_solver.ne.-99)

  select case(poisson_solver)
  case(-1)
    call poisson1d_solve(m, pot, rho)
  case(-2)
    call poisson2d_solve(m, pot, rho)
  case(CG)
    call poisson_cg(m, f_der%der_discr, pot, rho)
#ifdef HAVE_FFT
  case(FFT_SPH,FFT_CYL,FFT_PLA,FFT_NOCUT)
    call poisson_fft(m, pot, rho)
#endif
  end select

  call pop_sub()
end subroutine poisson_solve

#include "poisson1D.F90"
#include "poisson2D.F90"
#include "poisson3D.F90"

end module poisson
