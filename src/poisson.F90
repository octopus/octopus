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
  use units
  use geometry
  use mesh
#ifdef HAVE_FFT
  use fft
  use cube_function
#endif
!!!!WARNING
  use mesh_function
!!!!END OF WARNING
  use functions
  use math
  use poisson_corrections
  use poisson_cg

  implicit none

  private
  public :: poisson_init,  &
            poisson_solve, &
            poisson_end

  integer :: poisson_solver = -99

#ifdef HAVE_FFT
  type(dcf) :: fft_cf
  FLOAT, pointer :: fft_coulb_FS(:,:,:)

  integer, parameter :: FFT_SPH       = 0, &
                        FFT_CYL       = 1, &
                        FFT_PLA       = 2, &
                        FFT_NOCUT     = 3, &
                        FFT_CORRECTED = 4
#endif

  integer, parameter :: CG            = 5, &
                        CG_CORRECTED  = 6



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
    call loct_parse_int("PoissonSolver", conf%periodic_dim, poisson_solver)
    if(poisson_solver < FFT_SPH .or. poisson_solver > CG_CORRECTED ) then
      write(message(1), '(a,i2,a)') "Input: '", poisson_solver, &
           "' is not a valid PoissonSolver"
      message(2) = 'PoissonSolver = 0 (fft spherical cutoff)   | '
      message(3) = '                1 (fft cylindrical cutoff) | '
      message(4) = '                2 (fft planar cutoff)      | '
      message(5) = '                3 (fft no cutoff)          | '
      message(6) = '                4 (fft + corrections)      | '
      message(7) = '                5 (conj grad)              | '
      message(8) = '                6 (corrected conj grad) '
      call write_fatal(8)
    end if
    if (poisson_solver /= conf%periodic_dim .and. &
        poisson_solver < CG .and. &
        poisson_solver /= FFT_CORRECTED) then
      write(message(1), '(a,i1,a)')'The System is periodic in ',conf%periodic_dim ,' dimension(s),'
      write(message(2), '(a,i1,a)')'but Poisson Solver is set for ',poisson_solver,' dimensions.'
      message(3) =                 'You know what you are doing, right?'
      call write_warning(3)
    end if
#else
    call loct_parse_int('PoissonSolver', CG, poisson_solver)
    if(poisson_solver < CG) then
      write(message(1), '(a,i2,a)') "Input: '", poisson_solver, &
           "' is not a valid PoissonSolver"
      message(2) = 'PoissonSolver = 5 (conj grad)              | '
      message(3) = '                6 (corrected conj grad) '
      message(4) = '[The code was compiled without FFT support ' 
      call write_fatal(3)
    endif
#endif

    if(m%use_curvlinear .and. (poisson_solver .ne. CG_CORRECTED) ) then
      message(1) = 'If curvilinear coordinates are used, then the only working'
      message(2) = 'Poisson solver is 5 ("corrected conjugate gradients")'
      call write_fatal(2)
    endif
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
  case(CG_CORRECTED)
    call poisson_cg2_end()
  end select

  call pop_sub()
  return
end subroutine poisson_end

subroutine poisson_solve(m, f_der, pot, rho)
  type(mesh_type), target,  intent(in)    :: m
  type(f_der_type), target, intent(in)    :: f_der
  FLOAT,            intent(inout) :: pot(:)  ! pot(m%np)
  FLOAT,            intent(in)    :: rho(:)  ! rho(m%np)

    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)

  call push_sub('poisson_solve')

  ASSERT(poisson_solver.ne.-99)

  select case(poisson_solver)
  case(-1)
    call poisson1d_solve(m, pot, rho)
  case(-2)
    call poisson2d_solve(m, pot, rho)
  case(CG)
    call poisson_cg1(m, f_der%der_discr, pot, rho)
  case(CG_CORRECTED)
    call poisson_cg2(m, f_der%der_discr, pot, rho)
#ifdef HAVE_FFT
  case(FFT_SPH,FFT_CYL,FFT_PLA,FFT_NOCUT)
    call poisson_fft(m, pot, rho)
  case(FFT_CORRECTED)
    allocate(rho_corrected(m%np), vh_correction(m%np))
    call correct_rho(m, maxl, rho, rho_corrected, vh_correction)
    call poisson_fft(m, pot, rho_corrected, average_to_zero = .true.)
    pot = pot + vh_correction
    deallocate(rho_corrected, vh_correction)
#endif
  end select

  call pop_sub()
end subroutine poisson_solve

#include "poisson1D.F90"
#include "poisson2D.F90"
#include "poisson3D.F90"

end module poisson
