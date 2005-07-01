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
!!
!! $Id$

#include "global.h"

module poisson
  use global
  use messages
  use lib_oct_parser
  use units
  use geometry
  use mesh
#ifdef HAVE_FFT
  use fft
  use cube_function
#endif
  use functions
  use math
  use poisson_corrections
  use poisson_cg
  use grid

  implicit none

  private
  public :: poisson_init,  &
            dpoisson_solve, zpoisson_solve, &
            poisson_end

  integer :: poisson_solver = -99

#ifdef HAVE_FFT
  type(dcf) :: fft_cf
  FLOAT, pointer :: fft_coulb_FS(:,:,:)

  integer, parameter :: DIRECT_SUM_1D = -1,&
                        DIRECT_SUM_2D = -2,&
                        FFT_SPH       = 0, &
                        FFT_CYL       = 1, &
                        FFT_PLA       = 2, &
                        FFT_NOCUT     = 3, &
                        FFT_CORRECTED = 4
#endif
  integer, parameter :: CG            = 5, &
                        CG_CORRECTED  = 6



contains

  !-----------------------------------------------------------------
  subroutine poisson_init(gr)
    type(grid_type), intent(inout) :: gr

    if(poisson_solver.ne.-99) return ! already initialized
  
    call push_sub('poisson_init')
  
    select case(conf%dim)
    case(1); call init_1D()
    case(2); call init_2D()
    case(3); call init_3D()
    end select

  contains

    !-----------------------------------------------------------------
    subroutine init_1D()
      poisson_solver = -conf%dim ! internal type
      message(1) = 'Info: Using direct integration method to solve poisson equation'
      call write_info(1)
    end subroutine init_1D


    !-----------------------------------------------------------------
    subroutine init_2D()
#if defined(HAVE_FFT)
      call loct_parse_int(check_inp('PoissonSolver'), gr%sb%periodic_dim, poisson_solver)
      if( (poisson_solver .ne. FFT_SPH) .and. (poisson_solver .ne. DIRECT_SUM_2D) ) then
        write(message(1), '(a,i2,a)') "Input: '", poisson_solver, &
           "' is not a valid PoissonSolver"
        message(2) = 'PoissonSolver = -1 (direct summation)     | '
        message(3) = '                 0 (fft spherical cutoff) '
        call write_fatal(3)
      endif
      
#else
      poisson_solver = -conf%dim ! internal type
      message(1) = 'Info: Using direct integration method to solve poisson equation'
      call write_info(1)
#endif

      if(gr%m%use_curvlinear .and. (poisson_solver .ne. -conf%dim) ) then
        message(1) = 'If curvilinear coordinates are used in 2D, then the only working'
        message(2) = 'Poisson solver is -2 ("direct summation in two dimensions")'
        call write_fatal(2)
      endif
      
      call poisson2D_init(gr)
    end subroutine init_2D


    !-----------------------------------------------------------------
    subroutine init_3D()
#ifdef HAVE_FFT
      call loct_parse_int(check_inp('PoissonSolver'), gr%sb%periodic_dim, poisson_solver)
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
      if (poisson_solver /= gr%sb%periodic_dim .and. &
         poisson_solver < CG .and. &
         poisson_solver /= FFT_CORRECTED) then
        write(message(1), '(a,i1,a)')'The System is periodic in ', gr%sb%periodic_dim ,' dimension(s),'
        write(message(2), '(a,i1,a)')'but Poisson Solver is set for ',poisson_solver,' dimensions.'
        message(3) =                 'You know what you are doing, right?'
        call write_warning(3)
      end if
#else
      call loct_parse_int(check_inp('PoissonSolver'), CG, poisson_solver)
      if(poisson_solver < CG) then
        write(message(1), '(a,i2,a)') "Input: '", poisson_solver, &
           "' is not a valid PoissonSolver"
        message(2) = 'PoissonSolver = 5 (conj grad)              | '
        message(3) = '                6 (corrected conj grad) '
        message(4) = '[The code was compiled without FFT support ' 
        call write_fatal(3)
      end if
#endif
      
      if(gr%m%use_curvlinear .and. (poisson_solver .ne. CG_CORRECTED) ) then
        message(1) = 'If curvilinear coordinates are used, then the only working'
        message(2) = 'Poisson solver is 5 ("corrected conjugate gradients")'
        call write_fatal(2)
      end if

      call poisson3D_init(gr)    
    end subroutine init_3D

  end subroutine poisson_init


  !-----------------------------------------------------------------
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
    poisson_solver = -99
    
    call pop_sub()
  end subroutine poisson_end


  !-----------------------------------------------------------------
  subroutine zpoisson_solve(gr, pot, rho)
    type(grid_type),  target, intent(in)    :: gr
    CMPLX,                    intent(inout) :: pot(:)  ! pot(m%np)
    CMPLX,                    intent(in)    :: rho(:)  ! rho(m%np)
    
    FLOAT, allocatable :: aux1(:), aux2(:)
    
    call push_sub('zpoisson_solve')

    allocate(aux1(gr%m%np), aux2(gr%m%np))
    
    ! first the real part
    aux1(:) = real(rho(:))
    aux2(:) = M_ZERO
    call dpoisson_solve(gr, aux2, aux1)
    pot(:)  = aux2(:)
    
    ! now the imaginary part
    aux1(:) = aimag(rho(:))
    aux2(:) = M_ZERO
    call dpoisson_solve(gr, aux2, aux1)
    pot(:)  = pot(:) + M_zI*aux2(:)

    deallocate(aux1, aux2)

    call pop_sub()
  end subroutine zpoisson_solve


  !-----------------------------------------------------------------
  subroutine dpoisson_solve(gr, pot, rho)
    type(grid_type), target, intent(in)    :: gr
    FLOAT,                   intent(inout) :: pot(:)  ! pot(m%np)
    FLOAT,                   intent(in)    :: rho(:)  ! rho(m%np)

    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)

    call push_sub('dpoisson_solve')

    ASSERT(poisson_solver.ne.-99)

    select case(poisson_solver)
    case(-1)
      call poisson1d_solve(gr%m, pot, rho)

    case(-2)
      call poisson2d_solve(gr%m, pot, rho)

    case(CG)
      call poisson_cg1(gr%m, gr%f_der%der_discr, pot, rho)

    case(CG_CORRECTED)
      call poisson_cg2(gr%m, gr%f_der%der_discr, pot, rho)

#ifdef HAVE_FFT
    case(FFT_SPH,FFT_CYL,FFT_PLA,FFT_NOCUT)
      call poisson_fft(gr%m, pot, rho)

    case(FFT_CORRECTED)
      allocate(rho_corrected(gr%m%np), vh_correction(gr%m%np))
      call correct_rho(gr%m, maxl, rho, rho_corrected, vh_correction)
      call poisson_fft(gr%m, pot, rho_corrected, average_to_zero = .true.)
      pot = pot + vh_correction
      deallocate(rho_corrected, vh_correction)
#endif
    end select

    call pop_sub()
  end subroutine dpoisson_solve

#if defined(HAVE_FFT)
  !-----------------------------------------------------------------
  subroutine poisson_fft(m, pot, rho, average_to_zero)
    type(mesh_type), intent(IN) :: m
    FLOAT, intent(out) :: pot(:) ! pot(m%np)
    FLOAT, intent(in)  :: rho(:) ! rho(m%np)
    logical, intent(in), optional :: average_to_zero

    integer :: k
    FLOAT :: average
    
    call push_sub('poisson_fft')
  
    call dcf_alloc_RS(fft_cf)          ! allocate the cube in real space
    call dcf_alloc_FS(fft_cf)          ! allocate the cube in Fourier space
    
    call dmf2cf(m, rho, fft_cf)        ! put the density in a cube
    call dcf_RS2FS(fft_cf)             ! Fourier transform
    
    ! multiply by the FS of the Coulomb interaction
    ! this works around a bug in Intel ifort 8
    do k = 1, fft_cf%n(3)
      fft_cf%FS(:,:,k) = fft_cf%FS(:,:,k)*fft_Coulb_FS(:,:,k)
    end do
    
    call dcf_FS2RS(fft_cf)             ! Fourier transform back
    if(present(average_to_zero)) then
      if(average_to_zero) then
        average = cf_surface_average(fft_cf)
        fft_cf%RS = fft_cf%RS - average
      end if
    end if
    call dcf2mf(m, fft_cf, pot)        ! put the density in a cube

    call dcf_free_RS(fft_cf)           ! memory is no longer needed
    call dcf_free_FS(fft_cf)

    call pop_sub()
  end subroutine poisson_fft
#endif


#include "poisson1D.F90"
#include "poisson2D.F90"
#include "poisson3D.F90"
  
end module poisson
