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

module poisson_m
  use global_m
  use messages_m
  use profiling_m
  use lib_oct_parser_m
  use units_m
  use geometry_m
  use mesh_m
#ifdef HAVE_FFT
  use fft_m
  use cube_function_m
#endif
  use functions_m
  use math_m
  use poisson_corrections_m
  use poisson_cg_m
  use grid_m
  use output_m
  use poisson_multigrid_m
  use mesh_function_m
  implicit none

  private
  public ::             &
    poisson_init,       &
    dpoisson_solve,     &
    zpoisson_solve,     &
    poisson_end

  integer, parameter :: &
    CG            =  5, &
    CG_CORRECTED  =  6, &
    MULTIGRILLA   =  7

#ifdef HAVE_FFT
  integer, parameter :: &
    DIRECT_SUM_1D = -1, &
    DIRECT_SUM_2D = -2, &
    FFT_SPH       =  0, &
    FFT_CYL       =  1, &
    FFT_PLA       =  2, &
    FFT_NOCUT     =  3, &
    FFT_CORRECTED =  4

  type(dcf) :: fft_cf
  FLOAT, pointer :: fft_coulb_FS(:,:,:)
#endif

  integer :: poisson_solver = -99


contains

  !-----------------------------------------------------------------
  subroutine poisson_init(gr)
    type(grid_t), intent(inout) :: gr

    if(poisson_solver.ne.-99) return ! already initialized

    call push_sub('poisson.poisson_init')

    select case(NDIM)
    case(1); call init_1D()
    case(2); call init_2D()
    case(3); call init_3D()
    end select

  contains

    !%Variable PoissonSolver
    !%Type integer
    !%Default fft
    !%Section Hamiltonian::Poisson
    !%Description
    !% Defines which method to use in order to solve the Poisson equation.
    !% The default for 1D and 2D is the direct evaluation of the Hartree potential.
    !%Option direct2D -2
    !% Direct evaluation of the Hartree potential (in 2D)
    !%Option direct1D -1
    !% Direct evaluation of the Hartree potential (in 1D)
    !%Option fft 0
    !% FFTs using spherical cutoff (in 2D or 3D; uses FFTW)
    !%Option fft_cyl 1
    !% FFTs using cylindrical cutoff (in 3D; uses FFTW)
    !%Option fft_pla 2
    !% FFTs using planar cutoff (in 3D; uses FFTW)
    !%Option fft_nocut 3
    !% FFTs without using a cutoff (in 3D; uses FFTW)
    !%Option fft_corrected 4
    !% FFTs + corrections
    !%Option cg 5
    !% Conjugated gradients
    !%Option cg_corrected 6
    !% Corrected conjugated gradients
    !%Option multigrid 7
    !% Multigrid method
    !%End

    !-----------------------------------------------------------------
    subroutine init_1D()
      poisson_solver = -NDIM ! internal type

      call messages_print_var_option(stdout, "PoissonSolver", poisson_solver)
    end subroutine init_1D


    !-----------------------------------------------------------------
    subroutine init_2D()
#if defined(HAVE_FFT)
      call loct_parse_int(check_inp('PoissonSolver'), gr%sb%periodic_dim, poisson_solver)
      if( (poisson_solver .ne. FFT_SPH) .and. (poisson_solver .ne. DIRECT_SUM_2D) ) then
        call input_error('PoissonSolver')
      end if

#else
      poisson_solver = -NDIM ! internal type
#endif

      if(gr%m%use_curvlinear .and. (poisson_solver .ne. -NDIM) ) then
        message(1) = 'If curvilinear coordinates are used in 2D, then the only working'
        message(2) = 'Poisson solver is -2 ("direct summation in two dimensions")'
        call write_fatal(2)
      end if

      call messages_print_var_option(stdout, "PoissonSolver", poisson_solver)
      call poisson2D_init(gr)
    end subroutine init_2D


    !-----------------------------------------------------------------
    subroutine init_3D()
#ifdef HAVE_FFT
      call loct_parse_int(check_inp('PoissonSolver'), gr%sb%periodic_dim, poisson_solver)
      if(poisson_solver < FFT_SPH .or. poisson_solver > MULTIGRILLA ) then
        call input_error('PoissonSolver')
      end if

      if(poisson_solver /= gr%sb%periodic_dim .and. &
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
        call input_error('PoissonSolver')
      end if
#endif

      if(gr%m%use_curvlinear .and. (poisson_solver.ne.CG_CORRECTED) .and. (poisson_solver.ne.MULTIGRILLA) ) then
        message(1) = 'If curvilinear coordinates are used, then the only working'
        message(2) = 'Poisson solvers are cg_corrected ("corrected conjugate gradients") and'
        message(3) = 'multigrid.'
        call write_fatal(3)
      end if

      if(gr%m%parallel_in_domains .and. (poisson_solver.ne.CG).and.(poisson_solver.ne.CG_CORRECTED)) then
        message(1) = 'When running in parallel in domains, you can only use'
        message(2) = 'PoissonSolver = cg | cg_corrected'
        call write_fatal(2)
      end if

      if( (gr%sb%box_shape .eq. MINIMUM) .and. (poisson_solver .eq. CG_CORRECTED) ) then
        message(1) = 'When using the "minimum" box shape and the "cg_corrected"'
        message(2) = 'Poisson solver, we have observed "sometimes" some non-'
        message(3) = 'negligible error. You may want to check that the "fft" or "cg"'
        message(4) = 'solver are providing, in your case, the same results.'
        call write_warning(4)
      end if

      call messages_print_var_option(stdout, "PoissonSolver", poisson_solver)
      call poisson3D_init(gr)
    end subroutine init_3D

  end subroutine poisson_init


  !-----------------------------------------------------------------
  subroutine poisson_end()
    call push_sub('poisson.poisson_end')

    select case(poisson_solver)
#ifdef HAVE_FFT
    case(FFT_SPH,FFT_CYL,FFT_PLA,FFT_NOCUT)
      call dcf_free(fft_cf)
      deallocate(fft_coulb_FS); nullify(fft_coulb_FS)
#endif
    case(CG_CORRECTED)
      call poisson_cg2_end()
    case(MULTIGRILLA)
      call poisson_multigrid_end()
    end select
    poisson_solver = -99

    call pop_sub()
  end subroutine poisson_end


  !-----------------------------------------------------------------
  subroutine zpoisson_solve(gr, pot, rho)
    type(grid_t),  target, intent(inout) :: gr
    CMPLX,                    intent(inout) :: pot(:)  ! pot(m%np)
    CMPLX,                    intent(in)    :: rho(:)  ! rho(m%np)

    FLOAT, allocatable :: aux1(:), aux2(:)

    call push_sub('poisson.zpoisson_solve')

    ALLOCATE(aux1(gr%m%np), gr%m%np)
    ALLOCATE(aux2(gr%m%np), gr%m%np)

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
    type(grid_t), target, intent(inout) :: gr
    FLOAT,                   intent(inout) :: pot(:)  ! pot(m%np)
    FLOAT,                   intent(in)    :: rho(:)  ! rho(m%np)

    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)

    call profiling_in(C_PROFILING_POISSON_SOLVE)
    call push_sub('poisson.dpoisson_solve')

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

    case(MULTIGRILLA)
      call poisson_multigrid_solver(gr, pot, rho)

#ifdef HAVE_FFT
    case(FFT_SPH,FFT_CYL,FFT_PLA,FFT_NOCUT)
      call poisson_fft(gr%m, pot, rho)

    case(FFT_CORRECTED)
      ALLOCATE(rho_corrected(gr%m%np), gr%m%np)
      ALLOCATE(vh_correction(gr%m%np), gr%m%np)

      call correct_rho(gr%m, maxl, rho, rho_corrected, vh_correction)
      call poisson_fft(gr%m, pot, rho_corrected, average_to_zero = .true.)

      pot = pot + vh_correction
      deallocate(rho_corrected, vh_correction)
#endif
    end select

    call pop_sub()
    call profiling_out(C_PROFILING_POISSON_SOLVE)
  end subroutine dpoisson_solve

#if defined(HAVE_FFT)
  !-----------------------------------------------------------------
  subroutine poisson_fft(m, pot, rho, average_to_zero)
    type(mesh_t), intent(in) :: m
    FLOAT, intent(out) :: pot(:) ! pot(m%np)
    FLOAT, intent(in)  :: rho(:) ! rho(m%np)
    logical, intent(in), optional :: average_to_zero

    integer :: k
    FLOAT :: average

    call push_sub('poisson.poisson_fft')

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
      if(average_to_zero) average = cf_surface_average(fft_cf)
    end if
    call dcf2mf(m, fft_cf, pot)        ! put the density in a cube
    if(present(average_to_zero)) then
      if(average_to_zero) pot = pot - average
    end if

    call dcf_free_RS(fft_cf)           ! memory is no longer needed
    call dcf_free_FS(fft_cf)

    call pop_sub()
  end subroutine poisson_fft
#endif

!!$  !-----------------------------------------------------------------
!!$  ! This routine checks the Hartree solver selected in the input
!!$  ! file by calculating numerically and analytically the Hartree
!!$  ! potential originated by a Gaussian distribution of charge.
!!$  ! This only makes sense for finite systems.
!!$  subroutine poisson_test(gr)
!!$    type(grid_t), intent(inout) :: gr
!!$
!!$    FLOAT, allocatable :: rho(:), vh(:), vh_exact(:)
!!$    FLOAT :: alpha, beta, r
!!$    integer :: i, ierr
!!$
!!$    call push_sub('poisson.poisson_test')
!!$
!!$    if(calc_dim.eq.1) then
!!$      write(message(1),'(a)') 'The Hartree integrator test is not implemented for the one dimensional case.'
!!$      call write_warning(1)
!!$      call pop_sub()
!!$      return
!!$    endif
!!$
!!$!   /* 
!!$    ALLOCATE(     rho(NP), NP)
!!$    ALLOCATE(      vh(NP), NP)
!!$    ALLOCATE(vh_exact(NP), NP)
!!$!   */
!!$    rho = M_ZERO; vh = M_ZERO; vh_exact = M_ZERO
!!$
!!$    ! This builds a normalized Gaussian charge
!!$    alpha = CNST(5.0) * gr%m%h(1)
!!$    beta = M_ONE / ( alpha**calc_dim * sqrt(M_PI)**calc_dim )
!!$    do i = 1, NP
!!$      call mesh_r(gr%m, i, r)
!!$      rho(i) = beta*exp(-(r/alpha)**2)
!!$    end do
!!$    write(message(1), '(a,f14.6)') 'Total charge of the Gaussian distribution', dmf_integrate(gr%m, rho)
!!$
!!$    ! This builds analytically its potential
!!$    do i = 1, NP
!!$      call mesh_r(gr%m, i, r)
!!$      select case(calc_dim)
!!$      case(3)
!!$        if(r > r_small) then
!!$          vh_exact(i) = loct_erf(r/alpha)/r
!!$        else
!!$          vh_exact(i) = (M_TWO/sqrt(M_PI))/alpha
!!$        end if
!!$      case(2)
!!$        vh_exact(i) = beta * (M_PI)**(M_THREE*M_HALF) * alpha * exp(-r**2/(M_TWO*alpha**2)) * &
!!$          loct_bessel_in(0, r**2/(M_TWO*alpha**2))
!!$      end select
!!$    end do
!!$
!!$    ! This calculates the numerical potential
!!$    call dpoisson_solve(gr, vh, rho)
!!$
!!$    ! And this compares.
!!$    write(message(2), '(a,f14.6)') 'Difference between exact and numerical result:', &
!!$      dmf_integrate(gr%m, vh-vh_exact)
!!$    ! Output
!!$    call write_info(2)
!!$    call doutput_function (output_fill_how('AxisX'), ".", "poisson_test_rho", gr%m, gr%sb, rho, M_ONE, ierr)
!!$    call doutput_function (output_fill_how('AxisX'), ".", "poisson_test_exact", gr%m, gr%sb, vh_exact, M_ONE, ierr)
!!$    call doutput_function (output_fill_how('AxisX'), ".", "poisson_test_numerical", gr%m, gr%sb, vh, M_ONE, ierr)
!!$
!!$    deallocate(rho, vh, vh_exact)
!!$    call pop_sub()
!!$  end subroutine poisson_test


#include "poisson1D.F90"
#include "poisson2D.F90"
#include "poisson3D.F90"

end module poisson_m
