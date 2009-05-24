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

module poisson_m
  use datasets_m
  use geometry_m
  use global_m
  use io_function_m
  use io_m
  use loct_math_m
  use loct_parser_m
  use mesh_m
  use messages_m
  use mpi_m
  use profiling_m
  use simul_box_m
  use units_m
  use math_m
  use poisson_corrections_m
  use poisson_cg_m
  use poisson_isf_m
  use poisson_fft_m
  use grid_m
  use poisson_multigrid_m
  use mesh_function_m
  use par_vec_m

  implicit none

  private
  public ::             &
    poisson_init,       &
    dpoisson_solve,     &
    zpoisson_solve,     &
    poisson_solver_is_iterative, &
    poisson_end,        &
    poisson_test

  integer, parameter :: &
    DIRECT_SUM_1D = -1, &
    DIRECT_SUM_2D = -2, &
    CG            =  5, &
    CG_CORRECTED  =  6, &
    MULTIGRID     =  7, &
    ISF           =  8

  integer :: poisson_solver = -99
  FLOAT   :: poisson_soft_coulomb_param = M_ONE
  type(mg_solver_t) :: mg
  logical :: all_nodes_default

  type hartree_t
    private
    logical                 :: increase_box
    type(grid_t), pointer   :: grid
  end type hartree_t

  type(hartree_t) :: hartree_integrator
  type(poisson_corr_t) :: corrector

contains

  !-----------------------------------------------------------------
  subroutine poisson_init(gr, geo)
    type(grid_t), intent(inout) :: gr
    type(geometry_t), intent(in) :: geo

    if(poisson_solver.ne.-99) return ! already initialized

    call push_sub('poisson.poisson_init')

    call messages_print_stress(stdout, "Hartree")

#ifdef HAVE_MPI
    !%Variable ParallelizationPoissonAllNodes
    !%Type logical
    !%Default true
    !%Section Execution::Parallelization
    !%Description
    !% When running in parallel, this variable selects whether the
    !% Poisson solver should divide the work among all nodes or only
    !% among the parallelization in domains groups.
    !%End

    call loct_parse_logical(datasets_check('ParallelizationPoissonAllNodes'), .true., all_nodes_default)
#endif

    select case(gr%mesh%sb%dim)
    case(1); call init_1D()
    case(2); call init_2D()
    case(3); call init_3D()
    end select

    call messages_print_stress(stdout)

    call pop_sub()
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
    !% Conjugate gradients
    !%Option cg_corrected 6
    !% Corrected conjugate gradients
    !%Option multigrid 7
    !% Multigrid method
    !%Option isf 8
    !% Interpolating Scaling Functions poisson solver.
    !%End

    !------------------------------------------------------------------
    subroutine init_1D()
      integer :: default_solver

      call push_sub('poisson.init_1D')

      if(gr%sb%periodic_dim.eq.0) then
        default_solver = FFT_SPH
      else
        default_solver = FFT_NOCUT
      end if
      call loct_parse_int(datasets_check('PoissonSolver'), default_solver, poisson_solver)

      select case(gr%sb%periodic_dim)
      case(0)
        if( (poisson_solver.ne.FFT_SPH)       .and. &
            (poisson_solver.ne.DIRECT_SUM_1D)) call input_error('PoissonSolver')
      case(1)
        if( (poisson_solver.ne.FFT_NOCUT) ) call input_error('PoissonSolver')
      end select

      if(gr%mesh%use_curvilinear.and.poisson_solver.ne.DIRECT_SUM_1D) then
        message(1) = 'If curvilinear coordinates are used in 1D, then the only working'
        message(2) = 'Poisson solver is -1 ("direct summation in one dimension").'
        call write_fatal(2)
      end if

      call messages_print_var_option(stdout, "PoissonSolver", poisson_solver)
      call poisson1d_init(gr)

      call pop_sub()
    end subroutine init_1D


    !-----------------------------------------------------------------
    subroutine init_2D()
      integer :: default_solver
      call push_sub('poisson.init_2D')

      if (gr%sb%periodic_dim > 0) then 
        default_solver = gr%sb%periodic_dim
      else
        default_solver = FFT_SPH
      end if

      call loct_parse_int(datasets_check('PoissonSolver'), gr%sb%periodic_dim, poisson_solver)
      if( (poisson_solver .ne. FFT_SPH)         .and. &
          (poisson_solver .ne. DIRECT_SUM_2D)   .and. &
          (poisson_solver .ne. 1)               .and. &
          (poisson_solver .ne. 2)               .and. &
          (poisson_solver .ne. 3)                       ) then
        call input_error('PoissonSolver')
      end if

      ! In 2D, periodic in two dimensions means no cut-off at all.
      if(poisson_solver .eq. 2) poisson_solver = 3

      if(gr%mesh%use_curvilinear .and. (poisson_solver .ne. -gr%mesh%sb%dim) ) then
        message(1) = 'If curvilinear coordinates are used in 2D, then the only working'
        message(2) = 'Poisson solver is -2 ("direct summation in two dimensions")'
        call write_fatal(2)
      end if

      call messages_print_var_option(stdout, "PoissonSolver", poisson_solver)
      call poisson2D_init(gr)

      call pop_sub()
    end subroutine init_2D


    !-----------------------------------------------------------------
    subroutine init_3D()
      integer :: default_solver

      call push_sub('poisson.init_3D')

#ifndef SINGLE_PRECISION
      default_solver = ISF
#else
      default_solver = FFT_SPH
#endif

      if (gr%mesh%use_curvilinear) default_solver = CG_CORRECTED
      if (gr%sb%periodic_dim > 0) default_solver = gr%sb%periodic_dim
      
      call loct_parse_int(datasets_check('PoissonSolver'), default_solver, poisson_solver)
      if(poisson_solver < FFT_SPH .or. poisson_solver > ISF ) then
        call input_error('PoissonSolver')
      end if

      if(gr%sb%periodic_dim > 0 .and. &
           poisson_solver /= gr%sb%periodic_dim .and. &
           poisson_solver < CG .and. &
           poisson_solver /= FFT_CORRECTED .and. poisson_solver /= FFT_CYL) then
        write(message(1), '(a,i1,a)')'The system is periodic in ', gr%sb%periodic_dim ,' dimension(s),'
        write(message(2), '(a,i1,a)')'but Poisson solver is set for ',poisson_solver,' dimensions.'
        message(3) =                 'You know what you are doing, right?'
        call write_warning(3)
      end if

      if(gr%mesh%use_curvilinear .and. (poisson_solver.ne.CG_CORRECTED)) then
        message(1) = 'If curvilinear coordinates are used, then the only working'
        message(2) = 'Poisson solver is cg_corrected'
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
      call poisson3D_init(gr, geo)

      call pop_sub()
    end subroutine init_3D

  end subroutine poisson_init


  !-----------------------------------------------------------------
  subroutine poisson_end()
    call push_sub('poisson.poisson_end')

    select case(poisson_solver)

    case(FFT_SPH,FFT_CYL,FFT_PLA,FFT_NOCUT)
      call poisson_fft_end()
    case(FFT_CORRECTED)
      call poisson_fft_end()
      call poisson_corrections_end(corrector)

    case(CG_CORRECTED, CG)
      call poisson_cg_end()
      call poisson_corrections_end(corrector)

    case(MULTIGRID)
      call poisson_multigrid_end(mg)
    case(ISF)
      call poisson_isf_end()

    end select
    poisson_solver = -99

    if(hartree_integrator%increase_box) then
      SAFE_DEALLOCATE_P(hartree_integrator%grid)
    end if

    call pop_sub()
  end subroutine poisson_end


  !-----------------------------------------------------------------
  subroutine zpoisson_solve(gr, pot, rho, all_nodes)
    type(grid_t),         intent(inout) :: gr
    CMPLX,                intent(inout) :: pot(:)  ! pot(m%np)
    CMPLX,                intent(in)    :: rho(:)  ! rho(m%np)
    logical, optional,    intent(in)    :: all_nodes

    FLOAT, allocatable :: aux1(:), aux2(:)

    logical :: all_nodes_value

    call push_sub('poisson.zpoisson_solve')

    if(present(all_nodes)) then
      all_nodes_value = all_nodes
    else
      all_nodes_value = all_nodes_default
    end if

    SAFE_ALLOCATE(aux1(1:gr%mesh%np))
    SAFE_ALLOCATE(aux2(1:gr%mesh%np))

    ! first the real part
    aux1(1:gr%mesh%np) = real(rho(1:gr%mesh%np))
    aux2(1:gr%mesh%np) = real(pot(1:gr%mesh%np))
    call dpoisson_solve(gr, aux2, aux1, all_nodes=all_nodes_value)
    pot(1:gr%mesh%np)  = aux2(1:gr%mesh%np)

    ! now the imaginary part
    aux1(1:gr%mesh%np) = aimag(rho(1:gr%mesh%np))
    aux2(1:gr%mesh%np) = aimag(pot(1:gr%mesh%np))
    call dpoisson_solve(gr, aux2, aux1, all_nodes=all_nodes_value)
    pot(1:gr%mesh%np) = pot(1:gr%mesh%np) + M_zI*aux2(1:gr%mesh%np)

    SAFE_DEALLOCATE_A(aux1)
    SAFE_DEALLOCATE_A(aux2)

    call pop_sub()
  end subroutine zpoisson_solve


  !-----------------------------------------------------------------
  subroutine dpoisson_solve(gr, pot, rho, all_nodes)
    type(grid_t),         intent(inout) :: gr
    FLOAT,                intent(inout) :: pot(:)    ! pot(m%np)
    FLOAT,                intent(in)    :: rho(:)    ! rho(m%np)
    logical, optional,    intent(in)    :: all_nodes ! Is the poisson solver allowed to utilise
                                                     ! all nodes or only the domain nodes for
                                                     ! its calculations? (Defaults to .true.)

    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)
    FLOAT, allocatable :: rhop(:), potp(:)

    logical               :: all_nodes_value
    type(profile_t), save :: prof
    
    call profiling_in(prof, 'POISSON_SOLVE')
    call push_sub('poisson.dpoisson_solve')

    ! Check optional argument and set to default if necessary.
    if(present(all_nodes)) then
      all_nodes_value = all_nodes
    else
      all_nodes_value = all_nodes_default
    end if

    ASSERT(poisson_solver.ne.-99)

    select case(poisson_solver)
    case(DIRECT_SUM_1D)
      call poisson1d_solve(gr%mesh, pot, rho)

    case(DIRECT_SUM_2D)
      call poisson2d_solve(gr%mesh, pot, rho)

    case(CG)
      call poisson_cg1(gr%mesh, corrector, gr%der, pot, rho)

    case(CG_CORRECTED)
      if(hartree_integrator%increase_box) then
        SAFE_ALLOCATE(potp(1:hartree_integrator%grid%mesh%np))
        SAFE_ALLOCATE(rhop(1:hartree_integrator%grid%mesh%np))
        SAFE_ALLOCATE(rho_corrected(1:gr%mesh%np))
        SAFE_ALLOCATE(vh_correction(1:gr%mesh%np))

        potp = M_ZERO; rhop = M_ZERO; rho_corrected = M_ZERO; vh_correction = M_ZERO
        call correct_rho(corrector, gr%mesh, rho, rho_corrected, vh_correction)

        pot(1:gr%mesh%np) = pot(1:gr%mesh%np) - vh_correction(1:gr%mesh%np)

        call dmf_interpolate(gr%mesh, hartree_integrator%grid%mesh, full_interpolation = .false., u = rho_corrected, f = rhop)
        call dmf_interpolate(gr%mesh, hartree_integrator%grid%mesh, full_interpolation = .false., u = pot, f = potp)

        call poisson_cg2(hartree_integrator%grid%mesh, hartree_integrator%grid%der, potp, rhop)

        call dmf_interpolate(hartree_integrator%grid%mesh, gr%mesh, full_interpolation = .false., u = potp, f = pot)
        pot(1:gr%mesh%np) = pot(1:gr%mesh%np) + vh_correction(1:gr%mesh%np)

        SAFE_DEALLOCATE_A(rho_corrected)
        SAFE_DEALLOCATE_A(vh_correction)
        SAFE_DEALLOCATE_A(potp)
        SAFE_DEALLOCATE_A(rhop)
      else
        SAFE_ALLOCATE(rho_corrected(1:gr%mesh%np))
        SAFE_ALLOCATE(vh_correction(1:gr%mesh%np))

        call correct_rho(corrector, gr%mesh, rho, rho_corrected, vh_correction)

        pot(1:gr%mesh%np) = pot(1:gr%mesh%np) - vh_correction(1:gr%mesh%np)
        call poisson_cg2(gr%mesh, gr%der, pot, rho_corrected)
        pot(1:gr%mesh%np) = pot(1:gr%mesh%np) + vh_correction(1:gr%mesh%np)

        SAFE_DEALLOCATE_A(rho_corrected)
        SAFE_DEALLOCATE_A(vh_correction)
      end if

    case(MULTIGRID)
      call poisson_multigrid_solver(mg, gr, pot, rho)

    case(FFT_SPH,FFT_CYL,FFT_PLA,FFT_NOCUT)
      call poisson_fft(gr%mesh, pot, rho)

    case(FFT_CORRECTED)
      SAFE_ALLOCATE(rho_corrected(1:gr%mesh%np))
      SAFE_ALLOCATE(vh_correction(1:gr%mesh%np))

      call correct_rho(corrector, gr%mesh, rho, rho_corrected, vh_correction)
      call poisson_fft(gr%mesh, pot, rho_corrected, average_to_zero = .true.)

      pot(1:gr%mesh%np) = pot(1:gr%mesh%np) + vh_correction(1:gr%mesh%np)
      SAFE_DEALLOCATE_A(rho_corrected)
      SAFE_DEALLOCATE_A(vh_correction)

    case(ISF)
      call poisson_isf_solve(gr%mesh, pot, rho, all_nodes_value)
      
    end select

    call pop_sub()
    call profiling_out(prof)
  end subroutine dpoisson_solve



  !-----------------------------------------------------------------
  ! This routine checks the Hartree solver selected in the input
  ! file by calculating numerically and analytically the Hartree
  ! potential originated by a Gaussian distribution of charge.
  ! This only makes sense for finite systems.
  subroutine poisson_test(gr)
    type(grid_t), intent(inout) :: gr

    FLOAT, allocatable :: rho(:), vh(:), vh_exact(:), rhop(:), x(:, :)
    FLOAT :: alpha, beta, r, delta, norm, rnd, ralpha, range
    integer :: i, k, ierr, iunit, n, n_gaussians

    call push_sub('poisson.poisson_test')

    if(calc_dim.eq.1) then
      write(message(1),'(a)') 'The Hartree integrator test is not implemented for the one dimensional case.'
      call write_warning(1)
      call pop_sub(); return
    endif

    n_gaussians = 4

    SAFE_ALLOCATE(     rho(1:gr%mesh%np))
    SAFE_ALLOCATE(    rhop(1:gr%mesh%np))
    SAFE_ALLOCATE(      vh(1:gr%mesh%np))
    SAFE_ALLOCATE(vh_exact(1:gr%mesh%np))
    SAFE_ALLOCATE(x(1:gr%mesh%sb%dim, 1:n_gaussians))

    rho = M_ZERO; vh = M_ZERO; vh_exact = M_ZERO

    alpha = CNST(4.0) * gr%mesh%h(1)
    beta = M_ONE / ( alpha**calc_dim * sqrt(M_PI)**calc_dim )

    write(message(1), '(a)') 'Building the Gaussian distribution of charge...'
    call write_info(1)

    range = CNST(8.0)

    rho = M_ZERO
    do n = 1, n_gaussians
      norm = M_ZERO
      do while(abs(norm-M_ONE)> CNST(1.0e-4))
        do k = 1, gr%mesh%sb%dim
          call random_number(rnd)
          x(k, n) = range * rnd 
        end do
        r = sqrt(sum(x(:, n)*x(:,n)))
        do i = 1, gr%mesh%np
          call mesh_r(gr%mesh, i, r, a = x(:, n))
          rhop(i) = beta*exp(-(r/alpha)**2)
        end do
        norm = dmf_integrate(gr%mesh, rhop)
      end do
      rhop = (-1)**n * rhop
      rho = rho + rhop
    end do
    write(message(1), '(a,f14.6)') 'Total charge of the Gaussian distribution', dmf_integrate(gr%mesh, rho)
    call write_info(1)

    ! This builds analytically its potential
    vh_exact = M_ZERO
    do n = 1, n_gaussians
      do i = 1, gr%mesh%np
        call mesh_r(gr%mesh, i, r, a = x(:, n))
        select case(calc_dim)
        case(3)
          if(r > r_small) then
            vh_exact(i) = vh_exact(i) + (-1)**n * loct_erf(r/alpha)/r
          else
            vh_exact(i) = vh_exact(i) + (-1)**n * (M_TWO/sqrt(M_PI))/alpha
          end if
        case(2)
          ralpha = r**2/(M_TWO*alpha**2)
          if(ralpha < CNST(100.0)) then
            vh_exact(i) = vh_exact(i) + (-1)**n * beta * (M_PI)**(M_THREE*M_HALF) * alpha * exp(-r**2/(M_TWO*alpha**2)) * &
              loct_bessel_in(0, r**2/(M_TWO*alpha**2))
          else
            vh_exact(i) = vh_exact(i) + (-1)**n * beta * (M_PI)**(M_THREE*M_HALF) * alpha * &
                          (M_ONE/sqrt(M_TWO*M_PI*ralpha)) 
          end if
        end select
      end do
    end do

    ! This calculates the numerical potential
    call dpoisson_solve(gr, vh, rho)

    ! And this compares.
    delta = dmf_nrm2(gr%mesh, vh-vh_exact)

    ! Output
    iunit = io_open("hartree_results", action='write')
    write(iunit, '(a,f10.2)' ) 'Hartree test = ', delta
    call io_close(iunit)
    call doutput_function (io_function_fill_how('AxisX'), ".", "poisson_test_rho", gr%mesh, gr%sb, rho, M_ONE, ierr)
    call doutput_function (io_function_fill_how('AxisX'), ".", "poisson_test_exact", gr%mesh, gr%sb, vh_exact, M_ONE, ierr)
    call doutput_function (io_function_fill_how('AxisX'), ".", "poisson_test_numerical", gr%mesh, gr%sb, vh, M_ONE, ierr)
    call doutput_function (io_function_fill_how('AxisY'), ".", "poisson_test_rho", gr%mesh, gr%sb, rho, M_ONE, ierr)
    call doutput_function (io_function_fill_how('AxisY'), ".", "poisson_test_exact", gr%mesh, gr%sb, vh_exact, M_ONE, ierr)
    call doutput_function (io_function_fill_how('AxisY'), ".", "poisson_test_numerical", gr%mesh, gr%sb, vh, M_ONE, ierr)

    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(rhop)
    SAFE_DEALLOCATE_A(vh)
    SAFE_DEALLOCATE_A(vh_exact)
    SAFE_DEALLOCATE_A(x)
    call pop_sub()
  end subroutine poisson_test

  logical pure function poisson_solver_is_iterative() result(iterative)

    iterative = poisson_solver == CG .or. poisson_solver == CG_CORRECTED .or. poisson_solver == MULTIGRID
    
  end function poisson_solver_is_iterative

#include "solver_1D.F90"
#include "solver_2D.F90"
#include "solver_3D.F90"

end module poisson_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
