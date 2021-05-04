!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio,
!! G. Bertsch, M. Oliveira
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

module poisson_oct_m
  use batch_oct_m
  use box_minimum_oct_m
  use comm_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use derivatives_oct_m
  use fft_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use index_oct_m
  use io_oct_m
  use io_function_oct_m
  use lattice_vectors_oct_m
  use loct_math_oct_m
  use mesh_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
#ifdef HAVE_OPENMP
  use omp_lib
#endif
  use par_vec_oct_m
  use parser_oct_m
  use partition_oct_m
  use photon_mode_oct_m
  use poisson_cg_oct_m
  use poisson_corrections_oct_m
  use poisson_isf_oct_m
  use poisson_fft_oct_m
  use poisson_fmm_oct_m
  use poisson_psolver_oct_m
  use poisson_multigrid_oct_m
  use poisson_no_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use submesh_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

#ifdef HAVE_POKE
  use poke
#endif

  implicit none

  private
  public ::                      &
    poisson_t,                   &
    poisson_fmm_t,               &
    poisson_get_solver,          &
    poisson_init,                &
    poisson_init_sm,             &
    dpoisson_solve,              &
    zpoisson_solve,              &
    dpoisson_solve_sm,           &
    zpoisson_solve_sm,           &
    poisson_solve_batch,         &
    poisson_solver_is_iterative, &
    poisson_solver_has_free_bc,  &
    poisson_end,                 &
    poisson_test,                &
    poisson_is_multigrid,        &
    poisson_slave_work,          &
    poisson_async_init,          &
    poisson_async_end,           &
    dpoisson_solve_start,        &
    dpoisson_solve_finish,       &
    zpoisson_solve_start,        &
    zpoisson_solve_finish,       &
    poisson_build_kernel,        &
    poisson_is_async

  integer, public, parameter ::         &
    POISSON_DIRECT_SUM    = -1,         &
    POISSON_FMM           = -4,         &
    POISSON_FFT           =  0,         &
    POISSON_CG            =  5,         &
    POISSON_CG_CORRECTED  =  6,         &
    POISSON_MULTIGRID     =  7,         &
    POISSON_ISF           =  8,         &
    POISSON_PSOLVER       = 10,         &
    POISSON_POKE          = 11,         &
    POISSON_NO            = -99,        &
    POISSON_NULL          = -999

  type poisson_t
    private
    type(derivatives_t), pointer, public :: der
    integer, public           :: method = POISSON_NULL
    integer, public           :: kernel
    type(cube_t), public      :: cube
    type(mesh_cube_parallel_map_t), public :: mesh_cube_map
    type(mg_solver_t) :: mg
    type(poisson_fft_t), public :: fft_solver
    FLOAT, public   :: poisson_soft_coulomb_param
    logical :: all_nodes_default
    type(poisson_corr_t) :: corrector
    type(poisson_isf_t)  :: isf_solver
    type(poisson_psolver_t) :: psolver_solver
    type(poisson_no_t) :: no_solver
    integer :: nslaves
    logical, public :: is_dressed = .false.
    type(photon_mode_t), public :: photons
    type(poisson_fmm_t)  :: params_fmm
#ifdef HAVE_MPI2
    integer         :: intercomm
    type(mpi_grp_t) :: local_grp
    logical         :: root
#endif
#ifdef HAVE_POKE
    type(PokeGrid)   :: poke_grid
    type(PokeSolver) :: poke_solver
#endif
    type(multigrid_t), allocatable, public  :: mgrid
  end type poisson_t

  integer, parameter ::             &
    CMD_FINISH = 1,                 &
    CMD_POISSON_SOLVE = 2

contains

  !-----------------------------------------------------------------
  subroutine poisson_init(this, namespace, space, der, mc, qtot, label, solver, verbose, force_serial, force_cmplx)
    type(poisson_t),             intent(inout) :: this
    type(space_t),               intent(in)    :: space
    type(namespace_t),           intent(in)    :: namespace
    type(derivatives_t), target, intent(in)    :: der
    type(multicomm_t),           intent(in)    :: mc
    FLOAT,                       intent(in)    :: qtot !< total charge
    character(len=*),  optional, intent(in)    :: label
    integer,           optional, intent(in)    :: solver
    logical,           optional, intent(in)    :: verbose
    logical,           optional, intent(in)    :: force_serial
    logical,           optional, intent(in)    :: force_cmplx

    logical :: need_cube, isf_data_is_parallel
    integer :: default_solver, default_kernel, box(MAX_DIM), fft_type, fft_library
    FLOAT :: fft_alpha
    character(len=60) :: str

    ! Make sure we do not try to initialize an already initialized solver
    ASSERT(this%method == POISSON_NULL)

    PUSH_SUB(poisson_init)

    if(optional_default(verbose,.true.)) then
      str = "Hartree"
      if(present(label)) str = trim(str) // trim(label)
      call messages_print_stress(stdout, trim(str))
    end if

    this%nslaves = 0
    this%der => der

    !%Variable DressedOrbitals
    !%Type logical
    !%Default false
    !%Section Hamiltonian::Poisson
    !%Description
    !% Allows for the calculation of coupled elecron-photon problems
    !% by applying the dressed orbital approach. Details can be found in
    !% https://arxiv.org/abs/1812.05562
    !% At the moment, N electrons in d (<=3) spatial dimensions, coupled
    !% to one photon mode can be described. The photon mode is included by
    !% raising the orbital dimension to d+1 and changing the particle interaction
    !% kernel and the local potential, where the former is included automatically,
    !% but the latter needs to by added by hand as a user_defined_potential!
    !% Coordinate 1-d: electron; coordinate d+1: photon.
    !%End
    call parse_variable(namespace, 'DressedOrbitals', .false., this%is_dressed)
    call messages_print_var_value(stdout, 'DressedOrbitals', this%is_dressed)
    if (this%is_dressed) then
      call messages_experimental('Dressed Orbitals')
      ASSERT(qtot > M_ZERO)
      call photon_mode_init(this%photons, namespace, der%mesh, der%dim-1, qtot)
      if (this%photons%nmodes > 1) then
        call messages_not_implemented('DressedOrbitals for more than one photon mode.')
      end if
    end if

#ifdef HAVE_MPI
    if(.not.optional_default(force_serial,.false.)) then
      !%Variable ParallelizationPoissonAllNodes
      !%Type logical
      !%Default true
      !%Section Execution::Parallelization
      !%Description
      !% When running in parallel, this variable selects whether the
      !% Poisson solver should divide the work among all nodes or only
      !% among the parallelization-in-domains groups.
      !%End

      call parse_variable(namespace, 'ParallelizationPoissonAllNodes', .true., this%all_nodes_default)
    else
      this%all_nodes_default = .false.
    end if
#endif

    !%Variable PoissonSolver
    !%Type integer
    !%Section Hamiltonian::Poisson
    !%Description
    !% Defines which method to use to solve the Poisson equation. Some incompatibilities apply depending on
    !% dimensionality, periodicity, etc.
    !% For a comparison of the accuracy and performance of the methods in Octopus, see P Garcia-Risue&ntilde;o,
    !% J Alberdi-Rodriguez <i>et al.</i>, <i>J. Comp. Chem.</i> <b>35</b>, 427-444 (2014)
    !% or <a href=http://arxiv.org/abs/1211.2092>arXiV</a>.
    !% Defaults:
    !% <br> 1D and 2D: <tt>fft</tt>.
    !% <br> 3D: <tt>cg_corrected</tt> if curvilinear, <tt>isf</tt> if not periodic, <tt>fft</tt> if periodic.
    !% <br> Dressed orbitals: <tt>direct_sum</tt>.
    !%Option NoPoisson -99
    !% Do not use a Poisson solver at all.
    !%Option FMM -4
    !% (Experimental) Fast multipole method. Requires FMM library.
    !%Option direct_sum -1
    !% Direct evaluation of the Hartree potential (only for finite systems).
    !%Option fft 0
    !% The Poisson equation is solved using FFTs. A cutoff technique
    !% for the Poisson kernel is selected so the proper boundary
    !% conditions are imposed according to the periodicity of the
    !% system. This can be overridden by the <tt>PoissonFFTKernel</tt>
    !% variable. To choose the FFT library use <tt>FFTLibrary</tt>
    !%Option cg 5
    !% Conjugate gradients (only for finite systems).
    !%Option cg_corrected 6
    !% Conjugate gradients, corrected for boundary conditions (only for finite systems).
    !%Option multigrid 7
    !% Multigrid method (only for finite systems).
    !%Option isf 8
    !% Interpolating Scaling Functions Poisson solver (only for finite systems).
    !%Option psolver 10
    !% Solver based on Interpolating Scaling Functions as implemented in the PSolver library.
    !% Parallelization in k-points requires <tt>PoissonSolverPSolverParallelData</tt> = no.
    !% Requires the PSolver external library.
    !%Option poke 11
    !% (Experimental) Solver from the Poke library.
    !%End

    default_solver = POISSON_FFT

    if (space%dim == 3 .and. .not. space%is_periodic()) default_solver = POISSON_ISF

    if (space%dim > 3) default_solver = POISSON_NO ! Kernel for higher dimensions is not implemented.

#ifdef HAVE_CLFFT
    ! this is disabled, since the difference between solvers are big
    ! enough to cause problems with the tests.
    ! if(accel_is_enabled()) default_solver = POISSON_FFT
#endif

    if(der%mesh%use_curvilinear) then
      select case (space%dim)
      case(1)
        default_solver = POISSON_DIRECT_SUM
      case(2)
        default_solver = POISSON_DIRECT_SUM
      case(3)
        default_solver = POISSON_CG_CORRECTED
      end select
    end if

    if (this%is_dressed) default_solver = POISSON_DIRECT_SUM

    if(.not.present(solver)) then
      call parse_variable(namespace, 'PoissonSolver', default_solver, this%method)
    else
      this%method = solver
    end if
    if(.not.varinfo_valid_option('PoissonSolver', this%method)) call messages_input_error(namespace, 'PoissonSolver')
    if(optional_default(verbose,.true.)) then
      select case(this%method)
      case (POISSON_DIRECT_SUM)
        str = "direct sum"
      case (POISSON_FMM)
        str = "fast multipole method"
      case (POISSON_FFT)
        str = "fast Fourier transform"
      case (POISSON_CG)
        str = "conjugate gradients"
      case (POISSON_CG_CORRECTED)
        str = "conjugate gradients, corrected"
      case (POISSON_MULTIGRID)
        str = "multigrid"
      case (POISSON_ISF)
        str = "interpolating scaling functions"
      case (POISSON_PSOLVER)
        str = "interpolating scaling functions (from BigDFT)"
      case (POISSON_NO)
        str = "no Poisson solver - Hartree set to 0"
      case (POISSON_POKE)
        str = "Poke library"
      end select
      write(message(1),'(a,a,a)') "The chosen Poisson solver is '", trim(str), "'"
      call messages_info(1)
    end if

    if (space%dim > 3 .and. this%method /= POISSON_NO) then
      call messages_input_error(namespace, 'PoissonSolver', 'Currently no Poisson solver is available for Dimensions > 3')
    end if

    if(this%method /= POISSON_FFT) then
      this%kernel = POISSON_FFT_KERNEL_NONE
    else

      ! Documentation in cube.F90
      call parse_variable(namespace, 'FFTLibrary', FFTLIB_FFTW, fft_library)

      !%Variable PoissonFFTKernel
      !%Type integer
      !%Section Hamiltonian::Poisson
      !%Description
      !% Defines which kernel is used to impose the correct boundary
      !% conditions when using FFTs to solve the Poisson equation. The
      !% default is selected depending on the dimensionality and
      !% periodicity of the system:
      !% <br>In 1D, <tt>spherical</tt> if finite, <tt>fft_nocut</tt> if periodic.
      !% <br>In 2D, <tt>spherical</tt> if finite, <tt>cylindrical</tt> if 1D-periodic, <tt>fft_nocut</tt> if 2D-periodic.
      !% <br>In 3D, <tt>spherical</tt> if finite, <tt>cylindrical</tt> if 1D-periodic, <tt>planar</tt> if 2D-periodic,
      !% <tt>fft_nocut</tt> if 3D-periodic.
      !% See C. A. Rozzi et al., <i>Phys. Rev. B</i> <b>73</b>, 205119 (2006) for 3D implementation and
      !% A. Castro et al., <i>Phys. Rev. B</i> <b>80</b>, 033102 (2009) for 2D implementation.
      !%Option spherical 0
      !% FFTs using spherical cutoff (in 2D or 3D).
      !%Option cylindrical 1
      !% FFTs using cylindrical cutoff (in 2D or 3D).
      !%Option planar 2
      !% FFTs using planar cutoff (in 3D).
      !%Option fft_nocut 3
      !% FFTs without using a cutoff (for fully periodic systems).
      !%Option multipole_correction 4
      !% The boundary conditions are imposed by using a multipole expansion. Only appropriate for finite systems.
      !% Further specification occurs with variables <tt>PoissonSolverBoundaries</tt> and <tt>PoissonSolverMaxMultipole</tt>.
      !%End

      select case (space%dim)
      case(1)
        if (.not. space%is_periodic()) then
          default_kernel = POISSON_FFT_KERNEL_SPH
        else
          default_kernel = POISSON_FFT_KERNEL_NOCUT
        end if
      case(2)
        if (space%periodic_dim == 2) then
          default_kernel = POISSON_FFT_KERNEL_NOCUT
        else if (space%is_periodic()) then
          default_kernel = space%periodic_dim
        else
          default_kernel = POISSON_FFT_KERNEL_SPH
        end if
      case(3)
        default_kernel = space%periodic_dim
      end select

      call parse_variable(namespace, 'PoissonFFTKernel', default_kernel, this%kernel)
      if(.not.varinfo_valid_option('PoissonFFTKernel', this%kernel)) call messages_input_error(namespace, 'PoissonFFTKernel')

      if(optional_default(verbose,.true.)) &
        call messages_print_var_option(stdout, "PoissonFFTKernel", this%kernel)

    end if

    !We assume the developer knows what he is doing by providing the solver option
    if(.not. present(solver)) then
      if (space%is_periodic() .and. this%method == POISSON_DIRECT_SUM) then
        message(1) = 'A periodic system may not use the direct_sum Poisson solver.'
        call messages_fatal(1)
      end if

      if (space%is_periodic() .and. this%method == POISSON_CG_CORRECTED) then
        message(1) = 'A periodic system may not use the cg_corrected Poisson solver.'
        call messages_fatal(1)
      end if

      if (space%is_periodic() .and. this%method == POISSON_CG) then
        message(1) = 'A periodic system may not use the cg Poisson solver.'
        call messages_fatal(1)
      end if

      if (space%is_periodic() .and. this%method == POISSON_MULTIGRID) then
        message(1) = 'A periodic system may not use the multigrid Poisson solver.'
        call messages_fatal(1)
      end if

      select case (space%dim)
      case(1)

        select case (space%periodic_dim)
        case(0)
          if( (this%method /= POISSON_FFT) .and. (this%method /= POISSON_DIRECT_SUM)) then
            message(1) = 'A finite 1D system may only use fft or direct_sum Poisson solvers.'
            call messages_fatal(1)
          end if
        case(1)
          if(this%method /= POISSON_FFT) then
            message(1) = 'A periodic 1D system may only use the fft Poisson solver.'
            call messages_fatal(1)
          end if
        end select

        if(der%mesh%use_curvilinear .and. this%method /= POISSON_DIRECT_SUM) then
          message(1) = 'If curvilinear coordinates are used in 1D, then the only working'
          message(2) = 'Poisson solver is direct_sum.'
          call messages_fatal(2)
        end if

      case(2)

        if ((this%method /= POISSON_FFT) .and. (this%method /= POISSON_DIRECT_SUM)) then
          message(1) = 'A 2D system may only use fft or direct_sum solvers.'
          call messages_fatal(1)
        end if

        if(der%mesh%use_curvilinear .and. (this%method /= POISSON_DIRECT_SUM) ) then
          message(1) = 'If curvilinear coordinates are used in 2D, then the only working'
          message(2) = 'Poisson solver is direct_sum.'
          call messages_fatal(2)
        end if

      case(3)

        if (space%is_periodic() .and. this%method == POISSON_FMM) then
          call messages_not_implemented('FMM for periodic systems')
        end if

        if (space%is_periodic() .and. this%method == POISSON_ISF) then
          call messages_write('The ISF solver can only be used for finite systems.')
          call messages_fatal()
        end if

        if (space%is_periodic() .and. this%method == POISSON_FFT .and. &
          this%kernel /= space%periodic_dim .and. this%kernel >=0 .and. this%kernel <=3) then
          write(message(1), '(a,i1,a)')'The system is periodic in ', space%periodic_dim ,' dimension(s),'
          write(message(2), '(a,i1,a)')'but Poisson solver is set for ', this%kernel, ' dimensions.'
          call messages_warning(2)
        end if

        if (space%is_periodic() .and. this%method == POISSON_FFT .and. this%kernel == POISSON_FFT_KERNEL_CORRECTED) then
          write(message(1), '(a,i1,a)')'PoissonFFTKernel = multipole_correction cannot be used for periodic systems.'
          call messages_fatal(1)
        end if

        if(der%mesh%use_curvilinear .and. (this%method/=POISSON_CG_CORRECTED)) then
          message(1) = 'If curvilinear coordinates are used, then the only working'
          message(2) = 'Poisson solver is cg_corrected.'
          call messages_fatal(2)
        end if

        select type (box => der%mesh%sb%box)
        type is (box_minimum_t)
          if (this%method == POISSON_CG_CORRECTED) then
            message(1) = 'When using the "minimum" box shape and the "cg_corrected"'
            message(2) = 'Poisson solver, we have observed "sometimes" some non-'
            message(3) = 'negligible error. You may want to check that the "fft" or "cg"'
            message(4) = 'solver are providing, in your case, the same results.'
            call messages_warning(4)
          end if
        end select

        if (this%method == POISSON_FMM) then
          call messages_experimental('FMM Poisson solver')
        end if
      end select
    end if

    if (this%method == POISSON_PSOLVER) then
#if !((defined HAVE_LIBISF) || (defined HAVE_PSOLVER))
      message(1) = "The PSolver Poisson solver cannot be used since the code was not compiled with the PSolver libary."
      call messages_fatal(1)
#endif
#ifdef HAVE_LIBISF
      message(1) = "The use of versions older than 1.8 of the PSolver library (previously known as LibISF)"
      message(2) = "are deprecated and will be removed in the next major release."
      call messages_warning(2)
#endif
    end if

    if(optional_default(verbose,.true.)) &
      call messages_print_stress(stdout)

    ! Now that we know the method, we check if we need a cube and its dimentions
    need_cube = .false.
    fft_type = FFT_REAL
    if(optional_default(force_cmplx, .false.)) fft_type = FFT_COMPLEX

    if (this%method == POISSON_ISF .or. this%method == POISSON_PSOLVER) then
      fft_type = FFT_NONE
      box(:) = der%mesh%idx%ll(:)
      need_cube = .true.
    end if

    if (this%method == POISSON_PSOLVER .and. multicomm_have_slaves(mc)) then
      call messages_not_implemented('Task parallelization with LibISF Poisson solver')
    end if

    if ( multicomm_strategy_is_parallel(mc, P_STRATEGY_KPOINTS) ) then
      ! Documentation in poisson_psolver.F90
      call parse_variable(namespace, 'PoissonSolverPSolverParallelData', .true., isf_data_is_parallel)
      if ( this%method == POISSON_PSOLVER .and. isf_data_is_parallel ) then
        call messages_not_implemented("k-point parallelization with PSolver library and PoissonSolverPSolverParallelData = yes")
      end if
      if ( this%method == POISSON_FFT .and. fft_library == FFTLIB_PFFT ) then
        call messages_not_implemented("k-point parallelization with PFFT library for Poisson solver")
      end if
    end if

    if (this%method == POISSON_FFT) then

      need_cube = .true.

      !%Variable DoubleFFTParameter
      !%Type float
      !%Default 2.0
      !%Section Mesh::FFTs
      !%Description
      !% For solving the Poisson equation in Fourier space, and for applying the local potential
      !% in Fourier space, an auxiliary cubic mesh is built. This mesh will be larger than
      !% the circumscribed cube of the usual mesh by a factor <tt>DoubleFFTParameter</tt>. See
      !% the section that refers to Poisson equation, and to the local potential for details
      !% [the default value of two is typically good].
      !%End
      call parse_variable(namespace, 'DoubleFFTParameter', M_TWO, fft_alpha)
      if (fft_alpha < M_ONE .or. fft_alpha > M_THREE ) then
        write(message(1), '(a,f12.5,a)') "Input: '", fft_alpha, &
          "' is not a valid DoubleFFTParameter"
        message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
        call messages_fatal(2)
      end if

      if (space%dim /= 3 .and. fft_library == FFTLIB_PFFT) then
        call messages_not_implemented('PFFT support for dimensionality other than 3')
      end if

      select case (space%dim)

      case (1)
        select case(this%kernel)
        case(POISSON_FFT_KERNEL_SPH)
          call mesh_double_box(space, der%mesh, fft_alpha, box)
        case(POISSON_FFT_KERNEL_NOCUT)
          box = der%mesh%idx%ll
        end select

      case (2)
        select case(this%kernel)
        case(POISSON_FFT_KERNEL_SPH)
          call mesh_double_box(space, der%mesh, fft_alpha, box)
          box(1:2) = maxval(box)
        case(POISSON_FFT_KERNEL_CYL)
          call mesh_double_box(space, der%mesh, fft_alpha, box)
        case(POISSON_FFT_KERNEL_NOCUT)
          box(:) = der%mesh%idx%ll(:)
        end select

      case (3)
        select case(this%kernel)
        case(POISSON_FFT_KERNEL_SPH)
          call mesh_double_box(space, der%mesh, fft_alpha, box)
          box(:) = maxval(box)
        case(POISSON_FFT_KERNEL_CYL)
          call mesh_double_box(space, der%mesh, fft_alpha, box)
          box(2) = maxval(box(2:3)) ! max of finite directions
          box(3) = maxval(box(2:3)) ! max of finite directions
        case(POISSON_FFT_KERNEL_CORRECTED)
          box(:) = der%mesh%idx%ll(:)
        case(POISSON_FFT_KERNEL_PLA, POISSON_FFT_KERNEL_NOCUT)
          call mesh_double_box(space, der%mesh, fft_alpha, box)
        end select

      end select

    end if

    if(this%method == POISSON_POKE) then
#ifndef HAVE_POKE
      call messages_write('Octopus was compiled without Poke support, you cannot use', new_line = .true.)
      call messages_write("  'PoissonSolver = poke'. ")
      call messages_fatal()
#endif

      call messages_experimental('Poke library')
      ASSERT(space%dim == 3)
      box(1:space%dim) = der%mesh%idx%ll(1:space%dim)
      need_cube = .true.
      fft_type = FFTLIB_NONE
    end if

    ! Create the cube
    if (need_cube) then
      call cube_init(this%cube, box, der%mesh%sb, namespace, fft_type = fft_type, &
                     need_partition=.not.der%mesh%parallel_in_domains)
      if (this%cube%parallel_in_domains .and. this%method == POISSON_FFT) then
        call mesh_cube_parallel_map_init(this%mesh_cube_map, der%mesh, this%cube)
      end if
    end if

    if(this%method == POISSON_POKE) then

#ifdef HAVE_POKE
      this%poke_grid = PokeGrid(der%mesh%spacing, this%cube%rs_n)
      if (space%is_periodic()) then
        call this%poke_grid%set_boundaries(POKE_BOUNDARIES_PERIODIC)
      else
        call this%poke_grid%set_boundaries(POKE_BOUNDARIES_FREE)
      end if
      this%poke_solver = PokeSolver(this%poke_grid)
      call this%poke_solver%build()
#endif
    end if

    if (this%is_dressed .and. .not. this%method == POISSON_DIRECT_SUM) then
      write(message(1), '(a)')'Dressed Orbital calculation currently only working with direct sum Poisson solver.'
      call messages_fatal(1)
    end if

    call poisson_kernel_init(this, namespace, space, mc%master_comm)

    POP_SUB(poisson_init)
  end subroutine poisson_init

  !-----------------------------------------------------------------
  subroutine poisson_end(this)
    type(poisson_t), intent(inout) :: this

    logical :: has_cube

    PUSH_SUB(poisson_end)

    has_cube = .false.

    select case(this%method)
    case(POISSON_FFT)
      call poisson_fft_end(this%fft_solver)
      if(this%kernel == POISSON_FFT_KERNEL_CORRECTED) call poisson_corrections_end(this%corrector)
      has_cube = .true.

    case(POISSON_CG_CORRECTED, POISSON_CG)
      call poisson_cg_end()
      call poisson_corrections_end(this%corrector)

    case(POISSON_MULTIGRID)
      call poisson_multigrid_end(this%mg)

    case(POISSON_ISF)
      call poisson_isf_end(this%isf_solver)
      has_cube = .true.

    case(POISSON_PSOLVER)
      call poisson_psolver_end(this%psolver_solver)
      has_cube = .true.

    case(POISSON_FMM)
      call poisson_fmm_end(this%params_fmm)

    case(POISSON_NO)
      call poisson_no_end(this%no_solver)

    case(POISSON_POKE)
#ifdef HAVE_POKE
      call this%poke_grid%end()
      call this%poke_solver%end()
#endif

    end select
    this%method = POISSON_NULL

    if (has_cube) then
      if (this%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_end(this%mesh_cube_map)
      end if
      call cube_end(this%cube)
    end if

    if (this%is_dressed) then
      call photon_mode_end(this%photons)
    end if
    this%is_dressed = .false.

    if (allocated(this%mgrid)) then
      call multigrid_end(this%mgrid)
      SAFE_DEALLOCATE_A(this%mgrid)
    end if

    POP_SUB(poisson_end)
  end subroutine poisson_end

  !-----------------------------------------------------------------

  subroutine zpoisson_solve_real_and_imag_separately(this, pot, rho, all_nodes, kernel)
    type(poisson_t),                    intent(in)    :: this
    CMPLX,                              intent(inout) :: pot(:)  !< pot(mesh%np)
    CMPLX,                              intent(in)    :: rho(:)  !< rho(mesh%np)
    logical, optional,                  intent(in)    :: all_nodes
    type(fourier_space_op_t), optional, intent(in)    :: kernel

    FLOAT, allocatable :: aux1(:), aux2(:)
    type(derivatives_t), pointer :: der
    logical :: all_nodes_value

    type(profile_t), save :: prof

    der => this%der

    PUSH_SUB(zpoisson_solve_real_and_imag_separately)

    call profiling_in(prof, 'POISSON_RE_IM_SOLVE')

    if(present(kernel)) then
      ASSERT(.not. any(abs(kernel%qq(:))>CNST(1e-8)))
    end if

    all_nodes_value = optional_default(all_nodes, this%all_nodes_default)

    SAFE_ALLOCATE(aux1(1:der%mesh%np))
    SAFE_ALLOCATE(aux2(1:der%mesh%np))
    ! first the real part
    aux1(1:der%mesh%np) = real(rho(1:der%mesh%np))
    aux2(1:der%mesh%np) = real(pot(1:der%mesh%np))
    call dpoisson_solve(this, aux2, aux1, all_nodes=all_nodes_value, kernel=kernel)
    pot(1:der%mesh%np)  = aux2(1:der%mesh%np)

    ! now the imaginary part
    aux1(1:der%mesh%np) = aimag(rho(1:der%mesh%np))
    aux2(1:der%mesh%np) = aimag(pot(1:der%mesh%np))
    call dpoisson_solve(this, aux2, aux1, all_nodes=all_nodes_value, kernel=kernel)
    pot(1:der%mesh%np) = pot(1:der%mesh%np) + M_zI*aux2(1:der%mesh%np)

    SAFE_DEALLOCATE_A(aux1)
    SAFE_DEALLOCATE_A(aux2)

    call profiling_out(prof)

    POP_SUB(zpoisson_solve_real_and_imag_separately)
  end subroutine zpoisson_solve_real_and_imag_separately

  !-----------------------------------------------------------------

  subroutine zpoisson_solve(this, pot, rho, all_nodes, kernel)
    type(poisson_t),                    intent(in)    :: this
    CMPLX,                              intent(inout) :: pot(:)  !< pot(mesh%np)
    CMPLX,                              intent(in)    :: rho(:)  !< rho(mesh%np)
    logical, optional,                  intent(in)    :: all_nodes
    type(fourier_space_op_t), optional, intent(in)    :: kernel

    logical :: all_nodes_value
    type(profile_t), save :: prof

    PUSH_SUB(zpoisson_solve)

    all_nodes_value = optional_default(all_nodes, this%all_nodes_default)

    ASSERT(ubound(pot, dim = 1) == this%der%mesh%np_part .or. ubound(pot, dim = 1) == this%der%mesh%np)
    ASSERT(ubound(rho, dim = 1) == this%der%mesh%np_part .or. ubound(rho, dim = 1) == this%der%mesh%np)

    ASSERT(this%method /= POISSON_NULL)

    if(this%method == POISSON_FFT .and. this%kernel /= POISSON_FFT_KERNEL_CORRECTED  &
          .and. .not. this%is_dressed) then
      !The default (real) Poisson solver is used for OEP and Sternheimer calls were we do not need
      !a complex-to-xomplex FFT as these parts use the normal Coulomb potential
      if(this%cube%fft%type == FFT_COMPLEX) then
        !We add the profiling here, as the other path uses dpoisson_solve
        call profiling_in(prof, 'ZPOISSON_SOLVE')
        call zpoisson_fft_solve(this%fft_solver, this%der%mesh, this%cube, pot, rho, this%mesh_cube_map, kernel=kernel)
        call profiling_out(prof)
      else
        call zpoisson_solve_real_and_imag_separately(this, pot, rho, all_nodes_value, kernel=kernel)
      end if
    else
      call zpoisson_solve_real_and_imag_separately(this, pot, rho, all_nodes_value, kernel = kernel)
    end if

    POP_SUB(zpoisson_solve)
  end subroutine zpoisson_solve


  !-----------------------------------------------------------------

  subroutine poisson_solve_batch(this, potb, rhob, all_nodes, kernel)
    type(poisson_t),                    intent(inout) :: this
    type(batch_t),                      intent(inout) :: potb
    type(batch_t),                      intent(inout) :: rhob
    logical, optional,                  intent(in)    :: all_nodes
    type(fourier_space_op_t), optional, intent(in)    :: kernel

    integer :: ii

    PUSH_SUB(poisson_solve_batch)

    ASSERT(potb%nst_linear == rhob%nst_linear)
    ASSERT(potb%type() == rhob%type())

    if(potb%type() == TYPE_FLOAT) then
      do ii = 1, potb%nst_linear
        call dpoisson_solve(this, potb%dff_linear(:, ii), rhob%dff_linear(:, ii), all_nodes, kernel=kernel)
      end do
    else
      do ii = 1, potb%nst_linear
        call zpoisson_solve(this, potb%zff_linear(:, ii), rhob%zff_linear(:, ii), all_nodes, kernel=kernel)
      end do
    end if

    POP_SUB(poisson_solve_batch)
  end subroutine poisson_solve_batch

  !-----------------------------------------------------------------

  !> Calculates the Poisson equation.
  !! Given the density returns the corresponding potential.
  !!
  !! Different solvers are available that can be chosen in the input file
  !! with the "PoissonSolver" parameter
  subroutine dpoisson_solve(this, pot, rho, all_nodes, kernel)
    type(poisson_t),                    intent(in)    :: this
    FLOAT,                              intent(inout) :: pot(:) !< Local size of the \b potential vector.
    FLOAT,                              intent(inout) :: rho(:) !< Local size of the \b density (rho) vector.
    !> Is the Poisson solver allowed to utilise
    !! all nodes or only the domain nodes for
    !! its calculations? (Defaults to .true.)
    logical, optional,                  intent(in)    :: all_nodes
    type(fourier_space_op_t), optional, intent(in)    :: kernel

    type(derivatives_t), pointer :: der
    type(cube_function_t) :: crho, cpot
    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)
    logical               :: all_nodes_value
    type(profile_t), save :: prof

    call profiling_in(prof, 'POISSON_SOLVE')
    PUSH_SUB(dpoisson_solve)

    der => this%der

    ASSERT(ubound(pot, dim = 1) == der%mesh%np_part .or. ubound(pot, dim = 1) == der%mesh%np)
    ASSERT(ubound(rho, dim = 1) == der%mesh%np_part .or. ubound(rho, dim = 1) == der%mesh%np)

    ! Check optional argument and set to default if necessary.
    all_nodes_value = optional_default(all_nodes, this%all_nodes_default)

    ASSERT(this%method /= POISSON_NULL)

    if(present(kernel)) then
      ASSERT(this%method == POISSON_FFT)
    end if

    select case(this%method)
    case(POISSON_DIRECT_SUM)
      if ( (this%is_dressed .and. this%der%dim - 1 > 3) .or. this%der%dim > 3) then
        message(1) = "Direct sum Poisson solver only available for 1, 2, or 3 dimensions."
        call messages_fatal(1)
      end if
      call poisson_solve_direct(this, pot, rho)

    case(POISSON_FMM)
      call poisson_fmm_solve(this%params_fmm, this%der, pot, rho)

    case(POISSON_CG)
      call poisson_cg1(der, this%corrector, pot, rho)

    case(POISSON_CG_CORRECTED)
      SAFE_ALLOCATE(rho_corrected(1:der%mesh%np))
      SAFE_ALLOCATE(vh_correction(1:der%mesh%np_part))

      call correct_rho(this%corrector, der, rho, rho_corrected, vh_correction)

      pot(1:der%mesh%np) = pot(1:der%mesh%np) - vh_correction(1:der%mesh%np)
      call poisson_cg2(der, pot, rho_corrected)
      pot(1:der%mesh%np) = pot(1:der%mesh%np) + vh_correction(1:der%mesh%np)

      SAFE_DEALLOCATE_A(rho_corrected)
      SAFE_DEALLOCATE_A(vh_correction)

    case(POISSON_MULTIGRID)
      call poisson_multigrid_solver(this%mg, der, pot, rho)

    case(POISSON_FFT)
      if(this%kernel /= POISSON_FFT_KERNEL_CORRECTED) then
        call dpoisson_fft_solve(this%fft_solver, der%mesh, this%cube, pot, rho, this%mesh_cube_map, kernel=kernel)
      else
        SAFE_ALLOCATE(rho_corrected(1:der%mesh%np))
        SAFE_ALLOCATE(vh_correction(1:der%mesh%np_part))

        call correct_rho(this%corrector, der, rho, rho_corrected, vh_correction)
        call dpoisson_fft_solve(this%fft_solver, der%mesh, this%cube, pot, rho_corrected, this%mesh_cube_map, &
          average_to_zero = .true., kernel=kernel)

        pot(1:der%mesh%np) = pot(1:der%mesh%np) + vh_correction(1:der%mesh%np)
        SAFE_DEALLOCATE_A(rho_corrected)
        SAFE_DEALLOCATE_A(vh_correction)
      end if

    case(POISSON_ISF)
      call poisson_isf_solve(this%isf_solver, der%mesh, this%cube, pot, rho, all_nodes_value)


    case(POISSON_PSOLVER)
      if (this%psolver_solver%datacode == "G") then
        ! Global version
        call poisson_psolver_global_solve(this%psolver_solver, der%mesh, this%cube, pot, rho)
      else ! "D" Distributed version
        call poisson_psolver_parallel_solve(this%psolver_solver, der%mesh, this%cube, pot, rho, this%mesh_cube_map)
      end if

    case(POISSON_POKE)
      call dcube_function_alloc_RS(this%cube, crho)
      call dcube_function_alloc_RS(this%cube, cpot)
      call dmesh_to_cube(der%mesh, rho, this%cube, crho)
#if HAVE_POKE
      call this%poke_solver%solve(crho%drs, cpot%drs)
#endif
      call dcube_to_mesh(this%cube, cpot, der%mesh, pot)
      call dcube_function_free_RS(this%cube, crho)
      call dcube_function_free_RS(this%cube, cpot)

    case(POISSON_NO)
      call poisson_no_solve(this%no_solver, der%mesh, this%cube, pot, rho)
    end select


    ! Add extra terms for dressed interaction
    if (this%is_dressed .and. this%method /= POISSON_NO) then
      call photon_mode_add_poisson_terms(this%photons, der%mesh, rho, pot)
    end if

    POP_SUB(dpoisson_solve)
    call profiling_out(prof)
  end subroutine dpoisson_solve

  !-----------------------------------------------------------------
  subroutine poisson_init_sm(this, namespace, space, main, der, sm, method, force_cmplx)
    type(poisson_t),             intent(inout) :: this
    type(namespace_t),           intent(in)    :: namespace
    type(space_t),               intent(in)    :: space
    type(poisson_t),             intent(in)    :: main
    type(derivatives_t), target, intent(in)    :: der
    type(submesh_t),             intent(inout) :: sm
    integer, optional,           intent(in)    :: method
    logical, optional,           intent(in)    :: force_cmplx

    integer :: default_solver, idir
    integer :: box(MAX_DIM)
    FLOAT   :: qq(1:MAX_DIM)

    if(this%method /= POISSON_NULL) return ! already initialized

    PUSH_SUB(poisson_init_sm)

    this%is_dressed = .false.
    !TODO: To be implemented as an option
    this%all_nodes_default = .false.

    this%nslaves = 0
    this%der => der

#ifdef HAVE_MPI
    this%all_nodes_default = main%all_nodes_default
#endif

    default_solver = POISSON_DIRECT_SUM
    this%method = default_solver
    if(present(method)) this%method = method

    if(der%mesh%use_curvilinear) then
      call messages_not_implemented("Submesh Poisson solver with curvilinear mesh")
    end if

    this%kernel = POISSON_FFT_KERNEL_NONE

    select case(this%method)
    case(POISSON_DIRECT_SUM)
      !Nothing to be done

    case(POISSON_ISF)
      !TODO: Add support for domain parrallelization
      ASSERT(.not. der%mesh%parallel_in_domains)
      call submesh_get_cube_dim(sm, box, der%dim)
      call submesh_init_cube_map(sm, der%dim)
      call cube_init(this%cube, box, der%mesh%sb, namespace, fft_type = FFT_NONE, &
                     need_partition=.not.der%mesh%parallel_in_domains)
      call poisson_isf_init(this%isf_solver, namespace, der%mesh, this%cube, mpi_world%comm, init_world = this%all_nodes_default)

    case(POISSON_PSOLVER)
      !TODO: Add support for domain parrallelization
      ASSERT(.not. der%mesh%parallel_in_domains)
      if(this%all_nodes_default) then
        this%cube%mpi_grp = mpi_world
      else
        this%cube%mpi_grp = this%der%mesh%mpi_grp
      end if
      call submesh_get_cube_dim(sm, box, der%dim)
      call submesh_init_cube_map(sm, der%dim)
      call cube_init(this%cube, box, der%mesh%sb, namespace, fft_type = FFT_NONE, &
                     need_partition=.not.der%mesh%parallel_in_domains)
      qq = M_ZERO
      call poisson_psolver_init(this%psolver_solver, namespace, space, this%der%mesh, this%cube, M_ZERO, qq, force_isolated=.true.)
      call poisson_psolver_get_dims(this%psolver_solver, this%cube)
    case(POISSON_FFT)
      !Here we impose zero boundary conditions
      this%kernel = POISSON_FFT_KERNEL_SPH
      !We need to parse this, in case this routine is called before poisson_init
      call parse_variable(namespace, 'FFTLibrary', FFTLIB_FFTW, fft_default_lib)

      call submesh_get_cube_dim(sm, box, der%dim)
      call submesh_init_cube_map(sm, der%dim)
      !We double the size of the cell
      !Maybe the factor of two should be controlled as a variable
      do idir = 1, der%dim
        box(idir) = nint(M_TWO * (box(idir) - 1)) + 1
      end do
      if(optional_default(force_cmplx, .false.)) then
        call cube_init(this%cube, box, der%mesh%sb, namespace, fft_type = FFT_COMPLEX, &
                       need_partition=.not.der%mesh%parallel_in_domains)
      else
        call cube_init(this%cube, box, der%mesh%sb, namespace, fft_type = FFT_REAL, &
                       need_partition=.not.der%mesh%parallel_in_domains)
      end if
      call poisson_fft_init(this%fft_solver, namespace, space, this%der%mesh, this%cube, this%kernel)
    end select

    POP_SUB(poisson_init_sm)
  end subroutine poisson_init_sm

  !-----------------------------------------------------------------
  !> This routine checks the Hartree solver selected in the input
  !! file by calculating numerically and analytically the Hartree
  !! potential originated by a Gaussian distribution of charge.
  !! For periodic systems, the periodic copies of the Gaussian
  !! are taken into account up to to a certain threshold that can
  !! be specified in the input file.
  subroutine poisson_test(this, space, mesh, namespace, repetitions)
    type(poisson_t),   intent(in) :: this
    type(space_t),     intent(in) :: space
    type(mesh_t),      intent(in) :: mesh
    type(namespace_t), intent(in) :: namespace
    integer,           intent(in) :: repetitions

    FLOAT, allocatable :: rho(:), vh(:), vh_exact(:), xx(:, :), xx_per(:)
    FLOAT :: alpha, beta, rr, delta, ralpha, hartree_nrg_num, &
         hartree_nrg_analyt, lcl_hartree_nrg
    FLOAT :: total_charge
    integer :: ip, ierr, iunit, nn, n_gaussians, itime, icell
    FLOAT :: threshold
    type(lattice_iterator_t) :: latt_iter

    PUSH_SUB(poisson_test)

    if(mesh%sb%dim == 1) then
      call messages_not_implemented('Poisson test for 1D case')
    end if

    !%Variable PoissonTestPeriodicThreshold
    !%Type float
    !%Default 1e-5
    !%Section Hamiltonian::Poisson
    !%Description
    !% This threshold determines the accuracy of the periodic copies of
    !% the Gaussian charge distribution that are taken into account when
    !% computing the analytical solution for periodic systems.
    !% Be aware that the default leads to good results for systems
    !% that are periodic in 1D - for 3D it is very costly because of the
    !% large number of copies needed.
    !%End
    call parse_variable(namespace, 'PoissonTestPeriodicThreshold', CNST(1e-5), threshold)

    ! Use two gaussians with different sign
    n_gaussians = 2

    SAFE_ALLOCATE(     rho(1:mesh%np))
    SAFE_ALLOCATE(      vh(1:mesh%np))
    SAFE_ALLOCATE(vh_exact(1:mesh%np))
    SAFE_ALLOCATE(xx(1:space%dim, 1:n_gaussians))
    SAFE_ALLOCATE(xx_per(1:space%dim))

    rho = M_ZERO; vh = M_ZERO; vh_exact = M_ZERO

    alpha = CNST(4.0)*mesh%spacing(1)
    write(message(1),'(a,f14.6)')  "Info: The alpha value is ", alpha
    write(message(2),'(a)')        "      Higher values of alpha lead to more physical densities and more reliable results."
    call messages_info(2)
    beta = M_ONE / ( alpha**space%dim * sqrt(M_PI)**space%dim )

    write(message(1), '(a)') 'Building the Gaussian distribution of charge...'
    call messages_info(1)

    ! Set the centers of the Gaussians by hand
    xx(1, 1) = M_ONE
    xx(2, 1) = -M_HALF
    if(space%dim == 3) xx(3, 1) = M_TWO
    xx(1, 2) = -M_TWO
    xx(2, 2) = M_ZERO
    if(space%dim == 3) xx(3, 2) = -M_ONE
    xx = xx * alpha

    ! Density as sum of Gaussians
    rho = M_ZERO
    do nn = 1, n_gaussians
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin = xx(:, nn))
        rho(ip) = rho(ip) + (-1)**nn * beta*exp(-(rr/alpha)**2)
      end do
    end do

    total_charge = dmf_integrate(mesh, rho)

    write(message(1), '(a,f14.6)') 'Total charge of the Gaussian distribution', total_charge
    call messages_info(1)

    write(message(1), '(a)') 'Computing exact potential.'
    call messages_info(1)

    ! This builds analytically its potential
    vh_exact = M_ZERO
    latt_iter = lattice_iterator_t(mesh%sb%latt, M_ONE/threshold)
    do nn = 1, n_gaussians
      ! sum over all periodic copies for each Gaussian
      write(message(1), '(a,i2,a,i9,a)') 'Computing Gaussian ', nn, ' for ', latt_iter%n_cells, ' periodic copies.'
      call messages_info(1)

      do icell = 1, latt_iter%n_cells
        xx_per = xx(:, nn) + latt_iter%get(icell)
        !$omp parallel do private(rr, ralpha)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin=xx_per)
          select case(space%dim)
          case(3)
            if(rr > R_SMALL) then
              vh_exact(ip) = vh_exact(ip) + (-1)**nn * loct_erf(rr/alpha)/rr
            else
              vh_exact(ip) = vh_exact(ip) + (-1)**nn * (M_TWO/sqrt(M_PI))/alpha
            end if
          case(2)
            ralpha = rr**2/(M_TWO*alpha**2)
            if(ralpha < CNST(100.0)) then
              vh_exact(ip) = vh_exact(ip) + (-1)**nn * beta * (M_PI)**(M_THREE*M_HALF) * alpha * exp(-rr**2/(M_TWO*alpha**2)) * &
                loct_bessel_in(0, rr**2/(M_TWO*alpha**2))
            else
              vh_exact(ip) = vh_exact(ip) + (-1)**nn * beta * (M_PI)**(M_THREE*M_HALF) * alpha * &
                            (M_ONE/sqrt(M_TWO*M_PI*ralpha))
            end if
          end select
        end do
      end do
    end do

    ! This calculates the numerical potential
    do itime = 1, repetitions
      call dpoisson_solve(this, vh, rho)
    end do

    ! Output results
    iunit = io_open("hartree_results", namespace, action='write')
    delta = dmf_nrm2(mesh, vh-vh_exact)
    write(iunit, '(a,f19.13)' ) 'Hartree test (abs.) = ', delta
    delta = delta/dmf_nrm2(mesh, vh_exact)
    write(iunit, '(a,f19.13)' ) 'Hartree test (rel.) = ', delta

    ! Calculate the numerical Hartree energy (serially)
    lcl_hartree_nrg = M_ZERO
    do ip = 1, mesh%np
      lcl_hartree_nrg = lcl_hartree_nrg + rho(ip) * vh(ip)
    end do
    lcl_hartree_nrg = lcl_hartree_nrg * mesh%spacing(1) * mesh%spacing(2) * mesh%spacing(3)/M_TWO
#ifdef HAVE_MPI
    call MPI_Reduce(lcl_hartree_nrg, hartree_nrg_num, 1, &
         MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    if(mpi_err /= 0) then
      write(message(1),'(a)') "MPI error in MPI_Reduce; subroutine poisson_test of file poisson.F90"
      call messages_warning(1)
    end if
#else
    hartree_nrg_num = lcl_hartree_nrg
#endif

    ! Calculate the analytical Hartree energy (serially, discrete - not exactly exact)
    lcl_hartree_nrg = M_ZERO
    do ip = 1, mesh%np
      lcl_hartree_nrg = lcl_hartree_nrg + rho(ip) * vh_exact(ip)
    end do
    lcl_hartree_nrg = lcl_hartree_nrg * mesh%spacing(1) * mesh%spacing(2) * mesh%spacing(3)/M_TWO
#ifdef HAVE_MPI
    call MPI_Reduce(lcl_hartree_nrg, hartree_nrg_analyt, 1, &
         MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    if(mpi_err /= 0) then
      write(message(1),'(a)') "MPI error in MPI_Reduce; subroutine poisson_test of file poisson.F90"
      call messages_warning(1)
    end if
#else
    hartree_nrg_analyt = lcl_hartree_nrg
#endif

    write(iunit, '(a,f19.13)' )

    if (mpi_world%rank == 0) then
      write(iunit,'(a,f19.13)') 'Hartree Energy (numerical) =',hartree_nrg_num,'Hartree Energy (analytical) =',hartree_nrg_analyt
    end if

    call io_close(iunit)

    call dio_function_output (io_function_fill_how('AxisX'), ".", "poisson_test_rho", namespace, &
      mesh, rho, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisX'), ".", "poisson_test_exact", namespace, &
      mesh, vh_exact, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisX'), ".", "poisson_test_numerical", namespace, &
      mesh, vh, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisY'), ".", "poisson_test_rho", namespace, &
      mesh, rho, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisY'), ".", "poisson_test_exact", namespace, &
      mesh, vh_exact, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisY'), ".", "poisson_test_numerical", namespace, &
      mesh, vh, unit_one, ierr)
    ! not dimensionless, but no need for unit conversion for a test routine

    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(vh)
    SAFE_DEALLOCATE_A(vh_exact)
    SAFE_DEALLOCATE_A(xx)

    POP_SUB(poisson_test)
  end subroutine poisson_test

  ! -----------------------------------------------------------------

  logical pure function poisson_solver_is_iterative(this) result(iterative)
    type(poisson_t), intent(in) :: this

    iterative = this%method == POISSON_CG .or. this%method == POISSON_CG_CORRECTED .or. this%method == POISSON_MULTIGRID
  end function poisson_solver_is_iterative

  ! -----------------------------------------------------------------

  logical pure function poisson_is_multigrid(this) result(is_multigrid)
    type(poisson_t), intent(in) :: this

    is_multigrid = (this%method == POISSON_MULTIGRID)

  end function poisson_is_multigrid

  ! -----------------------------------------------------------------

  logical pure function poisson_solver_has_free_bc(this) result(free_bc)
    type(poisson_t), intent(in) :: this

    free_bc = .true.

    if (this%method == POISSON_FFT .and. &
      this%kernel /= POISSON_FFT_KERNEL_SPH .and. this%kernel /= POISSON_FFT_KERNEL_CORRECTED) then
      free_bc = .false.
    end if

  end function poisson_solver_has_free_bc

  !-----------------------------------------------------------------

  integer pure function poisson_get_solver(this) result (solver)
    type(poisson_t), intent(in) :: this

    solver = this%method
  end function poisson_get_solver

  !-----------------------------------------------------------------

  subroutine poisson_async_init(this, mc)
    type(poisson_t), intent(inout) :: this
    type(multicomm_t), intent(in)  :: mc

    PUSH_SUB(poisson_async_init)

#ifdef HAVE_MPI2
    if(multicomm_have_slaves(mc)) then

      call mpi_grp_init(this%local_grp, mc%group_comm(P_STRATEGY_STATES))

      this%root = (this%local_grp%rank == 0)

      this%intercomm = mc%slave_intercomm
      call MPI_Comm_remote_size(this%intercomm, this%nslaves, mpi_err)

    end if
#endif

    POP_SUB(poisson_async_init)

  end subroutine poisson_async_init

  !-----------------------------------------------------------------

  subroutine poisson_async_end(this, mc)
    type(poisson_t), intent(inout) :: this
    type(multicomm_t), intent(in)  :: mc

#ifdef HAVE_MPI2
    integer :: islave
#endif

    PUSH_SUB(poisson_async_end)

#ifdef HAVE_MPI2
    if(multicomm_have_slaves(mc)) then

      ! send the finish signal
      do islave = this%local_grp%rank, this%nslaves - 1, this%local_grp%size
        call MPI_Send(M_ONE, 1, MPI_FLOAT, islave, CMD_FINISH, this%intercomm, mpi_err)
      end do

    end if
#endif

    POP_SUB(poisson_async_end)

  end subroutine poisson_async_end

  !-----------------------------------------------------------------

  subroutine poisson_slave_work(this)
    type(poisson_t), intent(inout) :: this

#ifdef HAVE_MPI2
    FLOAT, allocatable :: rho(:), pot(:)
    logical :: done
    integer :: status(MPI_STATUS_SIZE)
    type(profile_t), save :: prof, bcast_prof, wait_prof
    integer :: bcast_root

    PUSH_SUB(poisson_slave_work)
    call profiling_in(prof, "SLAVE_WORK")

    SAFE_ALLOCATE(rho(1:this%der%mesh%np))
    SAFE_ALLOCATE(pot(1:this%der%mesh%np))
    done = .false.

    do while(.not. done)

      call profiling_in(wait_prof, "SLAVE_WAIT")
      call MPI_Recv(rho(1), this%der%mesh%np, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, this%intercomm, status(1), mpi_err)
      call profiling_out(wait_prof)

      ! The tag of the message tells us what we have to do.
      select case(status(MPI_TAG))

      case(CMD_FINISH)
        done = .true.

      case(CMD_POISSON_SOLVE)
        call dpoisson_solve(this, pot, rho)

        call profiling_in(bcast_prof, "SLAVE_BROADCAST")
        bcast_root = MPI_PROC_NULL
        if(this%root) bcast_root = MPI_ROOT
        call MPI_Bcast(pot(1), this%der%mesh%np, MPI_FLOAT, bcast_root, this%intercomm, mpi_err)
        call profiling_out(bcast_prof)

      end select

    end do

    SAFE_DEALLOCATE_A(pot)
    SAFE_DEALLOCATE_A(rho)

    call profiling_out(prof)
    POP_SUB(poisson_slave_work)
#endif
  end subroutine poisson_slave_work

  !----------------------------------------------------------------

  logical pure function poisson_is_async(this) result(async)
    type(poisson_t),  intent(in) :: this

    async = (this%nslaves > 0)

  end function poisson_is_async

  !----------------------------------------------------------------

  subroutine poisson_build_kernel(this, namespace, space, coulb, qq, mu, singul)
    type(poisson_t),          intent(in)    :: this
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(fourier_space_op_t), intent(inout) :: coulb
    FLOAT,                    intent(in)    :: qq(:)
    FLOAT,                    intent(in)    :: mu
    FLOAT,          optional, intent(in)    :: singul

    PUSH_SUB(poisson_build_kernel)

    if (space%is_periodic()) then
      ASSERT(ubound(qq, 1) >= space%periodic_dim)
      ASSERT(this%method == POISSON_FFT)
    end if

    if(mu > M_EPSILON) then
      if(this%method /= POISSON_FFT) then
        write(message(1),'(a)') "Poisson solver with range separation is only implemented with FFT."
        call messages_fatal(1)
      end if
      coulb%mu = mu
    end if

    !TODO: this should be a select case supporting other kernels.
    ! This means that we need an abstract object for kernels.
    select case(this%method)
    case(POISSON_FFT)
      !We only reinitialize the poisson sover if needed
      if(any(abs(coulb%qq(1:space%periodic_dim) - qq(1:space%periodic_dim)) > M_EPSILON)) then
        call fourier_space_op_end(coulb)
        coulb%qq(1:space%periodic_dim) = qq(1:space%periodic_dim)
        !We must define the singularity if we specify a q vector and we do not use the short-range Coulomb potential
        coulb%singularity = optional_default(singul, M_ZERO)
        call poisson_fft_get_kernel(namespace, space, this%der%mesh, this%cube, coulb, this%kernel, &
          this%poisson_soft_coulomb_param)
      end if
    case default
      call messages_not_implemented("poisson_build_kernel with other methods than FFT")
    end select


    POP_SUB(poisson_build_kernel)
  end subroutine poisson_build_kernel

#include "poisson_init_inc.F90"
#include "poisson_direct_inc.F90"
#include "poisson_direct_sm_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "poisson_inc.F90"
#include "undef.F90"
#include "complex.F90"
#include "poisson_inc.F90"

end module poisson_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
