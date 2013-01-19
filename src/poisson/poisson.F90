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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module poisson_m
  use boundaries_m
  use cube_m
  use datasets_m
  use derivatives_m
  use fft_m
  use geometry_m
  use global_m
  use index_m
  use io_m
  use io_function_m
  use loct_math_m
  use math_m
  use mesh_m
  use mesh_cube_parallel_map_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multicomm_m
  use opencl_m
  use par_vec_m
  use parser_m
  use poisson_cg_m
  use poisson_corrections_m
  use poisson_isf_m
  use poisson_fft_m
  use poisson_fmm_m
  use poisson_multigrid_m
  use poisson_sete_m
  use profiling_m
  use simul_box_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use xyz_file_m      

  implicit none

  private
  public ::                      &
    poisson_t,                   &
    poisson_fmm_t,               &
    poisson_get_solver,          &
    poisson_init,                &
    dpoisson_solve,              &
    zpoisson_solve,              &
    poisson_solver_is_iterative, &
    poisson_solver_has_free_bc,  &
    poisson_end,                 &
    poisson_test,                &
    poisson_is_multigrid,        &
    poisson_energy,              &
    poisson_slave_work,          &
    poisson_async_init,          &
    poisson_async_end,           &
    dpoisson_solve_start,        &
    dpoisson_solve_finish,       &
    zpoisson_solve_start,        &
    zpoisson_solve_finish,       &
    poisson_is_async

  integer, public, parameter ::         &
    POISSON_DIRECT_SUM_1D = -1,         &
    POISSON_DIRECT_SUM_2D = -2,         &
    POISSON_DIRECT_SUM_3D = -3,         &  
    POISSON_FMM           = -4,         &
    POISSON_FFT           =  0,         &
    POISSON_CG            =  5,         &
    POISSON_CG_CORRECTED  =  6,         &
    POISSON_MULTIGRID     =  7,         &
    POISSON_ISF           =  8,         &
    POISSON_SETE          =  9
  
  type poisson_t
    type(derivatives_t), pointer :: der
    integer           :: method = -99
    integer           :: kernel
    type(cube_t)      :: cube
    type(mesh_cube_parallel_map_t) :: mesh_cube_map
    type(mg_solver_t) :: mg
    type(poisson_fft_t) :: fft_solver
    FLOAT   :: poisson_soft_coulomb_param = M_ONE
    logical :: all_nodes_default
    type(poisson_corr_t) :: corrector
    type(poisson_sete_t) :: sete_solver
    type(poisson_isf_t)  :: isf_solver
    type(poisson_fmm_t)  :: params_fmm
    integer :: nslaves
#ifdef HAVE_MPI2
    integer         :: intercomm
    type(mpi_grp_t) :: local_grp
    logical         :: root
#endif

  end type poisson_t

  type(poisson_t), target, save, public :: psolver

  integer, parameter ::             &
    CMD_FINISH = 1,                 &
    CMD_POISSON_SOLVE = 2

contains

  !-----------------------------------------------------------------
  subroutine poisson_init(this, der, geo, all_nodes_comm)
    type(poisson_t),             intent(out) :: this
    type(derivatives_t), target, intent(in)  :: der
    type(geometry_t),            intent(in)  :: geo
    integer,                     intent(in)  :: all_nodes_comm

    logical :: need_cube
    integer :: default_solver, default_kernel, box(MAX_DIM), fft_type

    if(this%method /= -99) return ! already initialized

    PUSH_SUB(poisson_init)

    call messages_print_stress(stdout, "Hartree")

    this%nslaves = 0
    this%der => der

#ifdef HAVE_MPI
    !%Variable ParallelizationPoissonAllNodes
    !%Type logical
    !%Default true
    !%Section Execution::Parallelization
    !%Description
    !% When running in parallel, this variable selects whether the
    !% Poisson solver should divide the work among all nodes or only
    !% among the parallelization-in-domains groups.
    !%End

    call parse_logical(datasets_check('ParallelizationPoissonAllNodes'), .true., this%all_nodes_default)
#endif

    !%Variable PoissonSolver
    !%Type integer
    !%Section Hamiltonian::Poisson
    !%Description
    !% Defines which method to use to solve the Poisson equation. Some incompatibilities apply depending on
    !% dimensionality, periodicity, etc.
    !% For comparison of accuracy and performance of the methods in Octopus, see http://arxiv.org/abs/1211.2092. 
    !% Defaults:
    !% <br> 1D and 2D: <tt>fft</tt>.
    !% <br> 3D: <tt>cg_corrected</tt> if curvilinear, <tt>isf</tt> if not periodic, <tt>fft</tt> if periodic.
    !%Option NoPoisson -999
    !% Do not use a Poisson solver at all.
    !%Option FMM -4
    !% Fast multipole method.                                  
    !%Option direct3D -3                                      
    !% Direct evaluation of the Hartree potential (in 3D).   
    !%Option direct2D -2
    !% Direct evaluation of the Hartree potential (in 2D).
    !%Option direct1D -1
    !% Direct evaluation of the Hartree potential (in 1D).
    !%Option fft 0
    !% The Poisson equation is solved using FFTs. A cutoff technique
    !% for the Poisson kernel is selected so the proper boundary
    !% conditions are imposed according to the periodicity of the
    !% system. This can be overridden by the <tt>PoissonFFTKernel</tt>
    !% variable.
    !%Option cg 5
    !% Conjugate gradients.
    !%Option cg_corrected 6
    !% Conjugate gradients, corrected for boundary conditions.
    !% Further specification occurs with variables <tt>PoissonSolverBoundaries</tt> and <tt>PoissonSolverMaxMultipole</tt>.
    !%Option multigrid 7
    !% Multigrid method.
    !%Option isf 8
    !% Interpolating Scaling Functions Poisson solver.
    !%Option sete 9
    !% (Experimental) SETE solver. Complex boundaries are required.
    !%End

    default_solver = POISSON_FFT

#ifndef SINGLE_PRECISION
    if(der%mesh%sb%dim == 3 .and. der%mesh%sb%periodic_dim == 0) default_solver = POISSON_ISF
#endif

    if(der%mesh%sb%dim > 3) default_solver = POISSON_CG_CORRECTED

#ifdef HAVE_CLAMDFFT
    ! this is disabled, since the difference between solvers are big
    ! enough to cause problems with the tests.
    ! if(opencl_is_enabled()) default_solver = POISSON_FFT
#endif

    if(der%mesh%use_curvilinear) then
      select case(der%mesh%sb%dim)
      case(1)
        default_solver = POISSON_DIRECT_SUM_1D
      case(2)
        default_solver = POISSON_DIRECT_SUM_2D
      case(3)
        default_solver = POISSON_CG_CORRECTED
      end select
    end if

    call parse_integer(datasets_check('PoissonSolver'), default_solver, this%method)
    if(.not.varinfo_valid_option('PoissonSolver', this%method)) call input_error('PoissonSolver')

    call messages_print_var_option(stdout, "PoissonSolver", this%method)

    if(this%method /= POISSON_FFT) then
      this%kernel = POISSON_FFT_KERNEL_NONE
    else

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

      select case(der%mesh%sb%dim)
      case(1)
        if(der%mesh%sb%periodic_dim == 0) then
          default_kernel = POISSON_FFT_KERNEL_SPH
        else
          default_kernel = POISSON_FFT_KERNEL_NOCUT
        end if
      case(2)
        if (der%mesh%sb%periodic_dim == 2) then
          default_kernel = POISSON_FFT_KERNEL_NOCUT
        else if (der%mesh%sb%periodic_dim > 0) then
          default_kernel = der%mesh%sb%periodic_dim
        else
          default_kernel = POISSON_FFT_KERNEL_SPH
        end if
      case(3)
        default_kernel = der%mesh%sb%periodic_dim
      end select

      call parse_integer(datasets_check('PoissonFFTKernel'), default_kernel, this%kernel)
      if(.not.varinfo_valid_option('PoissonFFTKernel', this%method)) call input_error('PoissonFFTKernel')

      call messages_print_var_option(stdout, "PoissonFFTKernel", this%kernel)

    end if

    select case(der%mesh%sb%dim)
    case(1)

      select case(der%mesh%sb%periodic_dim)
      case(0)
        if( (this%method /= POISSON_FFT) .and. (this%method /= POISSON_DIRECT_SUM_1D)) then
          message(1) = 'A finite 1D system may only use fft or direct1D Poisson solvers.'
          call messages_fatal(1)
        end if
      case(1)
        message(1) = 'A periodic 1D system may only use the fft Poisson solver.'
        call messages_fatal(1)
      end select

      if(der%mesh%use_curvilinear .and. this%method /= POISSON_DIRECT_SUM_1D) then
        message(1) = 'If curvilinear coordinates are used in 1D, then the only working'
        message(2) = 'Poisson solver is direct1D.'
        call messages_fatal(2)
      end if

    case(2)

      if( (this%method /= POISSON_FFT) .and. (this%method /= POISSON_DIRECT_SUM_2D) ) then
        message(1) = 'A 2D system may only use fft or direct2D Poisson solvers.'
        call messages_fatal(1)
      end if

      if(der%mesh%use_curvilinear .and. (this%method /= POISSON_DIRECT_SUM_2D) ) then
        message(1) = 'If curvilinear coordinates are used in 2D, then the only working'
        message(2) = 'Poisson solver is direct2D.'
        call messages_fatal(2)
      end if

    case(3)
      
      if(this%method == POISSON_DIRECT_SUM_1D .or. this%method == POISSON_DIRECT_SUM_2D) then
        message(1) = 'A 3D system may not use the direct1D or direct2D Poisson solvers.'
        call messages_fatal(1)
      end if

      if(der%mesh%sb%periodic_dim > 0 .and. this%method == POISSON_FMM) then
        write(message(1), '(a,i1,a)')'FMM is not ready to deal with periodic boundaries at present.'
        call messages_warning(1)
      end if

      if(der%mesh%sb%periodic_dim > 0 .and. this%method == POISSON_ISF) then
        call messages_write('The ISF solver should only be used for finite systems.')
        call messages_warning()
      end if

      if(der%mesh%sb%periodic_dim > 0 .and. this%method == POISSON_FFT .and. &
        this%kernel /= der%mesh%sb%periodic_dim .and. this%kernel >=0 .and. this%kernel <=3) then
        write(message(1), '(a,i1,a)')'The system is periodic in ', der%mesh%sb%periodic_dim ,' dimension(s),'
        write(message(2), '(a,i1,a)')'but Poisson solver is set for ', this%kernel, ' dimensions.'
        call messages_warning(2)
      end if

      if(der%mesh%sb%periodic_dim > 0 .and. this%method == POISSON_FFT .and. &
        this%kernel == POISSON_FFT_KERNEL_CORRECTED) then
        write(message(1), '(a,i1,a)')'PoissonFFTKernel = multipole_correction cannot be used for periodic systems.'
        call messages_fatal(1)
      end if

      if(der%mesh%use_curvilinear .and. (this%method/=POISSON_CG_CORRECTED)) then
        message(1) = 'If curvilinear coordinates are used, then the only working'
        message(2) = 'Poisson solver is cg_corrected.'
        call messages_fatal(2)
      end if

      if( (der%mesh%sb%box_shape == MINIMUM) .and. (this%method == POISSON_CG_CORRECTED) ) then
        message(1) = 'When using the "minimum" box shape and the "cg_corrected"'
        message(2) = 'Poisson solver, we have observed "sometimes" some non-'
        message(3) = 'negligible error. You may want to check that the "fft" or "cg"'
        message(4) = 'solver are providing, in your case, the same results.'
        call messages_warning(4)
      end if

      if (this%method == POISSON_SETE) then
        call messages_experimental('SETE poisson solver')

        if(.not. simul_box_complex_boundaries(der%mesh%sb)) then
          message(1) = 'Complex boundaries must be enabled to use the SETE Poisson solver.'
          message(2) = 'Use ComplexBoundaries = yes in your input file.'
          call messages_fatal(2)
        end if
      end if

      if (this%method == POISSON_FMM) then
        call messages_experimental('FMM Poisson solver')
      end if

    end select

    ! Now that we know the method, we check if we need a cube and its dimentions
    need_cube = .false.
    fft_type = FFT_REAL

    if (this%method == POISSON_ISF) then
      fft_type = FFT_NONE
      box(:) = der%mesh%idx%ll(:)
      need_cube = .true.
    end if

    if (this%method == POISSON_FFT) then

      need_cube = .true.

      select case (der%mesh%sb%dim)

      case (1)
        select case(this%kernel)
        case(POISSON_FFT_KERNEL_SPH)
          call mesh_double_box(der%mesh%sb, der%mesh, box)
        case(POISSON_FFT_KERNEL_NOCUT)
          box = der%mesh%idx%ll
        end select

      case (2)
        select case(this%kernel)
        case(POISSON_FFT_KERNEL_SPH)
          call mesh_double_box(der%mesh%sb, der%mesh, box)
          box(1:2) = maxval(box)
        case(POISSON_FFT_KERNEL_CYL)
          call mesh_double_box(der%mesh%sb, der%mesh, box)
        case(POISSON_FFT_KERNEL_NOCUT)
          box(:) = der%mesh%idx%ll(:)
        end select

      case (3)
        select case(this%kernel)
        case(POISSON_FFT_KERNEL_SPH) 
          call mesh_double_box(der%mesh%sb, der%mesh, box)
          box(:) = maxval(box)
        case(POISSON_FFT_KERNEL_CORRECTED)
          box(:) = der%mesh%idx%ll(:)
        case(POISSON_FFT_KERNEL_CYL, POISSON_FFT_KERNEL_PLA, POISSON_FFT_KERNEL_NOCUT)
          call mesh_double_box(der%mesh%sb, der%mesh, box)
        end select

      end select

    end if

    ! Create the cube
    if (need_cube) then
      call cube_init(this%cube, box, der%mesh%sb, fft_type = fft_type, verbose = .true.)
      if (der%mesh%parallel_in_domains .and. this%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_init(this%mesh_cube_map, der%mesh, this%cube)
      end if
    end if

    select case(der%mesh%sb%dim)
    case (1)
      call poisson1d_init(this)
    case (2)
      call poisson2D_init(this)
    case (3)
      call poisson3D_init(this, geo, all_nodes_comm)
    end select

    call messages_print_stress(stdout)

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

    case(POISSON_SETE)
      call poisson_sete_end(this%sete_solver)

    ! We should clean up and release
    ! case(POISSON_DIRECT_SUM_3D)  
    ! call poisson_direct_sum_end
    case(POISSON_FMM)
      call poisson_fmm_end(this%params_fmm)

    end select
    this%method = -99

    if (has_cube) then
      if (this%der%mesh%parallel_in_domains .and. this%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_end(this%mesh_cube_map)
      end if
      call cube_end(this%cube)
    end if

    POP_SUB(poisson_end)
  end subroutine poisson_end

  !-----------------------------------------------------------------

  subroutine zpoisson_solve_1(this, pot, rho, all_nodes)
    type(poisson_t),      intent(inout) :: this
    CMPLX,                intent(inout) :: pot(:)  !< pot(mesh%np)
    CMPLX,                intent(in)    :: rho(:)  !< rho(mesh%np)
    logical, optional,    intent(in)    :: all_nodes

    FLOAT, allocatable :: aux1(:), aux2(:)
    type(derivatives_t), pointer :: der
    logical :: all_nodes_value

    der => this%der

    PUSH_SUB(zpoisson_solve_1)

    if(present(all_nodes)) then
      all_nodes_value = all_nodes
    else
      all_nodes_value = this%all_nodes_default
    end if

    SAFE_ALLOCATE(aux1(1:der%mesh%np))
    SAFE_ALLOCATE(aux2(1:der%mesh%np))

    ! first the real part
    aux1(1:der%mesh%np) = real(rho(1:der%mesh%np))
    aux2(1:der%mesh%np) = real(pot(1:der%mesh%np))
    call dpoisson_solve(this, aux2, aux1, all_nodes=all_nodes_value)
    pot(1:der%mesh%np)  = aux2(1:der%mesh%np)

    ! now the imaginary part
    aux1(1:der%mesh%np) = aimag(rho(1:der%mesh%np))
    aux2(1:der%mesh%np) = aimag(pot(1:der%mesh%np))
    call dpoisson_solve(this, aux2, aux1, all_nodes=all_nodes_value)
    pot(1:der%mesh%np) = pot(1:der%mesh%np) + M_zI*aux2(1:der%mesh%np)

    SAFE_DEALLOCATE_A(aux1)
    SAFE_DEALLOCATE_A(aux2)

    POP_SUB(zpoisson_solve_1)
  end subroutine zpoisson_solve_1

  !-----------------------------------------------------------------

  subroutine zpoisson_solve(this, pot, rho, all_nodes, theta)
    type(poisson_t),      intent(inout) :: this
    CMPLX,                intent(inout) :: pot(:)  !< pot(mesh%np)
    CMPLX,                intent(in)    :: rho(:)  !< rho(mesh%np)
    logical, optional,    intent(in)    :: all_nodes
    FLOAT, optional,      intent(in)    :: theta   !< complex scaling angle

    logical :: all_nodes_value, cmplxscl
    FLOAT :: theta_

    PUSH_SUB(zpoisson_solve)

    if(present(all_nodes)) then
      all_nodes_value = all_nodes
    else
      all_nodes_value = this%all_nodes_default
    end if

    cmplxscl = .false.
    theta_ = M_ZERO
    if(present(theta))  then
      cmplxscl = .true.
      theta_ = theta
    end if

    ASSERT(this%method.ne.-99)
      
    select case(this%method)

    case(POISSON_DIRECT_SUM_1D)
      call zpoisson1D_solve(this, pot, rho, theta_)

    case(POISSON_DIRECT_SUM_2D)
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    case(POISSON_DIRECT_SUM_3D)
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    case(POISSON_FMM)
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    case(POISSON_CG)
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    case(POISSON_CG_CORRECTED)
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    case(POISSON_MULTIGRID)
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    case(POISSON_FFT)
      if(this%der%mesh%sb%dim == 1 .and. cmplxscl) then
        call messages_not_implemented('Complex scaled 1D soft coulomb and FFT poisson solver')
      end if
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    case(POISSON_ISF)
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    case(POISSON_SETE)
      call zpoisson_solve_1(this, pot, rho, all_nodes)
      if(cmplxscl) pot = pot * exp(- M_zI * theta_)

    end select



    POP_SUB(zpoisson_solve)
  end subroutine zpoisson_solve

  !-----------------------------------------------------------------

  !> Calculates the Poisson equation.
  !! Given the density returns the corresponding potential.
  !!
  !! Different solvers are available that could be choosen in the input file
  !! with the "PoissonSolver" parameter
  subroutine dpoisson_solve(this, pot, rho, all_nodes)
    type(poisson_t),      intent(inout) :: this
    FLOAT,                intent(inout) :: pot(:)
    FLOAT,                intent(inout) :: rho(:)
    !> Is the Poisson solver allowed to utilise
    !! all nodes or only the domain nodes for
    !! its calculations? (Defaults to .true.)
    logical, optional,    intent(in)    :: all_nodes 
    type(derivatives_t), pointer :: der
    integer :: counter
    integer :: nx_half, nx
    integer :: ny_half, ny
    integer :: nz_half, nz

    FLOAT :: xl, yl, zl

    FLOAT, allocatable :: vh0(:,:,:), rh0(:,:,:)

    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)

    logical               :: all_nodes_value
    type(profile_t), save :: prof
    integer               :: conversion(3)

    call profiling_in(prof, 'POISSON_SOLVE')
    PUSH_SUB(dpoisson_solve)

    der => this%der

    ASSERT(ubound(pot, dim = 1) == der%mesh%np_part .or. ubound(pot, dim = 1) == der%mesh%np)
    ASSERT(ubound(rho, dim = 1) == der%mesh%np_part .or. ubound(rho, dim = 1) == der%mesh%np)

    ! Check optional argument and set to default if necessary.
    if(present(all_nodes)) then
      all_nodes_value = all_nodes
    else
      all_nodes_value = this%all_nodes_default
    end if

    ASSERT(this%method /= -99)
      
    select case(this%method)
    case(POISSON_DIRECT_SUM_1D)
      call poisson1d_solve(this, pot, rho)

    case(POISSON_DIRECT_SUM_2D)
      call poisson2d_solve(this, pot, rho)

    case(POISSON_DIRECT_SUM_3D)
      call poisson3D_solve_direct(this, pot, rho)

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
        call poisson_fft_solve(this%fft_solver, der%mesh, this%cube, pot, rho, this%mesh_cube_map)
      else
        SAFE_ALLOCATE(rho_corrected(1:der%mesh%np))
        SAFE_ALLOCATE(vh_correction(1:der%mesh%np_part))
        
        call correct_rho(this%corrector, der, rho, rho_corrected, vh_correction)
        call poisson_fft_solve(this%fft_solver, der%mesh, this%cube, pot, rho_corrected, this%mesh_cube_map, &
          average_to_zero = .true.)
        
        pot(1:der%mesh%np) = pot(1:der%mesh%np) + vh_correction(1:der%mesh%np)
        SAFE_DEALLOCATE_A(rho_corrected)
        SAFE_DEALLOCATE_A(vh_correction)
      end if

    case(POISSON_ISF)
      call poisson_isf_solve(this%isf_solver, der%mesh, this%cube, pot, rho, all_nodes_value)
     
    case(POISSON_SETE)

      nx = der%mesh%idx%nr(2,1) - der%mesh%idx%nr(1,1) + 1 - 2*der%mesh%idx%enlarge(1)
      ny = der%mesh%idx%nr(2,2) - der%mesh%idx%nr(1,2) + 1 - 2*der%mesh%idx%enlarge(2)
      nz = der%mesh%idx%nr(2,3) - der%mesh%idx%nr(1,3) + 1 - 2*der%mesh%idx%enlarge(3)

      nx_half = (der%mesh%idx%nr(2,1) - der%mesh%idx%nr(1,1) - 2*der%mesh%idx%enlarge(1))/2 + 1
      ny_half = (der%mesh%idx%nr(2,2) - der%mesh%idx%nr(1,2) - 2*der%mesh%idx%enlarge(2))/2 + 1
      nz_half = (der%mesh%idx%nr(2,3) - der%mesh%idx%nr(1,3) - 2*der%mesh%idx%enlarge(3))/2 + 1

      SAFE_ALLOCATE(rh0(1:nx, 1:ny, 1:nz))
      SAFE_ALLOCATE(vh0(1:nx, 1:ny, 1:nz))

      ! The %size gives out half the size of a box. 
      xl = 2*der%mesh%sb%lsize(1) 
      yl = 2*der%mesh%sb%lsize(2) 
      zl = 2*der%mesh%sb%lsize(3)

      do counter = 1, der%mesh%np
        call  index_to_coords(der%mesh%idx,der%mesh%sb%dim,counter,conversion)
        conversion(1)=conversion(1)+nx_half
        conversion(2)=conversion(2)+ny_half
        conversion(3)=conversion(3)+nz_half
        vh0(conversion(1),conversion(2),conversion(3))=pot(counter)
        rh0(conversion(1),conversion(2),conversion(3))=rho(counter)
      end do

      call poisson_sete_solve(this%sete_solver, rh0, vh0, nx, ny, nz)

      do counter = 1, der%mesh%np

        call  index_to_coords(der%mesh%idx,der%mesh%sb%dim,counter,conversion)

        conversion(1) = conversion(1) + nx_half
        conversion(2) = conversion(2) + ny_half
        conversion(3) = conversion(3) + nz_half
        
        pot(counter)=vh0(conversion(1), conversion(2), conversion(3))

      end do

      SAFE_DEALLOCATE_A(vh0)
      SAFE_DEALLOCATE_A(rh0)

    end select

    POP_SUB(dpoisson_solve)
    call profiling_out(prof)
  end subroutine dpoisson_solve


  !-----------------------------------------------------------------
  !> This routine checks the Hartree solver selected in the input
  !! file by calculating numerically and analytically the Hartree
  !! potential originated by a Gaussian distribution of charge.
  !! This only makes sense for finite systems.
  subroutine poisson_test(mesh)
    type(mesh_t), intent(inout) :: mesh
    FLOAT, allocatable :: rho(:), vh(:), vh_exact(:), rhop(:), xx(:, :)
    FLOAT :: alpha, beta, rr, delta, norm, ralpha, hartree_nrg_num, &
         hartree_nrg_analyt, lcl_hartree_nrg 
    FLOAT :: total_charge
    integer :: ip, idir, ierr, iunit, nn, n_gaussians

    PUSH_SUB(poisson_test)

    if(mesh%sb%dim == 1) then
      call messages_not_implemented('Poisson test for 1D case')
    endif

    n_gaussians = 1 

    SAFE_ALLOCATE(     rho(1:mesh%np))
    SAFE_ALLOCATE(    rhop(1:mesh%np))
    SAFE_ALLOCATE(      vh(1:mesh%np))
    SAFE_ALLOCATE(vh_exact(1:mesh%np))
    SAFE_ALLOCATE(xx(1:mesh%sb%dim, 1:n_gaussians))

    rho = M_ZERO; vh = M_ZERO; vh_exact = M_ZERO; rhop = M_ZERO

    alpha = CNST(4.0)*mesh%spacing(1)
    write(message(1),'(a,f14.6)')  "Info: The alpha value is ", alpha
    write(message(2),'(a)')        "      Higher values of alpha lead to more physical densities and more reliable results."
    call messages_info(2)
    beta = M_ONE / ( alpha**mesh%sb%dim * sqrt(M_PI)**mesh%sb%dim )

    write(message(1), '(a)') 'Building the Gaussian distribution of charge...'
    call messages_info(1)

    rho = M_ZERO
    do nn = 1, n_gaussians
      norm = M_ZERO

      do idir = 1, mesh%sb%dim
        xx(idir, nn) = 0.000 
      end do

      rr = sqrt(sum(xx(:, nn)*xx(:,nn)))
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin = xx(:, nn))
        rhop(ip) = beta*exp(-(rr/alpha)**2)
      end do

      norm = dmf_integrate(mesh, rhop)
      
      rhop = (-1)**nn * rhop
      do ip = 1, mesh%np 
        rho(ip) = rho(ip) + rhop(ip)
      end do
    end do

    total_charge = dmf_integrate(mesh, rho)

    write(message(1), '(a,f14.6)') 'Total charge of the Gaussian distribution', total_charge
    call messages_info(1)

    ! This builds analytically its potential
    vh_exact = M_ZERO
    do nn = 1, n_gaussians
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin = xx(:, nn))
        select case(mesh%sb%dim)
        case(3)
          if(rr > r_small) then
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

    ! This calculates the numerical potential
    call dpoisson_solve(psolver, vh, rho)

    ! Output results
    iunit = io_open("hartree_results", action='write')
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

    ! Calculate the anallytical Hartree energy (serially, discrete - not exactly exact)
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
    
    call dio_function_output (io_function_fill_how('AxisX'), ".", "poisson_test_rho", mesh, rho, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisX'), ".", "poisson_test_exact", mesh, vh_exact, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisX'), ".", "poisson_test_numerical", mesh, vh, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisY'), ".", "poisson_test_rho", mesh, rho, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisY'), ".", "poisson_test_exact", mesh, vh_exact, unit_one, ierr)
    call dio_function_output (io_function_fill_how('AxisY'), ".", "poisson_test_numerical", mesh, vh, unit_one, ierr)
    ! not dimensionless, but no need for unit conversion for a test routine

    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(rhop)
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

    if (this%method == POISSON_SETE) free_bc = .false.

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
  
  FLOAT function poisson_energy(this) result(energy)
    type(poisson_t), intent(in) :: this

    PUSH_SUB(poisson_energy)

    energy = M_ZERO
    if(this%method == POISSON_SETE) energy = poisson_sete_energy(this%sete_solver)

    POP_SUB(poisson_energy)
  end function poisson_energy

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

    call profiling_out(prof)
    POP_SUB(poisson_slave_work)
#endif
  end subroutine poisson_slave_work

  !------------------------------------------------------------------

  subroutine dpoisson_solve_start(this, rho)
    type(poisson_t),      intent(inout) :: this
    FLOAT,                intent(in)    :: rho(:)

#ifdef HAVE_MPI2    
    integer :: islave
    type(profile_t), save :: prof

    PUSH_SUB(dpoisson_solve_start)
    call profiling_in(prof, "POISSON_START")

    ! we assume all nodes have a copy of the density
    do islave = this%local_grp%rank, this%nslaves - 1, this%local_grp%size !all nodes are used for communication
      call MPI_Send(rho(1), this%der%mesh%np, MPI_FLOAT, islave, CMD_POISSON_SOLVE, this%intercomm, mpi_err)
    end do
    
    call profiling_out(prof)
    POP_SUB(dpoisson_solve_start)
#endif
    
  end subroutine dpoisson_solve_start
  
  !----------------------------------------------------------------

  subroutine dpoisson_solve_finish(this, pot)
    type(poisson_t),  intent(inout) :: this
    FLOAT,            intent(inout) :: pot(:)

#ifdef HAVE_MPI2
    type(profile_t), save :: prof

    PUSH_SUB(dpoisson_solve_finish)
    call profiling_in(prof, "POISSON_FINISH")

    call MPI_Bcast(pot(1), this%der%mesh%np, MPI_FLOAT, 0, this%intercomm, mpi_err)

    call profiling_out(prof)
    POP_SUB(dpoisson_solve_finish)
#endif
  end subroutine dpoisson_solve_finish

  !----------------------------------------------------------------

  !------------------------------------------------------------------

  subroutine zpoisson_solve_start(this, rho)
    type(poisson_t),      intent(inout) :: this
    CMPLX,                intent(in)    :: rho(:)

#ifdef HAVE_MPI2    
    integer :: islave
    type(profile_t), save :: prof

    PUSH_SUB(zpoisson_solve_start)
    call profiling_in(prof, "ZPOISSON_START")

    ! we assume all nodes have a copy of the density
    do islave = this%local_grp%rank, this%nslaves - 1, this%local_grp%size !all nodes are used for communication
      call MPI_Send(rho(1), this%der%mesh%np, MPI_CMPLX, islave, CMD_POISSON_SOLVE, this%intercomm, mpi_err)
    end do
    
    call profiling_out(prof)
    POP_SUB(zpoisson_solve_start)
#endif
    
  end subroutine zpoisson_solve_start
  
  !----------------------------------------------------------------

  subroutine zpoisson_solve_finish(this, pot)
    type(poisson_t),  intent(inout) :: this
    CMPLX,            intent(inout) :: pot(:)

#ifdef HAVE_MPI2
    type(profile_t), save :: prof

    PUSH_SUB(zpoisson_solve_finish)
    call profiling_in(prof, "ZPOISSON_FINISH")

    call MPI_Bcast(pot(1), this%der%mesh%np, MPI_CMPLX, 0, this%intercomm, mpi_err)

    call profiling_out(prof)
    POP_SUB(zpoisson_solve_finish)
#endif
  end subroutine zpoisson_solve_finish

  !----------------------------------------------------------------


  logical pure function poisson_is_async(this) result(async)
    type(poisson_t),  intent(in) :: this
    
    async = (this%nslaves > 0)

  end function poisson_is_async

#include "solver_1d_inc.F90"
#include "solver_2d_inc.F90"
#include "solver_3d_inc.F90"

end module poisson_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
