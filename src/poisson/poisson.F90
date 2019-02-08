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
  use boundaries_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use derivatives_oct_m
  use fft_oct_m
  use global_oct_m
  use index_oct_m
  use io_oct_m
  use io_function_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
#ifdef HAVE_OPENMP
  use omp_lib
#endif
  use par_vec_oct_m
  use parser_oct_m
  use partition_oct_m
  use poisson_cg_oct_m
  use poisson_corrections_oct_m
  use poisson_isf_oct_m
  use poisson_fft_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                      &
    poisson_t,                   &
    poisson_get_solver,          &
    poisson_get_qpoint,          &
    poisson_init,                &
    poisson_kernel_reinit,       &
    dpoisson_solve,              &
    zpoisson_solve,              &
    poisson_solve_batch,         &
    poisson_solver_is_iterative, &
    poisson_solver_has_free_bc,  &
    poisson_end,                 &
    poisson_test,                &
    poisson_slave_work,          &
    poisson_async_init,          &
    poisson_async_end,           &
    dpoisson_solve_start,        &
    dpoisson_solve_finish,       &
    zpoisson_solve_start,        &
    zpoisson_solve_finish,       &
    poisson_is_async

  integer, public, parameter ::         &
    POISSON_FFT           =  0,         &
    POISSON_CG            =  5,         &
    POISSON_CG_CORRECTED  =  6,         &
    POISSON_ISF           =  8,         &
    POISSON_NULL          = -999
  
  type poisson_t
    type(derivatives_t), pointer :: der
    integer           :: method = POISSON_NULL
    integer           :: kernel
    type(cube_t)      :: cube
    type(mesh_cube_parallel_map_t) :: mesh_cube_map
    type(poisson_fft_t) :: fft_solver
    logical :: all_nodes_default
    type(poisson_corr_t) :: corrector
    type(poisson_isf_t)  :: isf_solver
    integer :: nslaves
    FLOAT :: qq(MAX_DIM) !< for exchange in periodic system
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
  subroutine poisson_init(this, der, mc, label, qq, solver)
    type(poisson_t),             intent(out) :: this
    type(derivatives_t), target, intent(in)  :: der
    type(multicomm_t),           intent(in)  :: mc
    character(len=*),  optional, intent(in)  :: label
    FLOAT,             optional, intent(in)  :: qq(:) !< (der%mesh%sb%periodic_dim)
    integer,           optional, intent(in)  :: solver

    logical :: need_cube
    integer :: default_solver, default_kernel, box(MAX_DIM), fft_type, fft_library
    FLOAT :: fft_alpha
    character(len=60) :: str

    if(this%method /= POISSON_NULL) return ! already initialized

    PUSH_SUB(poisson_init)

    str = "Hartree"
    if(present(label)) str = trim(str) // trim(label)
    call messages_print_stress(stdout, trim(str))

    this%nslaves = 0
    this%der => der

    this%qq = M_ZERO
    if(present(qq)  .and. simul_box_is_periodic(der%mesh%sb)) then
      ASSERT(ubound(qq, 1) >= der%mesh%sb%periodic_dim)
      ASSERT(this%method == POISSON_FFT)
      this%qq(1:der%mesh%sb%periodic_dim) = qq(1:der%mesh%sb%periodic_dim)
    end if

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

    call parse_variable('ParallelizationPoissonAllNodes', .true., this%all_nodes_default)
#endif

    !%Variable PoissonSolver
    !%Type integer
    !%Section Hamiltonian::Poisson
    !%Description
    !% Defines which method to use to solve the Poisson equation. Some incompatibilities apply depending on
    !% dimensionality, periodicity, etc.
    !% Defaults: <br> 3D: <tt>isf</tt> if not periodic, <tt>fft</tt> if periodic.
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
    !%Option isf 8
    !% Interpolating Scaling Functions Poisson solver (only for finite systems).
    !%End
    default_solver = POISSON_FFT
    if(der%mesh%sb%periodic_dim == 0) default_solver = POISSON_ISF

    
#ifdef HAVE_CLFFT
    ! this is disabled, since the difference between solvers are big
    ! enough to cause problems with the tests.
    ! if(accel_is_enabled()) default_solver = POISSON_FFT
#endif

    if(.not.present(solver)) then
      call parse_variable('PoissonSolver', default_solver, this%method)
    else
      this%method = solver
    end if
    if(.not.varinfo_valid_option('PoissonSolver', this%method)) call messages_input_error('PoissonSolver')
   
    select case(this%method)
    case (POISSON_FFT)
      str = "fast Fourier transform"
    case (POISSON_CG)
      str = "conjugate gradients"
    case (POISSON_CG_CORRECTED)
      str = "conjugate gradients, corrected"
    case (POISSON_ISF)
      str = "interpolating scaling functions"
    end select
    write(message(1),'(a,a,a)') "The chosen Poisson solver is '", trim(str), "'"
    call messages_info(1)

    if(this%method /= POISSON_FFT) then
      this%kernel = POISSON_FFT_KERNEL_NONE
    else

      ! Documentation in cube.F90
      call parse_variable('FFTLibrary', FFTLIB_FFTW, fft_library)
      
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

      call parse_variable('PoissonFFTKernel', default_kernel, this%kernel)
      if(.not.varinfo_valid_option('PoissonFFTKernel', this%kernel)) call messages_input_error('PoissonFFTKernel')

      call messages_print_var_option(stdout, "PoissonFFTKernel", this%kernel)

    end if

    if(.not. present(solver)) then 
      if(der%mesh%sb%periodic_dim > 0 .and. this%method == POISSON_CG_CORRECTED) then
        message(1) = 'A periodic system may not use the cg_corrected Poisson solver.'
        call messages_fatal(1)
      end if

      if(der%mesh%sb%periodic_dim > 0 .and. this%method == POISSON_CG) then
        message(1) = 'A periodic system may not use the cg Poisson solver.'
        call messages_fatal(1)
      end if

      if(der%mesh%sb%periodic_dim > 0 .and. this%method == POISSON_ISF) then
        call messages_write('The ISF solver can only be used for finite systems.')
        call messages_fatal()
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
      
      if( (der%mesh%sb%box_shape == MINIMUM) .and. (this%method == POISSON_CG_CORRECTED) ) then
        message(1) = 'When using the "minimum" box shape and the "cg_corrected"'
        message(2) = 'Poisson solver, we have observed "sometimes" some non-'
        message(3) = 'negligible error. You may want to check that the "fft" or "cg"'
        message(4) = 'solver are providing, in your case, the same results.'
        call messages_warning(4)
      end if
    end if

    call messages_print_stress(stdout)

    ! Now that we know the method, we check if we need a cube and its dimentions
    need_cube = .false.
    fft_type = FFT_REAL

    if (this%method == POISSON_ISF) then
      fft_type = FFT_NONE
      box(:) = der%mesh%idx%ll(:)
      need_cube = .true.
    end if

    if ( multicomm_strategy_is_parallel(mc, P_STRATEGY_KPOINTS) ) then
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
      call parse_variable('DoubleFFTParameter', M_TWO, fft_alpha)
      if (fft_alpha < M_ONE .or. fft_alpha > M_THREE ) then
        write(message(1), '(a,f12.5,a)') "Input: '", fft_alpha, &
          "' is not a valid DoubleFFTParameter"
        message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
        call messages_fatal(2)
      end if

      if (der%mesh%sb%periodic_dim /= 0 .and. fft_library == FFTLIB_PFFT) then
        call messages_not_implemented('PFFT support for periodic systems')
      end if

      select case(this%kernel)
      case(POISSON_FFT_KERNEL_SPH) 
        call mesh_double_box(der%mesh%sb, der%mesh, fft_alpha, box)
        box(:) = maxval(box)
      case(POISSON_FFT_KERNEL_CYL) 
        call mesh_double_box(der%mesh%sb, der%mesh, fft_alpha, box)
        box(2) = maxval(box(2:3)) ! max of finite directions
        box(3) = maxval(box(2:3)) ! max of finite directions
      case(POISSON_FFT_KERNEL_CORRECTED)
        box(:) = der%mesh%idx%ll(:)
      case(POISSON_FFT_KERNEL_PLA, POISSON_FFT_KERNEL_NOCUT)
        call mesh_double_box(der%mesh%sb, der%mesh, fft_alpha, box)
      end select
    end if

    ! Create the cube
    if (need_cube) then
      call cube_init(this%cube, box, der%mesh%sb, fft_type = fft_type, verbose = .true., &
                     need_partition=.not.der%mesh%parallel_in_domains)
      if (this%cube%parallel_in_domains .and. this%method == POISSON_FFT) then
        call mesh_cube_parallel_map_init(this%mesh_cube_map, der%mesh, this%cube)
      end if
    end if

    call poisson_kernel_init(this, mc%master_comm)

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

    case(POISSON_ISF)
      call poisson_isf_end(this%isf_solver)
      has_cube = .true.

    end select
    this%method = POISSON_NULL

    if (has_cube) then
      if (this%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_end(this%mesh_cube_map)
      end if
      call cube_end(this%cube)
    end if

    POP_SUB(poisson_end)
  end subroutine poisson_end

  !-----------------------------------------------------------------

  subroutine zpoisson_solve_real_and_imag_separately(this, pot, rho, all_nodes)
    type(poisson_t),      intent(in)    :: this
    CMPLX,                intent(inout) :: pot(:)  !< pot(mesh%np)
    CMPLX,                intent(in)    :: rho(:)  !< rho(mesh%np)
    logical, optional,    intent(in)    :: all_nodes

    FLOAT, allocatable :: aux1(:), aux2(:)
    type(derivatives_t), pointer :: der
    logical :: all_nodes_value

    der => this%der

    PUSH_SUB(zpoisson_solve_real_and_imag_separately)

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

    POP_SUB(zpoisson_solve_real_and_imag_separately)
  end subroutine zpoisson_solve_real_and_imag_separately

  !-----------------------------------------------------------------

  subroutine zpoisson_solve(this, pot, rho, all_nodes)
    type(poisson_t),      intent(in)    :: this
    CMPLX,                intent(inout) :: pot(:)  !< pot(mesh%np)
    CMPLX,                intent(in)    :: rho(:)  !< rho(mesh%np)
    logical, optional,    intent(in)    :: all_nodes

    logical :: all_nodes_value

    PUSH_SUB(zpoisson_solve)

    if(present(all_nodes)) then
      all_nodes_value = all_nodes
    else
      all_nodes_value = this%all_nodes_default
    end if

    ASSERT(this%method /= POISSON_NULL)

    call zpoisson_solve_real_and_imag_separately(this, pot, rho, all_nodes_value)

    POP_SUB(zpoisson_solve)
  end subroutine zpoisson_solve


  !-----------------------------------------------------------------

  subroutine poisson_solve_batch(this, potb, rhob, all_nodes)
    type(poisson_t),      intent(inout) :: this
    type(batch_t),        intent(inout) :: potb 
    type(batch_t),        intent(inout) :: rhob 
    logical, optional,    intent(in)    :: all_nodes

    integer :: ii

    PUSH_SUB(poisson_solve_batch)

    ASSERT(potb%nst_linear == rhob%nst_linear)
    ASSERT(batch_type(potb) == batch_type(rhob))

    if(batch_type(potb) == TYPE_FLOAT) then
      do ii = 1, potb%nst_linear
        call dpoisson_solve(this, potb%states_linear(ii)%dpsi, rhob%states_linear(ii)%dpsi, all_nodes)
      end do
    else
      do ii = 1, potb%nst_linear
        call zpoisson_solve(this, potb%states_linear(ii)%zpsi, rhob%states_linear(ii)%zpsi, all_nodes)
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
  subroutine dpoisson_solve(this, pot, rho, all_nodes)
    type(poisson_t),      intent(in)    :: this
    FLOAT,                intent(inout) :: pot(:) !< Local size of the \b potential vector. 
    FLOAT,                intent(inout) :: rho(:) !< Local size of the \b density (rho) vector.
    !> Is the Poisson solver allowed to utilise
    !! all nodes or only the domain nodes for
    !! its calculations? (Defaults to .true.)
    logical, optional,    intent(in)    :: all_nodes 
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
    if(present(all_nodes)) then
      all_nodes_value = all_nodes
    else
      all_nodes_value = this%all_nodes_default
    end if

    ASSERT(this%method /= POISSON_NULL)
      
    select case(this%method)
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
     
    end select

    POP_SUB(dpoisson_solve)
    call profiling_out(prof)
  end subroutine dpoisson_solve

  !-----------------------------------------------------------------
  !> This routine checks the Hartree solver selected in the input
  !! file by calculating numerically and analytically the Hartree
  !! potential originated by a Gaussian distribution of charge.
  !! This only makes sense for finite systems.
  subroutine poisson_test(mesh, repetitions)
    type(mesh_t), intent(in) :: mesh
    integer,      intent(in) :: repetitions

    FLOAT, allocatable :: rho(:), vh(:), vh_exact(:), rhop(:), xx(:, :)
    FLOAT :: alpha, beta, rr, delta, ralpha, hartree_nrg_num, &
         hartree_nrg_analyt, lcl_hartree_nrg 
    FLOAT :: total_charge
    integer :: ip, idir, ierr, iunit, nn, n_gaussians, itime

    PUSH_SUB(poisson_test)

    if(mesh%sb%dim == 1) then
      call messages_not_implemented('Poisson test for 1D case')
    end if

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
      do idir = 1, mesh%sb%dim
        xx(idir, nn) = M_ZERO 
      end do

      rr = sqrt(sum(xx(:, nn)*xx(:,nn)))
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin = xx(:, nn))
        rhop(ip) = beta*exp(-(rr/alpha)**2)
      end do

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
    do itime = 1, repetitions
      call dpoisson_solve(psolver, vh, rho)
    end do

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

    iterative = this%method == POISSON_CG .or. this%method == POISSON_CG_CORRECTED
  end function poisson_solver_is_iterative

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

  pure subroutine poisson_get_qpoint(this, qq)
    type(poisson_t), intent(in)  :: this
    FLOAT,           intent(out) :: qq(:)

    qq = this%qq
  end subroutine poisson_get_qpoint

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
  
  ! ---------------------------------------------------------
  subroutine poisson_kernel_init(this, all_nodes_comm)
    type(poisson_t),  intent(inout) :: this
    integer,          intent(in)    :: all_nodes_comm

    integer :: maxl, iter
    logical :: valid_solver

    PUSH_SUB(poisson_kernel_init)

    select case(this%method)
    case(POISSON_FFT, POISSON_CG, POISSON_CG_CORRECTED, POISSON_ISF)
      valid_solver = .true.
    case default
      valid_solver = .false.
    end select

    ASSERT(valid_solver)

    !%Variable PoissonSolverMaxMultipole
    !%Type integer
    !%Section Hamiltonian::Poisson
    !%Description
    !% Order of the multipolar expansion for boundary corrections. 
    !%
    !% The Poisson solvers <tt>multigrid</tt>, <tt>cg</tt>, and <tt>cg_corrected</tt>
    !% (and <tt>fft</tt> with <tt>PoissonFFTKernel = multipole_correction</tt>)
    !% do a multipolar expansion of the given
    !% charge density, such that <math>\rho = \rho_{multip.expansion}+\Delta
    !% \rho</math>. The Hartree potential due to the <math>\rho_{multip.expansion}</math> is
    !% calculated analytically, while the Hartree potential due to <math>\Delta \rho</math>
    !% is calculated with either a multigrid or cg solver.
    !% The order of the multipolar expansion is set by this variable.
    !%
    !% Default is 4 for <tt>PoissonSolver = cg_corrected</tt> and <tt>multigrid</tt>, and 2
    !% for <tt>fft</tt> with <tt>PoissonFFTKernel = multipole_correction</tt>.
    !%End

    !%Variable PoissonSolverMaxIter
    !%Type integer
    !%Section Hamiltonian::Poisson
    !%Default 400
    !%Description
    !% The maximum number of iterations for conjugate-gradient
    !% Poisson solvers.
    !%End

    !%Variable PoissonSolverThreshold
    !%Type float
    !%Section Hamiltonian::Poisson
    !%Default 1e-5
    !%Description
    !% The tolerance for the Poisson solution, used by the <tt>cg</tt>,
    !% <tt>cg_corrected</tt>.
    !%End

    !! This variable is disabled for the moment
    !!
    !!Variable PoissonSolverIncreaseBox
    !!Type logical
    !!Section Hamiltonian::Poisson
    !!Description
    !! (experimental) If the selected Poisson solver is
    !! <tt>cg_corrected</tt> the boundary conditions have to be
    !! calculated by means of performing a multipole
    !! expansion. Unfortunately, if the charge distribution is not
    !! contained in a simulation box of approximately spherical shape,
    !! the error can be quite large. Good cases are the spherical box,
    !! the parallelepiped when all dimensions are of similar magnitude,
    !! or the cylinder when the height is not too different to the
    !! diameter of the base. Bad cases are the rest, including the
    !! <tt>minimum</tt> box, when the geometry of the molecule is not
    !! compact enough.
    !!
    !! In order to cure this problem, the Hartree problem may be solved
    !! in an auxiliary simulation box, which will contain the original
    !! one, but which will be a sphere.  This implies some extra
    !! computational effort -- since the density and potential have to
    !! be transferred between boxes -- and extra memory consumption --
    !! since a new grid has to be stored.
    !!End

    select case(this%method)

    case(POISSON_CG)
      call parse_variable('PoissonSolverMaxMultipole', 4, maxl)
      write(message(1),'(a,i2)')'Info: Boundary conditions fixed up to L =',  maxl
      call messages_info(1)
      call parse_variable('PoissonSolverMaxIter', 400, iter)
      call parse_variable('PoissonSolverThreshold', CNST(1.0e-6), threshold)
      call poisson_corrections_init(this%corrector, maxl, this%der%mesh)
      call poisson_cg_init(threshold, iter)

    case(POISSON_CG_CORRECTED)
      call parse_variable('PoissonSolverMaxMultipole', 4, maxl)
      call parse_variable('PoissonSolverMaxIter', 400, iter)
      call parse_variable('PoissonSolverThreshold', CNST(1.0e-6), threshold)
      write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
      call messages_info(1)
      call poisson_corrections_init(this%corrector, maxl, this%der%mesh)
      call poisson_cg_init(threshold, iter)

    case(POISSON_ISF)
      call poisson_isf_init(this%isf_solver, this%der%mesh, this%cube, all_nodes_comm, init_world = this%all_nodes_default)

    case(POISSON_FFT)
      call poisson_fft_init(this%fft_solver, this%der%mesh, this%cube, this%kernel, qq = this%qq)
      ! soft parameter has no effect unless in 1D

      if (this%kernel == POISSON_FFT_KERNEL_CORRECTED) then
        call parse_variable('PoissonSolverMaxMultipole', 2, maxl)
        write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
        call messages_info(1)
        call poisson_corrections_init(this%corrector, maxl, this%der%mesh)
      end if

    end select

    POP_SUB(poisson_kernel_init)
  end subroutine poisson_kernel_init

  !-----------------------------------------------------------------
  subroutine poisson_kernel_reinit(this, qq)
    type(poisson_t), intent(inout) :: this
    FLOAT,           intent(in)    :: qq(:)

    PUSH_SUB(poisson_kernel_reinit)

    select case(this%method)
    case(POISSON_FFT)
      if(any(abs(this%qq(1:this%der%mesh%sb%periodic_dim) - qq(1:this%der%mesh%sb%periodic_dim)) > M_EPSILON)) then
        this%qq(1:this%der%mesh%sb%periodic_dim) = qq(1:this%der%mesh%sb%periodic_dim)
        call poisson_fft_end(this%fft_solver)
        call poisson_fft_init(this%fft_solver, this%der%mesh, this%cube, this%kernel, qq = this%qq)
      end if
    case default
      call messages_not_implemented("poisson_kernel_reinit with other methods than FFT")
    end select

    POP_SUB(poisson_kernel_reinit)
  end subroutine poisson_kernel_reinit

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
