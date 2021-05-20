!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012-2013 D. Strubbe
!! Copyright (C) 2017-2018 J. Flick, S. Ohlmann
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

module casida_oct_m
  use batch_oct_m
  use blacs_proc_grid_oct_m
  use calc_mode_par_oct_m
  use comm_oct_m
  use density_oct_m
#ifdef HAVE_ELPA
  use elpa
#endif
  use excited_states_oct_m
  use forces_oct_m
  use gauss_legendre_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use kpoints_oct_m
  use lalg_adv_oct_m
  use lda_u_oct_m
  use loct_oct_m
  use linear_response_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pblas_oct_m
  use pcm_oct_m
  use pert_oct_m
  use phonons_lr_oct_m
  use photon_mode_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scalapack_oct_m
  use sort_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use sternheimer_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                &
    casida_run,            &
    casida_run_init

  integer, parameter ::      &
    CASIDA_EPS_DIFF     = 1, &
    CASIDA_PETERSILKA   = 2, &
    CASIDA_TAMM_DANCOFF = 4, &
    CASIDA_VARIATIONAL  = 8, &
    CASIDA_CASIDA       = 16

  integer, parameter ::      &
    SOLVER_ELPA         = 1, &
    SOLVER_SCALAPACK    = 2

  type casida_t
    private
    integer :: type !< CASIDA_EPS_DIFF | CASIDA_PETERSILKA | CASIDA_TAMM_DANCOFF |
                    !< CASIDA_VARIATIONAL | CASIDA_CASIDA

    logical              :: states_are_real
    integer, allocatable :: n_occ(:)       !< number of occupied states
    integer, allocatable :: n_unocc(:)     !< number of unoccupied states
    integer              :: nst            !< total number of states
    integer              :: nik
    integer              :: space_dim         !< number of spatial dimensions
    integer              :: el_per_state
    character(len=80)    :: trandens
    character(len=80)    :: print_exst     !< excited states for which Casida coefficients will be printed
    FLOAT                :: weight_thresh  !< threshold for the Casida coefficients to be printed
    logical              :: triplet        !< use triplet kernel?
    logical              :: calc_forces    !< calculate excited-state forces
    logical              :: calc_forces_kernel    !< calculate excited-state forces with kernel
    logical              :: calc_forces_scf       !< calculate excited-state forces with SCF forces
    logical              :: herm_conj      !< use Hermitian conjugate of matrix
    type(restart_t)      :: restart_load
    type(restart_t)      :: restart_dump
    
    logical, allocatable :: is_included(:,:,:) !< (i, a, k) is in the basis?
    integer              :: n_pairs        !< number of pairs to take into account
    type(states_pair_t), allocatable :: pair(:)
    integer, allocatable :: index(:,:,:)   !< index(pair(j)%i, pair(j)%a, pair(j)%kk) = j
    integer, allocatable :: ind(:)         !< ordering in energy of solutions

    FLOAT, allocatable   :: dmat(:,:)      !< general-purpose matrix
    FLOAT, allocatable   :: dmat_save(:,:) !< to save mat when it gets turned into the eigenvectors
    CMPLX, allocatable   :: zmat(:,:)      !< general-purpose matrix
    CMPLX, allocatable   :: zmat_save(:,:) !< to save mat when it gets turned into the eigenvectors
    FLOAT, allocatable   :: w(:)           !< The excitation energies.
    FLOAT, allocatable   :: dtm(:, :)      !< The transition matrix elements (between the many-particle states)
    CMPLX, allocatable   :: ztm(:, :)      !< The transition matrix elements (between the many-particle states)
    FLOAT, allocatable   :: f(:)           !< The (dipole) strengths
    FLOAT, allocatable   :: s(:)           !< The diagonal part of the S-matrix

    FLOAT, allocatable   :: rho(:,:)       !< density
    FLOAT, allocatable   :: fxc(:,:,:)     !< derivative of xc potential
    FLOAT                :: kernel_lrc_alpha

    FLOAT, allocatable   :: dmat2(:,:)     !< matrix to diagonalize for forces
    CMPLX, allocatable   :: zmat2(:,:)     !< matrix to diagonalize for forces
    FLOAT, allocatable   :: dlr_hmat2(:,:) !< derivative of single-particle contribution to mat
    CMPLX, allocatable   :: zlr_hmat2(:,:) !< derivative of single-particle contribution to mat
    FLOAT, allocatable   :: forces(:,:,:)  !< excited-state forces
    FLOAT, allocatable   :: dw2(:)         !< perturbed excitation energies.
    FLOAT, allocatable   :: zw2(:)         !< perturbed excitation energies.

    ! variables for momentum-transfer-dependent calculation
    logical              :: qcalc
    FLOAT                :: qvector(MAX_DIM)
    FLOAT, allocatable   :: qf(:)
    FLOAT, allocatable   :: qf_avg(:)      !< Directionally averaged intensity
    integer              :: avg_order      !< Quadrature order for directional averaging (Gauss-Legendre scheme)

    logical              :: parallel_in_eh_pairs
    logical              :: parallel_in_domains
    logical              :: distributed_matrix
    logical              :: write_matrix
    integer              :: parallel_solver
    type(mpi_grp_t)      :: mpi_grp
    logical              :: fromScratch
    logical              :: has_photons
    integer              :: pt_nmodes
    type(photon_mode_t)  :: pt

    integer              :: n, nb_rows, nb_cols, block_size !< parallel matrix layout
    type(blacs_proc_grid_t) :: proc_grid       !< BLACS process grid type
    integer              :: desc(BLACS_DLEN)   !< descriptor for distributed matrix
    integer              :: darray             !< MPI IO type
  end type casida_t

  type casida_save_pot_t
    private
    integer :: qi                    !< previous mtxel calculated in K_term
    integer :: qa                    !< previous mtxel calculated in K_term
    integer :: qk                    !< previous mtxel calculated in K_term
    FLOAT, allocatable :: dpot(:)    !< previous exchange potential calculated in K_term
    CMPLX, allocatable :: zpot(:)    !< previous exchange potential calculated in K_term    
  end type casida_save_pot_t

contains

  subroutine casida_run_init()
    
    PUSH_SUB(casida_run_init)
    
    ! Pure 'other' parallelization is a bad idea. Trying to solve the Poisson equation separately on each node
    ! consumes excessive memory and time (easily more than is available). In principle, the line below would setup
    ! joint domain/other parallelization, but 'other' parallelization takes precedence, especially since
    ! multicomm_init does not know the actual problem size and uses a fictitious value of 10000, making it
    ! impossible to choose joint parallelization wisely, and generally resulting in a choice of only one domain
    ! group. FIXME! --DAS
    ! With the recent improvements, the 'other' parallelization over electron-hole
    ! pairs works quite well now. For smaller matrices, a combination
    ! seems to give the fastest run times. For larger matrices (more than a few
    ! thousand entries per dimension), CasidaDistributedMatrix is needed for
    ! the Casida matrix to fit into memory; this takes the cores from the
    ! 'other' strategy. For very large matrices (more than 100000), it is
    ! advisable to use only the 'other' strategy because the diagonalization
    ! uses most of the computation time.
    ! Thus you may want to enable this or a combination of other and domain to get
    ! better performance - STO

    call calc_mode_par_set_parallelization(P_STRATEGY_OTHER, default = .false.) ! enabled, but not default

    call calc_mode_par_unset_parallelization(P_STRATEGY_KPOINTS) ! disabled. FIXME: could be implemented.

    POP_SUB(casida_run_init)
  end subroutine casida_run_init

  ! ---------------------------------------------------------
  subroutine casida_run(system, from_scratch)
    class(*),        intent(inout) :: system
    logical,         intent(in)    :: from_scratch

    PUSH_SUB(casida_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = casida not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call casida_run_legacy(system, from_scratch)
    end select

    POP_SUB(casida_run)
  end subroutine casida_run

  ! ---------------------------------------------------------
  subroutine casida_run_legacy(sys, fromScratch)
    type(electrons_t), intent(inout) :: sys
    logical,           intent(in)    :: fromScratch

    type(casida_t) :: cas
    type(block_t) :: blk
    integer :: idir, theorylevel, iatom, ierr, default_int
    character(len=100) :: restart_filename
    type(profile_t), save :: prof
    logical :: is_frac_occ
    type(restart_t) :: gs_restart

    PUSH_SUB(casida_run_legacy)
    call profiling_in(prof, 'CASIDA')

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    if (sys%space%is_periodic()) then
      message(1) = "Casida oscillator strengths will be incorrect in periodic systems."
      call messages_warning(1)
    end if

    if(kpoints_number(sys%kpoints) > 1) then
      ! Hartree matrix elements may not be correct, not tested anyway. --DAS
      call messages_not_implemented("Casida with k-points")
    end if
    if (family_is_mgga_with_exc(sys%hm%xc)) then
      call messages_not_implemented("Casida with MGGA and non-local terms")
    end if
    if(sys%hm%lda_u_level /= DFT_U_NONE) then
      call messages_not_implemented("Casida with DFT+U")
    end if
    if(sys%hm%theory_level == HARTREE_FOCK) then
      call messages_not_implemented("Casida for Hartree-Fock")
    end if
    if(sys%hm%theory_level == GENERALIZED_KOHN_SHAM_DFT) then
      call messages_not_implemented("Casida for generalized Kohn-Sham")
    end if

    message(1) = 'Info: Starting Casida linear-response calculation.'
    call messages_info(1)

    call restart_init(gs_restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
    if(ierr == 0) then
      call states_elec_look_and_load(gs_restart, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints)
      call restart_end(gs_restart)
    else
      message(1) = "Previous gs calculation is required."
      call messages_fatal(1)
    end if

    cas%el_per_state = sys%st%smear%el_per_state
    cas%nst = sys%st%nst
    cas%nik = sys%st%d%nik
    cas%space_dim = sys%space%dim
    SAFE_ALLOCATE(cas%n_occ(1:sys%st%d%nik))
    SAFE_ALLOCATE(cas%n_unocc(1:sys%st%d%nik))

    call states_elec_count_pairs(sys%st, sys%namespace, cas%n_pairs, cas%n_occ, cas%n_unocc, cas%is_included, is_frac_occ)
    if(is_frac_occ) then
      call messages_not_implemented("Casida with partial occupations")
      ! Formulas are in Casida 1995 reference. The occupations are not used at all here currently.
    end if

    select case(sys%st%d%ispin)
    case(UNPOLARIZED, SPINORS)
      write(message(1),'(a,i4,a)') "Info: Found", cas%n_occ(1), " occupied states."
      write(message(2),'(a,i4,a)') "Info: Found", cas%n_unocc(1), " unoccupied states."
      call messages_info(2)
    case(SPIN_POLARIZED)
      write(message(1),'(a,i4,a)') "Info: Found", cas%n_occ(1), " occupied states with spin up."
      write(message(2),'(a,i4,a)') "Info: Found", cas%n_unocc(1), " unoccupied states with spin up."
      write(message(3),'(a,i4,a)') "Info: Found", cas%n_occ(2), " occupied states with spin down."
      write(message(4),'(a,i4,a)') "Info: Found", cas%n_unocc(2), " unoccupied states with spin down."
      call messages_info(4)
    end select


    ! setup Hamiltonian, without recalculating eigenvalues (use the ones from the restart information)
    message(1) = 'Info: Setting up Hamiltonian.'
    call messages_info(1)
    call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, calc_eigenval=.false.)

    !%Variable CasidaTheoryLevel
    !%Type flag
    !%Section Linear Response::Casida
    !%Default <tt>eps_diff + petersilka + lrtddft_casida</tt>
    !%Description
    !% Choose which electron-hole matrix-based theory levels to use in calculating excitation energies.
    !% More than one may be used to take advantage of the significant commonality between the calculations.
    !% <tt>variational</tt> and <tt>lrttdft_casida</tt> are not usable with complex wavefunctions.
    !% Note the restart data saved by each theory level is compatible with all the others.
    !%Option eps_diff 1
    !% Difference of eigenvalues, <i>i.e.</i> independent-particle approximation.
    !%Option petersilka 2
    !% The Petersilka approximation uses only elements of the Tamm-Dancoff matrix between degenerate
    !% transitions (if no degeneracy, this is just the diagonal elements). Also called the "single-pole" approximation.
    !% This is acceptable if there is little mixing between single-particle transitions.
    !% Ref: M Petersilka, UJ Gossmann, and EKU Gross, <i>Phys. Rev. Lett.</i> <b>76</b>, 1212 (1996);
    !% T Grabo, M Petersilka,and EKU Gross, <i>Theochem</i> <b>501-502</b> 353 (2000).
    !%Option tamm_dancoff 4
    !% The Tamm-Dancoff approximation uses only occupied-unoccupied transitions and not
    !% unoccupied-occupied transitions.
    !% Ref: S Hirata and M Head-Gordon, <i>Chem. Phys. Lett.</i> <b>314</b>, 291 (1999).
    !%Option variational 8
    !% Second-order constrained variational theory CV(2)-DFT. Only applies to real wavefunctions.
    !% Ref: T Ziegler, M Seth, M Krykunov, J Autschbach, and F Wang,
    !% <i>J. Chem. Phys.</i> <b>130</b>, 154102 (2009).
    !%Option lrtddft_casida 16
    !% The full Casida method. Only applies to real wavefunctions.
    !% Ref: C Jamorski, ME Casida, and DR Salahub, <i>J. Chem. Phys.</i> <b>104</b>, 5134 (1996)
    !% and ME Casida, "Time-dependent density functional response theory for molecules,"
    !% in <i>Recent Advances in Density Functional Methods</i>, edited by DE Chong, vol. 1
    !% of <i>Recent Advances in Computational Chemistry</i>, pp. 155-192 (World Scientific,
    !% Singapore, 1995).
    !%End

    call parse_variable(sys%namespace, 'CasidaTheoryLevel', CASIDA_EPS_DIFF + CASIDA_PETERSILKA + CASIDA_CASIDA, theorylevel)

    if (states_are_complex(sys%st)) then
      if((bitand(theorylevel, CASIDA_VARIATIONAL) /= 0 &
        .or. bitand(theorylevel, CASIDA_CASIDA) /= 0)) then
        message(1) = "Variational and full Casida theory levels do not apply to complex wavefunctions."
        call messages_fatal(1, only_root_writes = .true.)
        ! see section II.D of CV(2) paper regarding this assumption. Would be Eq. 30 with complex wfns.
      end if
    end if

    ! This variable is documented in xc_oep_init.
    call parse_variable(sys%namespace, 'EnablePhotons', .false., cas%has_photons)
    cas%pt_nmodes = 0
    if (cas%has_photons) then
      if(cas%has_photons) call messages_experimental('EnablePhotons = yes')
      call photon_mode_init(cas%pt, sys%namespace, sys%gr%mesh, sys%space%dim, sys%st%qtot)
      write(message(1), '(a,i7,a)') 'INFO: Solving Casida equation with ', &
        cas%pt%nmodes, ' photon modes.'
      write(message(2), '(a)') 'as described in ACS Photonics 2019, 6, 11, 2757-2778.'
      call messages_info(2)
      cas%pt_nmodes = cas%pt%nmodes
    end if

    !%Variable CasidaTransitionDensities
    !%Type string
    !%Section Linear Response::Casida
    !%Default write none
    !%Description
    !% Specifies which transition densities are to be calculated and written down. The
    !% transition density for the many-body state <i>n</i> will be written to a file called
    !% <tt>rho_0n</tt> prefixed by the theory level. Format is set by <tt>OutputFormat</tt>.
    !%
    !% This variable is a string in list form, <i>i.e.</i> expressions such as "1,2-5,8-15" are
    !% valid.
    !%End
    call parse_variable(sys%namespace, 'CasidaTransitionDensities', "0", cas%trandens)

    if (cas%trandens /= "0") then 
      call io_function_read_what_how_when(sys%namespace, sys%space, sys%outp%what,&
        sys%outp%how, sys%outp%output_interval)
    end if

    !%Variable CasidaMomentumTransfer
    !%Type block
    !%Section Linear Response::Casida
    !%Default 0
    !%Description
    !% Momentum-transfer vector for the calculation of the dynamic structure
    !% factor. When this variable is set, the transition rates are determined
    !% using an exponential operator instead of the normal dipole one.
    !%End

    if(parse_block(sys%namespace, 'CasidaMomentumTransfer', blk)==0) then
      do idir = 1, cas%space_dim
        call parse_block_float(blk, 0, idir - 1, cas%qvector(idir))
        cas%qvector(idir) = units_to_atomic(unit_one / units_inp%length, cas%qvector(idir))
      end do
      call parse_block_end(blk)
      call messages_experimental("IXS/EELS transition rate calculation")
      message(1) = "Info: Calculating IXS/EELS transition rates."
      call messages_info(1)
      cas%qcalc = .true.

      !%Variable CasidaQuadratureOrder
      !%Type integer
      !%Section Linear Response::Casida
      !%Default 5
      !%Description
      !% Only applies if <tt>CasidaMomentumTransfer</tt> is nonzero.
      !% Directionally averaged dynamic structure factor is calculated by
      !% averaging over the results from a set of <math>\vec{q}</math>-vectors. The vectors
      !% are generated using Gauss-Legendre quadrature scheme [see <i>e.g.</i>
      !% K. Atkinson, <i>J. Austral. Math. Soc.</i> <b>23</b>, 332 (1982)], and this
      !% variable determines the order of the scheme.
      !%End
      call parse_variable(sys%namespace, 'CasidaQuadratureOrder', 5, cas%avg_order)
    else
      cas%qvector(:) = M_ZERO
      cas%qcalc = .false.
    end if

    !%Variable CasidaCalcTriplet
    !%Type logical
    !%Section Linear Response::Casida
    !%Default false
    !%Description
    !% For a non-spin-polarized ground state, singlet or triplet excitations can be calculated
    !% using different matrix elements. Default is to calculate singlets. This variable has no
    !% effect for a spin-polarized calculation.
    !%End
    if(sys%st%d%ispin == UNPOLARIZED) then
      call parse_variable(sys%namespace, 'CasidaCalcTriplet', .false., cas%triplet)
    else
      cas%triplet = .false.
    end if

    if(cas%triplet) then
      call messages_experimental("Casida triplet calculation")
      message(1) = "Info: Using triplet kernel. Oscillator strengths will be for spin magnetic-dipole field."
      call messages_info(1)
    end if

    !%Variable CasidaHermitianConjugate
    !%Type logical
    !%Section Linear Response::Casida
    !%Default false
    !%Description
    !% The Casida matrix is Hermitian, so it should not matter whether we calculate the upper or
    !% lower diagonal. Numerical issues may cause small differences however. Use this variable to
    !% calculate the Hermitian conjugate of the usual matrix, for testing.
    !%End
    call parse_variable(sys%namespace, 'CasidaHermitianConjugate', .false., cas%herm_conj)

    !%Variable CasidaDistributedMatrix
    !%Type logical
    !%Section Linear Response::Casida
    !%Default false
    !%Description
    !% Large matrices with more than a few thousand rows and columns usually do
    !% not fit into the memory of one processor anymore. With this option, the
    !% Casida matrix is distributed in block-cyclic fashion over all cores in the
    !% ParOther group. The diagonalization is done in parallel using ScaLAPACK
    !% or ELPA, if available. For very large matrices (>100000), only the
    !% ParOther strategy should be used because the diagonalization dominates
    !% the run time of the computation.
    !%End
    call parse_variable(sys%namespace, 'CasidaDistributedMatrix', .false., cas%distributed_matrix)
#ifndef HAVE_SCALAPACK
    if(cas%distributed_matrix) then
      message(1) = "ScaLAPACK layout requested, but code not compiled with ScaLAPACK"
      call messages_fatal(1)
    end if
#endif

    !%Variable CasidaWriteDistributedMatrix
    !%Type logical
    !%Section Linear Response::Casida
    !%Default false
    !%Description
    !% Set to true to write out the full distributed Casida matrix to a file
    !% using MPI-IO.
    !%End
    call parse_variable(sys%namespace, 'CasidaWriteDistributedMatrix', .false., cas%write_matrix)
    if(.not.cas%distributed_matrix .and. cas%write_matrix) then
      message(1) = "CasidaWriteDistributedMatrix con only be used with CasidaDistributedMatrix"
      call messages_fatal(1)
    end if

    !%Variable CasidaParallelEigensolver
    !%Type integer
    !%Section Linear Response::Casida
    !%Description
    !% Choose library to use for solving the parallel eigenproblem
    !% of the Casida problem. This options is only relevant if a
    !% distributed matrix is used (CasidaDistributedMatrix=true).
    !% By default, elpa is chosen if available.
    !%Option casida_elpa 1
    !% Use ELPA library as parallel eigensolver
    !%Option casida_scalapack 2
    !% Use Scalapack as parallel eigensolver
    !%End
#ifdef HAVE_ELPA
    default_int = SOLVER_ELPA
#else
    default_int = SOLVER_SCALAPACK
#endif
    call parse_variable(sys%namespace, 'CasidaParallelEigensolver', default_int, cas%parallel_solver)
    if(.not.varinfo_valid_option('CasidaParallelEigensolver', cas%parallel_solver)) then
      call messages_input_error(sys%namespace, 'CasidaParallelEigensolver')
    end if
#ifndef HAVE_ELPA
    if(cas%distributed_matrix .and. cas%parallel_solver == SOLVER_ELPA) then
      message(1) = "ELPA solver requested, but code not compiled with ELPA"
      call messages_fatal(1)
    end if
#endif

    !%Variable CasidaPrintExcitations
    !%Type string
    !%Section Linear Response::Casida
    !%Default write all
    !%Description
    !% Specifies which excitations are written at the end of the calculation. 
    !%
    !% This variable is a string in list form, <i>i.e.</i> expressions such as "1,2-5,8-15" are
    !% valid.
    !%End
    call parse_variable(sys%namespace, 'CasidaPrintExcitations', "all", cas%print_exst)
    if(cas%distributed_matrix) then
      ! do not print excited states -> too many files generated!
      cas%print_exst = "none"
      message(1) = "Using ScaLAPACK layout, thus disabling output of excited states."
      message(2) = "This options creates too many files for large Casida matrices."
      call messages_info(2)
    end if

    !%Variable CasidaWeightThreshold
    !%Type float
    !%Section Linear Response::Casida
    !%Default -1.
    !%Description
    !% Specifies the threshold value for which the individual excitations are printed. 
    !% i.e. juste-h pairs with weight larger than this threshold will be printed. 
    !% 
    !% If a negative value (default) is set, all coefficients will be printed.
    !% For many case, a 0.01 value is a valid option.
    !%End
    call parse_variable(sys%namespace, 'CasidaWeightThreshold', -M_ONE, cas%weight_thresh)
    if (cas%weight_thresh > M_ONE) then
      message(1) = 'Casida coefficients have values between 0 and 1'
      message(2) = 'Threshold values reset to default value'
      call messages_warning(2)
      cas%weight_thresh = -M_ONE
    end if

    !%Variable CasidaCalcForces
    !%Type logical
    !%Section Linear Response::Casida
    !%Default false
    !%Description
    !% (Experimental) Enable calculation of excited-state forces. Requires previous <tt>vib_modes</tt> calculation.
    !%End
    call parse_variable(sys%namespace, 'CasidaCalcForces', .false., cas%calc_forces)
    if(cas%calc_forces) then
      call messages_experimental("Excited-state forces calculation")

      !%Variable CasidaCalcForcesKernel
      !%Type logical
      !%Section Linear Response::Casida
      !%Default true
      !%Description
      !% If false, the derivative of the kernel will not be included in the excited-state force calculation.
      !%End
      call parse_variable(sys%namespace, 'CasidaCalcForcesKernel', .true., cas%calc_forces_kernel)

      !%Variable CasidaCalcForcesSCF
      !%Type logical
      !%Section Linear Response::Casida
      !%Default false
      !%Description
      !% If true, the ground-state forces will be included in the excited-state forces, so they are total forces.
      !% If false, the excited-state forces that are produced are only the gradients of the excitation energy.
      !%End
      call parse_variable(sys%namespace, 'CasidaCalcForcesSCF', .false., cas%calc_forces_scf)

      if(cas%distributed_matrix) then
        message(1) = "Info: Forces calculation not compatible with ScaLAPACK layout."
        message(2) = "Using normal layout."
        call messages_info(2)
        cas%distributed_matrix = .false.
      end if
    end if

    ! Initialize structure
    call casida_type_init(cas, sys)

    cas%fromScratch = fromScratch

    if(cas%fromScratch) then ! remove old restart files
      if(cas%triplet) then
        call restart_rm(cas%restart_dump, 'kernel_triplet')
      else
        call restart_rm(cas%restart_dump, 'kernel')
      end if

      if(cas%calc_forces) then
        do iatom = 1, sys%ions%natoms
          do idir = 1, cas%space_dim
            write(restart_filename,'(a,i6.6,a,i1)') 'lr_kernel_', iatom, '_', idir
            if (cas%triplet) restart_filename = trim(restart_filename)//'_triplet'
            call restart_rm(cas%restart_dump, restart_filename)

            write(restart_filename,'(a,i6.6,a,i1)') 'lr_hmat1_', iatom, '_', idir
            call restart_rm(cas%restart_dump, restart_filename)
          end do
        end do
      end if
    end if

    ! First, print the differences between KS eigenvalues (first approximation to the excitation energies).
    if(bitand(theorylevel, CASIDA_EPS_DIFF) /= 0) then
      message(1) = "Info: Approximating resonance energies through KS eigenvalue differences"
      call messages_info(1)
      cas%type = CASIDA_EPS_DIFF
      call casida_work(sys, cas)
    end if

    if (sys%st%d%ispin /= SPINORS) then

      if(bitand(theorylevel, CASIDA_TAMM_DANCOFF) /= 0) then
        call messages_experimental("Tamm-Dancoff calculation")
        message(1) = "Info: Calculating matrix elements in the Tamm-Dancoff approximation"
        call messages_info(1)
        cas%type = CASIDA_TAMM_DANCOFF
        call casida_work(sys, cas)
      end if

      if(bitand(theorylevel, CASIDA_VARIATIONAL) /= 0) then
        call messages_experimental("CV(2)-DFT calculation")
        message(1) = "Info: Calculating matrix elements with the CV(2)-DFT theory"
        call messages_info(1)
        cas%type = CASIDA_VARIATIONAL
        call casida_work(sys, cas)
      end if

      if(bitand(theorylevel, CASIDA_CASIDA) /= 0) then
        message(1) = "Info: Calculating matrix elements with the full Casida method"
        call messages_info(1)
        cas%type = CASIDA_CASIDA
        call casida_work(sys, cas)
      end if

      ! Doing this first, if doing the others later, takes longer, because we would use
      ! each Poisson solution for only one matrix element instead of a whole column.
      if(bitand(theorylevel, CASIDA_PETERSILKA) /= 0) then
        message(1) = "Info: Calculating resonance energies via the Petersilka approximation"
        call messages_info(1)
        cas%type = CASIDA_PETERSILKA
        call casida_work(sys, cas)
      end if

    end if

    call casida_type_end(cas)

    call profiling_out(prof)
    POP_SUB(casida_run_legacy)
  end subroutine casida_run_legacy

  ! ---------------------------------------------------------
  !> allocates stuff, and constructs the arrays pair_i and pair_j
  subroutine casida_type_init(cas, sys)
    type(casida_t),      intent(inout) :: cas
    type(electrons_t),   intent(in)    :: sys

    integer :: ist, ast, jpair, ik, ierr
#ifdef HAVE_SCALAPACK
    integer :: np, np_rows, np_cols, ii, info
#endif

    PUSH_SUB(casida_type_init)

    cas%kernel_lrc_alpha = sys%ks%xc%kernel_lrc_alpha
    cas%states_are_real = states_are_real(sys%st)
    if(cas%distributed_matrix .and. .not. cas%states_are_real) then
      call messages_not_implemented("Complex wavefunctions with ScaLAPACK layout")
    end if

    write(message(1), '(a,i9)') "Number of occupied-unoccupied pairs: ", cas%n_pairs
    call messages_info(1)

    if(cas%n_pairs < 1) then
      message(1) = "No Casida pairs -- maybe there are no unoccupied states?"
      call messages_fatal(1, only_root_writes = .true.)
    end if

    if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

    ! now let us take care of initializing the parallel stuff
    cas%parallel_in_eh_pairs = multicomm_strategy_is_parallel(sys%mc, P_STRATEGY_OTHER)
    if(cas%parallel_in_eh_pairs) then
      call mpi_grp_init(cas%mpi_grp, sys%mc%group_comm(P_STRATEGY_OTHER))
    else
      call mpi_grp_init(cas%mpi_grp, -1)
    end if
    cas%parallel_in_domains = multicomm_strategy_is_parallel(sys%mc, P_STRATEGY_DOMAINS)

    if(cas%distributed_matrix .and. .not. cas%parallel_in_eh_pairs) then
      message(1) = "ScaLAPACK layout requested, but 'Other' parallelization strategy not available."
      message(2) = "Please set ParOther to use the ScaLAPACK layout."
      message(3) = "Continuing without ScaLAPACK layout."
      call messages_info(3)
      cas%distributed_matrix = .false.
    end if

    ! dimension of matrix
    cas%n = cas%n_pairs + cas%pt_nmodes

    ! initialize block-cyclic matrix
    if(cas%distributed_matrix) then
#ifdef HAVE_SCALAPACK
      ! processor layout: always use more processors for rows, this leads to
      ! better load balancing when computing the matrix elements
      np = cas%mpi_grp%size
      np_cols = 1
      if(np > 3) then
        do ii = floor(sqrt(real(np))), 2, -1
          if(mod(np, ii) == 0) then
            np_cols = ii
            exit
          end if
        end do
      end if
      np_rows = np / np_cols

      ! recommended block size: 64, take smaller value for smaller matrices for
      ! better load balancing
      cas%block_size = min(64, cas%n / np_rows)
      ! limit to a minimum block size of 5 for diagonalization efficiency
      cas%block_size = max(5, cas%block_size)
      write(message(1), '(A,I5,A,I5,A,I5,A)') 'Parallel layout: using block size of ',&
        cas%block_size, ' and a processor grid with ', np_rows, 'x', np_cols, &
        ' processors (rows x cols)'
      call messages_info(1)

      call blacs_proc_grid_init(cas%proc_grid, cas%mpi_grp, procdim = (/np_rows, np_cols/))

      ! get size of local matrices
      cas%nb_rows = numroc(cas%n, cas%block_size, cas%proc_grid%myrow, 0, cas%proc_grid%nprow)
      cas%nb_cols = numroc(cas%n, cas%block_size, cas%proc_grid%mycol, 0, cas%proc_grid%npcol)

      ! get ScaLAPACK descriptor
      call descinit(cas%desc(1), cas%n, cas%n, cas%block_size, cas%block_size, 0, 0, &
        cas%proc_grid%context, cas%nb_rows, info)
#endif
    else
      ! set to full size
      cas%nb_rows = cas%n
      cas%nb_cols = cas%n
    end if


    ! allocate stuff
    SAFE_ALLOCATE(cas%pair(1:cas%n))
    if(cas%states_are_real) then
      SAFE_ALLOCATE( cas%dmat(1:cas%nb_rows, 1:cas%nb_cols))
      SAFE_ALLOCATE(  cas%dtm(1:cas%n, 1:cas%space_dim))
    else
      ! caution: ScaLAPACK layout not yet tested for complex wavefunctions!
      SAFE_ALLOCATE( cas%zmat(1:cas%nb_rows, 1:cas%nb_cols))
      SAFE_ALLOCATE(  cas%ztm(1:cas%n, 1:cas%space_dim))
    end if
    SAFE_ALLOCATE(   cas%f(1:cas%n))
    SAFE_ALLOCATE(   cas%s(1:cas%n_pairs))
    SAFE_ALLOCATE(   cas%w(1:cas%n))
    SAFE_ALLOCATE(cas%index(1:maxval(cas%n_occ), cas%nst - maxval(cas%n_unocc) + 1:cas%nst, cas%nik))
    SAFE_ALLOCATE( cas%ind(1:cas%n))

    if(cas%calc_forces) then
      if(cas%states_are_real) then
        SAFE_ALLOCATE(cas%dmat_save(1:cas%n_pairs, 1:cas%n_pairs))
      else
        SAFE_ALLOCATE(cas%zmat_save(1:cas%n_pairs, 1:cas%n_pairs))
      end if
      SAFE_ALLOCATE(cas%forces(1:cas%space_dim, 1:sys%ions%natoms, 1:cas%n_pairs))
    end if

    if(cas%qcalc) then
      SAFE_ALLOCATE( cas%qf    (1:cas%n_pairs))
      SAFE_ALLOCATE( cas%qf_avg(1:cas%n_pairs))
    end if

    cas%index(:,:,:) = 0

    ! create pairs
    jpair = 1
    do ik = 1, cas%nik
      do ast = cas%n_occ(ik) + 1, cas%nst
        do ist = 1, cas%n_occ(ik)
          if(cas%is_included(ist, ast, ik)) then
            cas%index(ist, ast, ik) = jpair
            cas%pair(jpair)%i = ist
            cas%pair(jpair)%a = ast
            cas%pair(jpair)%kk = ik
            jpair = jpair + 1
          end if
        end do
      end do
    end do

    if (cas%has_photons) then
    ! create pairs for photon modes (negative number refers to photonic excitation)
      do ik = 1, cas%pt_nmodes
         cas%pair(cas%n_pairs + ik)%i = 1
         cas%pair(cas%n_pairs + ik)%a = -ik
         cas%pair(cas%n_pairs + ik)%kk = -ik
      end do
    end if

    SAFE_DEALLOCATE_A(cas%is_included)

    call restart_init(cas%restart_dump, sys%namespace, RESTART_CASIDA, RESTART_TYPE_DUMP, sys%mc, ierr)
    call restart_init(cas%restart_load, sys%namespace, RESTART_CASIDA, RESTART_TYPE_LOAD, sys%mc, ierr)

    POP_SUB(casida_type_init)
  end subroutine casida_type_init


  ! ---------------------------------------------------------
  subroutine casida_type_end(cas)
    type(casida_t), intent(inout) :: cas

    PUSH_SUB(casida_type_end)

    ASSERT(allocated(cas%pair))
    SAFE_DEALLOCATE_A(cas%pair)
    SAFE_DEALLOCATE_A(cas%index)
    if(cas%states_are_real) then
      SAFE_DEALLOCATE_A(cas%dmat)
      SAFE_DEALLOCATE_A(cas%dtm)
    else
      SAFE_DEALLOCATE_A(cas%zmat)
      SAFE_DEALLOCATE_A(cas%ztm)
    end if
    SAFE_DEALLOCATE_A(cas%s)
    SAFE_DEALLOCATE_A(cas%f)
    SAFE_DEALLOCATE_A(cas%w)
    SAFE_DEALLOCATE_A(cas%ind)

    if(cas%qcalc) then
      SAFE_DEALLOCATE_A(cas%qf)
      SAFE_DEALLOCATE_A(cas%qf_avg)
    end if

    SAFE_DEALLOCATE_A(cas%n_occ)
    SAFE_DEALLOCATE_A(cas%n_unocc)

    if(cas%calc_forces) then
      if(cas%states_are_real) then
        SAFE_DEALLOCATE_A(cas%dmat_save)
      else
        SAFE_DEALLOCATE_A(cas%zmat_save)
      end if
      SAFE_DEALLOCATE_A(cas%forces)
    end if

    call restart_end(cas%restart_dump)
    call restart_end(cas%restart_load)

    if (cas%has_photons) then
      call photon_mode_end(cas%pt)
    end if
    if(cas%distributed_matrix) then
#ifdef HAVE_SCALAPACK
      call blacs_proc_grid_end(cas%proc_grid)
#endif
    end if

    POP_SUB(casida_type_end)
  end subroutine casida_type_end


  ! ---------------------------------------------------------
  !> this subroutine calculates electronic excitation energies using
  !! the matrix formulation of M. Petersilka, or of M. Casida
  subroutine casida_work(sys, cas)
    type(electrons_t),   target, intent(inout) :: sys
    type(casida_t),              intent(inout) :: cas

    type(states_elec_t), pointer :: st
    type(mesh_t),   pointer :: mesh

    FLOAT, allocatable :: rho_spin(:, :)
    FLOAT, allocatable :: fxc_spin(:,:,:)
    character(len=100) :: restart_filename

    PUSH_SUB(casida_work)

    ! sanity checks
    ASSERT(cas%type >= CASIDA_EPS_DIFF .and. cas%type <= CASIDA_CASIDA)

    ! some shortcuts
    st => sys%st
    mesh => sys%gr%mesh

    ! initialize stuff
    if(cas%states_are_real) then
      cas%dmat = M_ZERO
      cas%dtm  = M_ZERO
    else
      cas%zmat = M_ZERO
      cas%ztm  = M_ZERO
    end if
    cas%f   = M_ZERO
    cas%w   = M_ZERO
    cas%s   = M_ZERO
    if(cas%qcalc) then
      cas%qf     = M_ZERO
      cas%qf_avg = M_ZERO
    end if

    if (cas%type /= CASIDA_EPS_DIFF .or. cas%calc_forces) then
      ! We calculate here the kernel, since it will be needed later.
      SAFE_ALLOCATE(cas%rho(1:mesh%np, 1:st%d%nspin))
      SAFE_ALLOCATE(cas%fxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
      cas%fxc = M_ZERO

      call states_elec_total_density(st, mesh, cas%rho)
      if(cas%triplet) then
        SAFE_ALLOCATE(rho_spin(1:mesh%np, 1:2))
        SAFE_ALLOCATE(fxc_spin(1:mesh%np, 1:2, 1:2))

        fxc_spin = M_ZERO
        rho_spin(:, 1) = M_HALF * cas%rho(:, 1)
        rho_spin(:, 2) = M_HALF * cas%rho(:, 1)

        call xc_get_fxc(sys%ks%xc, mesh, sys%namespace, rho_spin, SPIN_POLARIZED, fxc_spin)
        cas%fxc(:, 1, 1) = M_HALF * (fxc_spin(:, 1, 1) - fxc_spin(:, 1, 2))

        SAFE_DEALLOCATE_A(rho_spin)
        SAFE_DEALLOCATE_A(fxc_spin)
      else
        call xc_get_fxc(sys%ks%xc, mesh, sys%namespace, cas%rho, st%d%ispin, cas%fxc)
      end if

      if (sys%ks%sic_type == SIC_ADSIC) then
        call fxc_add_adsic(sys%namespace, sys%ks, st, mesh, cas)
      end if

    end if

    restart_filename = 'kernel'
    if(cas%triplet) restart_filename = trim(restart_filename)//'_triplet'

    select case(cas%type)
    case(CASIDA_EPS_DIFF)
      call solve_eps_diff()
    case(CASIDA_TAMM_DANCOFF,CASIDA_VARIATIONAL,CASIDA_CASIDA,CASIDA_PETERSILKA)
      if(cas%states_are_real) then
        call dcasida_get_matrix(cas, sys%hm, st, sys%ks, mesh, cas%dmat, cas%fxc, restart_filename)
        call dcasida_solve(cas, sys)
      else
        call zcasida_get_matrix(cas, sys%hm, st, sys%ks, mesh, cas%zmat, cas%fxc, restart_filename)
        call zcasida_solve(cas, sys)
      end if
    end select

    ! compute oscillator strengths on all processes for the ScaLAPACK layout
    if (mpi_grp_is_root(cas%mpi_grp) .or. cas%distributed_matrix) then
      if(cas%states_are_real) then
        call doscillator_strengths(cas, mesh, st)
      else
        call zoscillator_strengths(cas, mesh, st)
      end if
    end if

    if(cas%calc_forces) then
      if(cas%states_are_real) then
        call dcasida_forces(cas, sys, mesh, st)
      else
        call zcasida_forces(cas, sys, mesh, st)
      end if
    end if

    if(cas%states_are_real) then
      call dcasida_write(cas, sys)
    else
      call zcasida_write(cas, sys)
    end if

    ! clean up
    if(cas%type /= CASIDA_EPS_DIFF .or. cas%calc_forces) then
      SAFE_DEALLOCATE_A(cas%fxc)
      SAFE_DEALLOCATE_A(cas%rho)
    end if

    POP_SUB(casida_work)

  contains

    ! ---------------------------------------------------------
    subroutine solve_eps_diff

      integer :: ia
      FLOAT, allocatable :: w(:)

      PUSH_SUB(casida_work.solve_eps_diff)

      ! initialize progress bar
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, cas%n_pairs)

      do ia = 1, cas%n_pairs
        cas%w(ia) = st%eigenval(cas%pair(ia)%a, cas%pair(ia)%kk) - &
                    st%eigenval(cas%pair(ia)%i, cas%pair(ia)%kk)
        if(cas%w(ia) < -M_EPSILON) then
          message(1) = "There is a negative unocc-occ KS eigenvalue difference for"
          write(message(2),'("states ",I5," and ",I5," of k-point ",I5,".")') cas%pair(ia)%i, cas%pair(ia)%a, cas%pair(ia)%kk
          message(3) = "This indicates an inconsistency between gs, unocc, and/or casida calculations."
          call messages_fatal(3, only_root_writes = .true.)
        end if
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ia, cas%n_pairs)
      end do

      SAFE_ALLOCATE(w(1:cas%n_pairs))
      w = cas%w
      call sort(w, cas%ind)
      SAFE_DEALLOCATE_A(w)

      if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

      POP_SUB(casida_work.solve_eps_diff)
    end subroutine solve_eps_diff

    ! ---------------------------------------------------------
    subroutine fxc_add_adsic(namespace, ks, st, mesh, cas)
      type(namespace_t),   intent(in)    :: namespace
      type(v_ks_t),        intent(in)    :: ks
      type(states_elec_t), intent(in)    :: st
      type(mesh_t),        intent(in)    :: mesh
      type(casida_t),      intent(inout) :: cas

      FLOAT, allocatable :: rho(:, :)
      FLOAT, allocatable :: fxc_sic(:,:,:)

      PUSH_SUB(casida_work.fxc_add_adsic)

      !Check spin and triplets
      if (st%d%ispin /= UNPOLARIZED) then
        message(1) = "Casida calculation with ADSIC not implemented for spin-polarized calculations."
        call messages_fatal(1, namespace=namespace)
      end if
      if (cas%triplet) then
        message(1) = "Casida calculation with ADSIC not implemented for triplet excitations."
        call messages_fatal(1, namespace=namespace)
      end if

      SAFE_ALLOCATE(fxc_sic(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
      SAFE_ALLOCATE(rho(1:mesh%np, 1:st%d%nspin))
      fxc_sic = M_ZERO
      rho = cas%rho/st%qtot

      call xc_get_fxc(ks%xc, mesh, namespace, rho, 1, fxc_sic)

      cas%fxc = cas%fxc - fxc_sic/st%qtot

      SAFE_DEALLOCATE_A(rho)
      SAFE_DEALLOCATE_A(fxc_sic)

      POP_SUB(casida_work.fxc_add_adsic)
    end subroutine fxc_add_adsic

  end subroutine casida_work

  ! ---------------------------------------------------------
  FLOAT function casida_matrix_factor(cas, sys)
    type(casida_t),      intent(in)    :: cas
    type(electrons_t),   intent(in)    :: sys
    
    PUSH_SUB(casida_matrix_factor)
    
    casida_matrix_factor = M_ONE
    
    if(cas%type == CASIDA_VARIATIONAL) then
      casida_matrix_factor = M_TWO * casida_matrix_factor
    end if
    
    if(sys%st%d%ispin == UNPOLARIZED) then
      casida_matrix_factor = M_TWO * casida_matrix_factor
    end if
    
    POP_SUB(casida_matrix_factor)
    
  end function casida_matrix_factor

  ! ---------------------------------------------------------
  subroutine qcasida_write(cas, namespace)
    type(casida_t),    intent(in) :: cas
    type(namespace_t), intent(in) :: namespace

    integer :: iunit, ia

    if(.not.mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(qcasida_write)

    call io_mkdir(CASIDA_DIR, namespace)
    iunit = io_open(CASIDA_DIR//'q'//trim(theory_name(cas)), namespace, action='write')
    write(iunit, '(a1,a14,1x,a24,1x,a24,1x,a10,3es15.8,a2)') '#','E' , '|<f|exp(iq.r)|i>|^2', &
                                                             '<|<f|exp(iq.r)|i>|^2>','; q = (',cas%qvector(1:cas%space_dim),')'
    write(iunit, '(a1,a14,1x,a24,1x,a24,1x,10x,a15)')        '#', trim(units_abbrev(units_out%energy)), &
                                                                  trim('-'), &
                                                                  trim('-'), &
                                                                  trim('a.u.')

    if(cas%avg_order == 0) then
      do ia = 1, cas%n_pairs
        write(iunit, '(es15.8,es15.8)') units_from_atomic(units_out%energy, cas%w(cas%ind(ia))), cas%qf(cas%ind(ia))
      end do
    else
      do ia = 1, cas%n_pairs
        write(iunit, '(3es15.8)') units_from_atomic(units_out%energy, cas%w(cas%ind(ia))), &
                                  cas%qf    (cas%ind(ia)), &
                                  cas%qf_avg(cas%ind(ia))
      end do
    end if

    call io_close(iunit)

    POP_SUB(qcasida_write)

  end subroutine qcasida_write

  ! ---------------------------------------------------------
  character(len=80) pure function theory_name(cas)
    type(casida_t), intent(in) :: cas

    select case(cas%type)
      case(CASIDA_EPS_DIFF)
        theory_name = "eps_diff"
      case(CASIDA_PETERSILKA)
        theory_name = "petersilka"
      case(CASIDA_TAMM_DANCOFF)
        theory_name = "tamm_dancoff"
      case(CASIDA_VARIATIONAL)
        theory_name = "variational"
      case(CASIDA_CASIDA)
        theory_name = "casida"
      case default
        theory_name = "unknown"
    end select

  end function theory_name

  logical function isnt_degenerate(cas, st, ia, jb)
    type(casida_t),      intent(in) :: cas
    type(states_elec_t), intent(in) :: st
    integer,             intent(in) :: ia
    integer,             intent(in) :: jb

    PUSH_SUB(isnt_degenerate)

    isnt_degenerate = (abs((st%eigenval(cas%pair(ia)%a, cas%pair(ia)%kk) - st%eigenval(cas%pair(ia)%i, cas%pair(ia)%kk)) &
      - (st%eigenval(cas%pair(jb)%a, cas%pair(jb)%kk) - st%eigenval(cas%pair(jb)%i, cas%pair(jb)%kk))) > CNST(1e-8))

    POP_SUB(isnt_degenerate)
  end function isnt_degenerate

  integer function get_global_row(cas, jb_local) result(jb)
    implicit none
    type(casida_t), intent(inout) :: cas
    integer, intent(in) :: jb_local

    if(.not. cas%distributed_matrix) then
      jb = jb_local
    else
#ifdef HAVE_SCALAPACK
      jb = indxl2g(jb_local, cas%block_size, cas%proc_grid%myrow, 0, cas%proc_grid%nprow)
#endif
    end if
  end function get_global_row

  integer function get_global_col(cas, ia_local) result(ia)
    implicit none
    type(casida_t), intent(inout) :: cas
    integer, intent(in) :: ia_local

    if(.not. cas%distributed_matrix) then
      ia = ia_local
    else
#ifdef HAVE_SCALAPACK
      ia = indxl2g(ia_local, cas%block_size, cas%proc_grid%mycol, 0, cas%proc_grid%npcol)
#endif
    end if
  end function get_global_col

  subroutine local_indices(cas, ia, jb, on_this_processor, ia_local, jb_local)
    implicit none
    type(casida_t), intent(in) :: cas
    integer, intent(in) :: ia, jb
    logical, intent(out) :: on_this_processor
    integer, intent(out) :: ia_local, jb_local
#ifdef HAVE_SCALAPACK
    integer :: ia_proc, jb_proc
#endif

    if(.not. cas%distributed_matrix) then
      on_this_processor = .true.
      ia_local = ia
      jb_local = jb
    else
#ifdef HAVE_SCALAPACK
      ia_proc = indxg2p(ia, cas%block_size, cas%proc_grid%mycol, 0, cas%proc_grid%npcol)
      jb_proc = indxg2p(jb, cas%block_size, cas%proc_grid%myrow, 0, cas%proc_grid%nprow)
      if(cas%proc_grid%mycol == ia_proc .and. cas%proc_grid%myrow == jb_proc) then
        on_this_processor = .true.
        ia_local = indxg2l(ia, cas%block_size, cas%proc_grid%mycol, 0, cas%proc_grid%npcol)
        jb_local = indxg2l(jb, cas%block_size, cas%proc_grid%myrow, 0, cas%proc_grid%nprow)
      else
        on_this_processor = .false.
        ia_local = -1
        jb_local = -1
      end if
#endif
    end if
  end subroutine local_indices

#include "undef.F90"
#include "real.F90"
#include "casida_inc.F90"
#include "undef.F90"
#include "complex.F90"
#include "casida_inc.F90"

end module casida_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
