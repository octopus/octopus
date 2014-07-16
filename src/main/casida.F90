!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012-2013 D. Strubbe
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
!! $Id$

#include "global.h"

module casida_m
  use batch_m
  use batch_ops_m
  use calc_mode_m
  use comm_m
  use datasets_m
  use density_m
  use excited_states_m
  use forces_m
  use gauss_legendre_m
  use global_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kpoints_m
  use lalg_adv_m
  use loct_m
  use linear_response_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multicomm_m
  use output_m
  use parser_m
  use pert_m
  use poisson_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_dim_m
  use states_restart_m
  use sternheimer_m
  use system_m
  use unit_m
  use unit_system_m
  use utils_m
  use phonons_lr_m
  use xc_m
  
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

  type casida_t
    integer :: type !< CASIDA_EPS_DIFF | CASIDA_PETERSILKA | CASIDA_TAMM_DANCOFF |
                    !< CASIDA_VARIATIONAL | CASIDA_CASIDA

    logical           :: states_are_real
    integer, pointer  :: n_occ(:)       !< number of occupied states
    integer, pointer  :: n_unocc(:)     !< number of unoccupied states
    integer           :: nst            !< total number of states
    integer           :: nik
    integer           :: sb_dim         !< number of spatial dimensions
    integer           :: el_per_state
    character(len=80) :: trandens
    logical           :: triplet        !< use triplet kernel?
    logical           :: calc_forces    !< calculate excited-state forces
    logical           :: calc_forces_kernel    !< calculate excited-state forces with kernel
    logical           :: calc_forces_scf       !< calculate excited-state forces with SCF forces
    logical           :: herm_conj      !< use Hermitian conjugate of matrix
    type(restart_t)   :: restart_load
    type(restart_t)   :: restart_dump
    
    logical, pointer  :: is_included(:,:,:) !< (i, a, k) is in the basis?
    integer           :: n_pairs        !< number of pairs to take into account
    type(states_pair_t), pointer :: pair(:)
    integer, pointer  :: index(:,:,:)   !< index(pair(j)%i, pair(j)%a, pair(j)%kk) = j
    integer, pointer  :: ind(:)         !< ordering in energy of solutions

    FLOAT,   pointer  :: dmat(:,:)      !< general-purpose matrix
    FLOAT,   pointer  :: dmat_save(:,:) !< to save mat when it gets turned into the eigenvectors
    CMPLX,   pointer  :: zmat(:,:)      !< general-purpose matrix
    CMPLX,   pointer  :: zmat_save(:,:) !< to save mat when it gets turned into the eigenvectors
    FLOAT,   pointer  :: w(:)           !< The excitation energies.
    FLOAT,   pointer  :: dtm(:, :)      !< The transition matrix elements (between the many-particle states)
    CMPLX,   pointer  :: ztm(:, :)      !< The transition matrix elements (between the many-particle states)
    FLOAT,   pointer  :: f(:)           !< The (dipole) strengths
    FLOAT,   pointer  :: s(:)           !< The diagonal part of the S-matrix

    FLOAT,   pointer  :: rho(:,:)       !< density
    FLOAT,   pointer  :: fxc(:,:,:)     !< derivative of xc potential
    FLOAT             :: kernel_lrc_alpha

    FLOAT,   pointer  :: dmat2(:,:)     !< matrix to diagonalize for forces
    CMPLX,   pointer  :: zmat2(:,:)     !< matrix to diagonalize for forces
    FLOAT,   pointer  :: dlr_hmat2(:,:) !< derivative of single-particle contribution to mat
    CMPLX,   pointer  :: zlr_hmat2(:,:) !< derivative of single-particle contribution to mat
    FLOAT,   pointer  :: forces(:,:,:)  !< excited-state forces
    FLOAT,   pointer  :: dw2(:)         !< perturbed excitation energies.
    FLOAT,   pointer  :: zw2(:)         !< perturbed excitation energies.

    ! variables for momentum-transfer-dependent calculation
    logical           :: qcalc
    FLOAT             :: qvector(MAX_DIM)
    FLOAT,   pointer  :: qf(:)
    FLOAT,   pointer  :: qf_avg(:)      !< Directionally averaged intensity
    integer           :: avg_order      !< Quadrature order for directional averaging (Gauss-Legendre scheme) 

    logical           :: parallel_in_eh_pairs
    type(mpi_grp_t)   :: mpi_grp
    logical           :: fromScratch
  end type casida_t

  type casida_save_pot_t
    integer :: qi                    !< previous mtxel calculated in K_term
    integer :: qa                    !< previous mtxel calculated in K_term
    integer :: mu                    !< previous mtxel calculated in K_term
    FLOAT,   pointer  :: dpot(:)     !< previous exchange potential calculated in K_term
    CMPLX,   pointer  :: zpot(:)     !< previous exchange potential calculated in K_term    
  end type casida_save_pot_t

  type(profile_t), save :: prof

contains

  subroutine casida_run_init()
    
    PUSH_SUB(casida_run_init)
    
    ! Pure 'other' parallelization is a bad idea. Trying to solve the Poisson equation separately on each node
    ! consumes excessive memory and time (easily more than is available). In principle, the line below would setup
    ! joint domain/other parallelization, but 'other' parallelization takes precedence, especially since
    ! multicomm_init does not know the actual problem size and uses a fictitious value of 10000, making it
    ! impossible to choose joint parallelization wisely, and generally resulting in a choice of only one domain
    ! group. FIXME! --DAS

    ! call calc_mode_set_parallelization(P_STRATEGY_OTHER, default = .true).
    call calc_mode_set_parallelization(P_STRATEGY_OTHER, default = .false.) ! enabled, but not default

    POP_SUB(casida_run_init)
  end subroutine casida_run_init

  ! ---------------------------------------------------------
  subroutine casida_run(sys, hm, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    logical,             intent(inout) :: fromScratch

    type(casida_t) :: cas
    type(block_t) :: blk
    integer :: idir, theorylevel, iatom
    character(len=100) :: restart_filename
    type(profile_t), save :: prof
    logical :: is_frac_occ
    type(restart_t) :: gs_restart

    PUSH_SUB(casida_run)
    call profiling_in(prof, 'CASIDA')

    if (simul_box_is_periodic(sys%gr%sb)) then
      message(1) = "Casida oscillator strengths will be incorrect in periodic systems."
      call messages_warning(1)
    end if

    if(kpoints_number(sys%gr%sb%kpoints) > 1) then
      ! Hartree matrix elements may not be correct, not tested anyway. --DAS
      call messages_not_implemented("Casida with k-points")
    endif

    message(1) = 'Info: Starting Casida linear-response calculation.'
    call messages_info(1)

    call restart_init(gs_restart, RESTART_GS, RESTART_TYPE_LOAD, sys%st%dom_st_kpt_mpi_grp, &
                      mesh=sys%gr%mesh, sb=sys%gr%sb, exact=.true.)
    call states_look_and_load(gs_restart, sys%st, sys%gr)
    call restart_end(gs_restart)

    cas%el_per_state = sys%st%smear%el_per_state
    cas%nst = sys%st%nst
    cas%nik = sys%st%d%nik
    cas%sb_dim = sys%gr%sb%dim
    SAFE_ALLOCATE(cas%n_occ(1:sys%st%d%nik))
    SAFE_ALLOCATE(cas%n_unocc(1:sys%st%d%nik))

    call states_count_pairs(sys%st, cas%n_pairs, cas%n_occ, cas%n_unocc, cas%is_included, is_frac_occ)
    if(is_frac_occ) then
      call messages_not_implemented("Casida with partial occupations")
      ! Formulas are in Casida 1995 reference. The occupations are not used at all here currently.
    endif

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


    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call messages_info(1)
    call system_h_setup(sys, hm)

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

    call parse_integer(datasets_check('CasidaTheoryLevel'), &
      CASIDA_EPS_DIFF + CASIDA_PETERSILKA + CASIDA_CASIDA, theorylevel)

    if (states_are_complex(sys%st)) then
      if((iand(theorylevel, CASIDA_VARIATIONAL) /= 0 &
        .or. iand(theorylevel, CASIDA_CASIDA) /= 0)) then
        message(1) = "Variational and full Casida theory levels do not apply to complex wavefunctions."
        call messages_fatal(1, only_root_writes = .true.)
        ! see section II.D of CV(2) paper regarding this assumption. Would be Eq. 30 with complex wfns.
      endif
    end if

    !%Variable CasidaTransitionDensities
    !%Type string
    !%Section Linear Response::Casida
    !%Default write none
    !%Description
    !% Specifies which transition densities are to be calculated and written down. The
    !% transition density for the many-body state <i>n</i> will be written to a file called
    !% <tt>rho0n</tt> prefixed by the theory level.
    !%
    !% This variable is a string in list form, <i>i.e.</i> expressions such as "1,2-5,8-15" are
    !% valid.
    !%End
    call parse_string(datasets_check('CasidaTransitionDensities'), "0", cas%trandens)

    if(cas%trandens /= "0") call io_function_read_how(sys%gr%sb, sys%outp%how)

    !%Variable CasidaMomentumTransfer
    !%Type block
    !%Section Linear Response::Casida
    !%Default 0
    !%Description
    !% Momentum-transfer vector for the calculation of the dynamic structure
    !% factor. When this variable is set, the transition rates are determined
    !% using an exponential operator instead of the normal dipole one.
    !%End

    if(parse_block(datasets_check('CasidaMomentumTransfer'), blk)==0) then
      do idir = 1, cas%sb_dim
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
      !% averaging over the results from a set of <i>q</i>-vectors. The vectors
      !% are generated using Gauss-Legendre quadrature scheme [see <i>e.g.</i>
      !% K. Atkinson, <i>J. Austral. Math. Soc.</i> <b>23</b>, 332 (1982)], and this
      !% variable determines the order of the scheme.
      !%End
      call parse_integer(datasets_check('CasidaQuadratureOrder'), 5, cas%avg_order)
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
      call parse_logical(datasets_check('CasidaCalcTriplet'), .false., cas%triplet)
    else
      cas%triplet = .false.
    endif

    if(cas%triplet) then
      call messages_experimental("Casida triplet calculation")
      message(1) = "Info: Using triplet kernel. Oscillator strengths will be for spin magnetic-dipole field."
      call messages_info(1)
    endif

    !%Variable CasidaHermitianConjugate
    !%Type logical
    !%Section Linear Response::Casida
    !%Default false
    !%Description
    !% The Casida matrix is Hermitian, so it should not matter whether we calculate the upper or
    !% lower diagonal. Numerical issues may cause small differences however. Use this variable to
    !% calculate the Hermitian conjugate of the usual matrix, for testing.
    !%End
    call parse_logical(datasets_check('CasidaHermitianConjugate'), .false., cas%herm_conj)

    !%Variable CasidaCalcForces
    !%Type logical
    !%Section Linear Response::Casida
    !%Default false
    !%Description
    !% (Experimental) Enable calculation of excited-state forces. Requires previous <tt>vib_modes</tt> calculation.
    !%End
    call parse_logical(datasets_check('CasidaCalcForces'), .false., cas%calc_forces)
    if(cas%calc_forces) then
      call messages_experimental("Excited-state forces calculation")

      !%Variable CasidaCalcForcesKernel
      !%Type logical
      !%Section Linear Response::Casida
      !%Default true
      !%Description
      !% If false, the derivative of the kernel will not be included in the excited-state force calculation.
      !%End
      call parse_logical(datasets_check('CasidaCalcForcesKernel'), .true., cas%calc_forces_kernel)

      !%Variable CasidaCalcForcesSCF
      !%Type logical
      !%Section Linear Response::Casida
      !%Default false
      !%Description
      !% If true, the ground-state forces will be included in the excited-state forces, so they are total forces.
      !% If false, the excited-state forces that are produced are only the gradients of the excitation energy.
      !%End
      call parse_logical(datasets_check('CasidaCalcForcesSCF'), .false., cas%calc_forces_scf)
    endif

    ! Initialize structure
    call casida_type_init(cas, sys)

    cas%fromScratch = fromScratch

    if(cas%fromScratch) then ! remove old restart files
      if(cas%triplet) then
        call restart_rm(cas%restart_dump, 'kernel_triplet')
      else
        call restart_rm(cas%restart_dump, 'kernel')
      endif

      if(cas%calc_forces) then
        do iatom = 1, sys%geo%natoms
          do idir = 1, cas%sb_dim
            write(restart_filename,'(a,i6.6,a,i1)') 'lr_kernel_', iatom, '_', idir
            if(cas%triplet) restart_filename = trim(restart_filename)//'_triplet'
            call restart_rm(cas%restart_dump, restart_filename)

            write(restart_filename,'(a,i6.6,a,i1)') 'lr_hmat1_', iatom, '_', idir
            call restart_rm(cas%restart_dump, restart_filename)
          enddo
        enddo
      endif
    endif

    ! First, print the differences between KS eigenvalues (first approximation to the excitation energies).
    if(iand(theorylevel, CASIDA_EPS_DIFF) /= 0) then
      message(1) = "Info: Approximating resonance energies through KS eigenvalue differences"
      call messages_info(1)
      cas%type = CASIDA_EPS_DIFF
      call casida_work(sys, hm, cas)
    endif

    if (sys%st%d%ispin /= SPINORS) then

      if(iand(theorylevel, CASIDA_TAMM_DANCOFF) /= 0) then
        call messages_experimental("Tamm-Dancoff calculation")
        message(1) = "Info: Calculating matrix elements in the Tamm-Dancoff approximation"
        call messages_info(1)
        cas%type = CASIDA_TAMM_DANCOFF
        call casida_work(sys, hm, cas)
      endif

      if(iand(theorylevel, CASIDA_VARIATIONAL) /= 0) then
        call messages_experimental("CV(2)-DFT calculation")
        message(1) = "Info: Calculating matrix elements with the CV(2)-DFT theory"
        call messages_info(1)
        cas%type = CASIDA_VARIATIONAL
        call casida_work(sys, hm, cas)
      endif

      if(iand(theorylevel, CASIDA_CASIDA) /= 0) then
        message(1) = "Info: Calculating matrix elements with the full Casida method"
        call messages_info(1)
        cas%type = CASIDA_CASIDA
        call casida_work(sys, hm, cas)
      endif

      ! Doing this first, if doing the others later, takes longer, because we would use
      ! each Poisson solution for only one matrix element instead of a whole column.
      if(iand(theorylevel, CASIDA_PETERSILKA) /= 0) then
        message(1) = "Info: Calculating resonance energies via the Petersilka approximation"
        call messages_info(1)
        cas%type = CASIDA_PETERSILKA
        call casida_work(sys, hm, cas)
      endif

    end if

    call casida_type_end(cas)

    call profiling_out(prof)
    POP_SUB(casida_run)
  end subroutine casida_run

  ! ---------------------------------------------------------
  !> allocates stuff, and constructs the arrays pair_i and pair_j
  subroutine casida_type_init(cas, sys)
    type(casida_t),    intent(inout) :: cas
    type(system_t),    intent(in)    :: sys

    integer :: ist, ast, jpair, ik

    PUSH_SUB(casida_type_init)

    cas%kernel_lrc_alpha = sys%ks%xc%kernel_lrc_alpha
    cas%states_are_real = states_are_real(sys%st)

    write(message(1), '(a,i9)') "Number of occupied-unoccupied pairs: ", cas%n_pairs
    call messages_info(1)

    if(cas%n_pairs < 1) then
      message(1) = "No Casida pairs -- maybe there are no unoccupied states?"
      call messages_fatal(1, only_root_writes = .true.)
    end if

    if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

    ! allocate stuff
    SAFE_ALLOCATE(cas%pair(1:cas%n_pairs))
    if(cas%states_are_real) then
      SAFE_ALLOCATE( cas%dmat(1:cas%n_pairs, 1:cas%n_pairs))
      SAFE_ALLOCATE(  cas%dtm(1:cas%n_pairs, 1:cas%sb_dim))
    else
      SAFE_ALLOCATE( cas%zmat(1:cas%n_pairs, 1:cas%n_pairs))
      SAFE_ALLOCATE(  cas%ztm(1:cas%n_pairs, 1:cas%sb_dim))
    endif
    SAFE_ALLOCATE(   cas%f(1:cas%n_pairs))
    SAFE_ALLOCATE(   cas%s(1:cas%n_pairs))
    SAFE_ALLOCATE(   cas%w(1:cas%n_pairs))
    SAFE_ALLOCATE(cas%index(1:maxval(cas%n_occ), cas%nst - maxval(cas%n_unocc) + 1:cas%nst, cas%nik))
    SAFE_ALLOCATE( cas%ind(1:cas%n_pairs))

    if(cas%calc_forces) then
      if(cas%states_are_real) then
        SAFE_ALLOCATE(cas%dmat_save(1:cas%n_pairs, 1:cas%n_pairs))
      else
        SAFE_ALLOCATE(cas%zmat_save(1:cas%n_pairs, 1:cas%n_pairs))
      endif
      SAFE_ALLOCATE(cas%forces(1:sys%geo%natoms, 1:cas%sb_dim, 1:cas%n_pairs))
    endif

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

    SAFE_DEALLOCATE_P(cas%is_included)

    ! now let us take care of initializing the parallel stuff
    cas%parallel_in_eh_pairs = multicomm_strategy_is_parallel(sys%mc, P_STRATEGY_OTHER)
    if(cas%parallel_in_eh_pairs) then
      call mpi_grp_init(cas%mpi_grp, sys%mc%group_comm(P_STRATEGY_OTHER))
    else
      call mpi_grp_init(cas%mpi_grp, -1)
    end if

    call restart_init(cas%restart_dump, RESTART_CASIDA, RESTART_TYPE_DUMP, mpi_world, mesh=sys%gr%mesh, sb=sys%gr%sb)
    call restart_init(cas%restart_load, RESTART_CASIDA, RESTART_TYPE_LOAD, mpi_world, mesh=sys%gr%mesh, sb=sys%gr%sb)

    POP_SUB(casida_type_init)
  end subroutine casida_type_init


  ! ---------------------------------------------------------
  subroutine casida_type_end(cas)
    type(casida_t), intent(inout) :: cas

    PUSH_SUB(casida_type_end)

    ASSERT(associated(cas%pair))
    SAFE_DEALLOCATE_P(cas%pair)
    SAFE_DEALLOCATE_P(cas%index)
    if(cas%states_are_real) then
      SAFE_DEALLOCATE_P(cas%dmat)
      SAFE_DEALLOCATE_P(cas%dtm)
    else
      SAFE_DEALLOCATE_P(cas%zmat)
      SAFE_DEALLOCATE_P(cas%ztm)
    endif
    SAFE_DEALLOCATE_P(cas%s)
    SAFE_DEALLOCATE_P(cas%f)
    SAFE_DEALLOCATE_P(cas%w)
    SAFE_DEALLOCATE_P(cas%ind)

    if(cas%qcalc) then
      SAFE_DEALLOCATE_P(cas%qf)
      SAFE_DEALLOCATE_P(cas%qf_avg)
    end if

    SAFE_DEALLOCATE_P(cas%n_occ)
    SAFE_DEALLOCATE_P(cas%n_unocc)

    if(cas%calc_forces) then
      if(cas%states_are_real) then
        SAFE_DEALLOCATE_P(cas%dmat_save)
      else
        SAFE_DEALLOCATE_P(cas%zmat_save)
      endif
      SAFE_DEALLOCATE_P(cas%forces)
    endif

    call restart_end(cas%restart_dump)
    call restart_end(cas%restart_load)

    POP_SUB(casida_type_end)
  end subroutine casida_type_end


  ! ---------------------------------------------------------
  !> this subroutine calculates electronic excitation energies using
  !! the matrix formulation of M. Petersilka, or of M. Casida
  subroutine casida_work(sys, hm, cas)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    type(casida_t),         intent(inout) :: cas

    type(states_t), pointer :: st
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
    endif
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

      call states_total_density(st, mesh, cas%rho)
      if(cas%triplet) then
        SAFE_ALLOCATE(rho_spin(1:mesh%np, 1:2))
        SAFE_ALLOCATE(fxc_spin(1:mesh%np, 1:2, 1:2))

        fxc_spin = M_ZERO
        rho_spin(:, 1) = M_HALF * cas%rho(:, 1)
        rho_spin(:, 2) = M_HALF * cas%rho(:, 1)

        call xc_get_fxc(sys%ks%xc, mesh, rho_spin, SPIN_POLARIZED, fxc_spin)
        cas%fxc(:, 1, 1) = M_HALF * (fxc_spin(:, 1, 1) - fxc_spin(:, 1, 2))

        SAFE_DEALLOCATE_A(rho_spin)
        SAFE_DEALLOCATE_A(fxc_spin)
      else
        call xc_get_fxc(sys%ks%xc, mesh, cas%rho, st%d%ispin, cas%fxc)
      endif
    end if

    restart_filename = 'kernel'
    if(cas%triplet) restart_filename = trim(restart_filename)//'_triplet'

    select case(cas%type)
    case(CASIDA_EPS_DIFF)
      call solve_eps_diff()
    case(CASIDA_TAMM_DANCOFF,CASIDA_VARIATIONAL,CASIDA_CASIDA,CASIDA_PETERSILKA)
      if(cas%states_are_real) then
        call dcasida_get_matrix(cas, hm, st, mesh, cas%dmat, cas%fxc, restart_filename)
        cas%dmat = cas%dmat * casida_matrix_factor(cas, sys)
        call dcasida_solve(cas, st)
      else
        call zcasida_get_matrix(cas, hm, st, mesh, cas%zmat, cas%fxc, restart_filename)
        cas%zmat = cas%zmat * casida_matrix_factor(cas, sys)
        call zcasida_solve(cas, st)
      endif
    end select

    if (mpi_grp_is_root(cas%mpi_grp)) then
      if(cas%states_are_real) then
        call doscillator_strengths(cas, mesh, st)
      else
        call zoscillator_strengths(cas, mesh, st)
      endif
    endif

    if(cas%calc_forces) then
      if(cas%states_are_real) then
        call dcasida_forces(cas, sys, mesh, st, hm)
      else
        call zcasida_forces(cas, sys, mesh, st, hm)
      endif
    endif

    if(cas%states_are_real) then
      call dcasida_write(cas, sys)
    else
      call zcasida_write(cas, sys)
    endif

    ! clean up
    if(cas%type /= CASIDA_EPS_DIFF .or. cas%calc_forces) then
      SAFE_DEALLOCATE_P(cas%fxc)
      SAFE_DEALLOCATE_P(cas%rho)
    endif

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
          message(1) = "There are negative unocc-occ KS eigenvalue differences."
          message(2) = "Probably this indicates an inconsistency in occupations between gs and unocc calculations."
          call messages_fatal(2, only_root_writes = .true.)
        endif
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ia, cas%n_pairs)
      end do

      SAFE_ALLOCATE(w(1:cas%n_pairs))
      w = cas%w
      call sort(w, cas%ind)
      SAFE_DEALLOCATE_A(w)

      if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

      POP_SUB(casida_work.solve_eps_diff)
    end subroutine solve_eps_diff

  end subroutine casida_work

  ! ---------------------------------------------------------
  FLOAT function casida_matrix_factor(cas, sys)
    type(casida_t), intent(in)    :: cas
    type(system_t), intent(in)    :: sys
    
    PUSH_SUB(casida_matrix_factor)
    
    casida_matrix_factor = M_ONE
    
    if(cas%type == CASIDA_VARIATIONAL) then
      casida_matrix_factor = M_TWO * casida_matrix_factor
    endif
    
    if(sys%st%d%ispin == UNPOLARIZED) then
      casida_matrix_factor = M_TWO * casida_matrix_factor
    endif
    
    POP_SUB(casida_matrix_factor)
    
  end function casida_matrix_factor

  ! ---------------------------------------------------------
  subroutine qcasida_write(cas)
    type(casida_t), intent(in) :: cas

    integer :: iunit, ia

    if(.not.mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(qcasida_write)

    call io_mkdir(CASIDA_DIR)
    iunit = io_open(CASIDA_DIR//'q'//trim(theory_name(cas)), action='write')
    write(iunit, '(a1,a14,1x,a24,1x,a24,1x,a10,3es15.8,a2)') '#','E' , '|<f|exp(iq.r)|i>|^2', &
                                                             '<|<f|exp(iq.r)|i>|^2>','; q = (',cas%qvector(1:cas%sb_dim),')'
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
    type(casida_t), intent(in) :: cas
    type(states_t), intent(in) :: st
    integer,        intent(in) :: ia
    integer,        intent(in) :: jb

    type(states_pair_t), pointer :: pp, qq

    PUSH_SUB(isnt_degenerate)

    pp => cas%pair(ia)
    qq => cas%pair(jb)

    isnt_degenerate = (abs((st%eigenval(pp%a, pp%kk) - st%eigenval(pp%i, pp%kk)) &
      - (st%eigenval(qq%a, qq%kk) - st%eigenval(qq%i, qq%kk))) > CNST(1e-8))

    POP_SUB(isnt_degenerate)
  end function isnt_degenerate

#include "undef.F90"
#include "real.F90"
#include "casida_inc.F90"
#include "undef.F90"
#include "complex.F90"
#include "casida_inc.F90"

end module casida_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
