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
    integer           :: el_per_state
    character(len=80) :: wfn_list
    character(len=80) :: trandens
    logical           :: triplet        !< use triplet kernel?
    logical           :: calc_forces    !< calculate excited-state forces
    logical           :: calc_forces_kernel    !< calculate excited-state forces with kernel
    logical           :: calc_forces_scf       !< calculate excited-state forces with SCF forces
    logical           :: herm_conj      !< use Hermitian conjugate of matrix
    character(len=80) :: restart_dir
    
    integer           :: n_pairs        !< number of pairs to take into account
    type(states_pair_t), pointer :: pair(:)
    integer, pointer  :: index(:,:,:)   !< index(pair(j)%i, pair(j)%a, pair(j)%sigma) = j
    integer, pointer  :: ind(:)         !< ordering in energy of solutions

    integer :: qi_old, qa_old, mu_old   !< previous mtxel calculated in K_term
    FLOAT,   pointer  :: dpot(:)         !< previous exchange potential calculated in K_term
    CMPLX,   pointer  :: zpot(:)         !< previous exchange potential calculated in K_term

    FLOAT,   pointer  :: dmat(:,:)      !< general-purpose matrix
    FLOAT,   pointer  :: dmat_save(:,:) !< to save mat when it gets turned into the eigenvectors
    CMPLX,   pointer  :: zmat(:,:)      !< general-purpose matrix
    CMPLX,   pointer  :: zmat_save(:,:) !< to save mat when it gets turned into the eigenvectors
    FLOAT,   pointer  :: w(:)           !< The excitation energies.
    FLOAT,   pointer  :: dtm(:, :)      !< The transition matrix elements (between the many-particle states)
    CMPLX,   pointer  :: ztm(:, :)      !< The transition matrix elements (between the many-particle states)
    FLOAT,   pointer  :: f(:)           !< The (dipole) strengths
    FLOAT,   pointer  :: s(:)           !< The diagonal part of the S-matrix

    FLOAT, pointer :: rho(:, :)
    FLOAT, pointer :: fxc(:,:,:)
    FLOAT, pointer :: xc(:,:,:)

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
    integer :: idir, ik, n_filled, n_partially_filled, n_half_filled, theorylevel, iatom
    character(len=80) :: nst_string, default
    character(len=100) :: restart_filename

    PUSH_SUB(casida_run)

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

    call restart_look_and_read(sys%st, sys%gr)

    cas%el_per_state = sys%st%smear%el_per_state
    cas%nst = sys%st%nst
    cas%nik = sys%st%d%nik

    SAFE_ALLOCATE(  cas%n_occ(1:cas%nik))
    SAFE_ALLOCATE(cas%n_unocc(1:cas%nik))

    cas%n_occ(:) = 0
    do ik = 1, cas%nik
      call occupied_states(sys%st, ik, n_filled, n_partially_filled, n_half_filled)
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        call messages_not_implemented("Casida with partial occupations")
        ! Formulas are in Casida 1995 reference. The occupations are not used at all here currently.
      endif
      cas%n_occ(ik) = n_filled + n_partially_filled + n_half_filled
      cas%n_unocc(ik) = cas%nst - n_filled
      ! when we implement occupations, partially occupied levels need to be counted as both occ and unocc.
    end do

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
    !% Only <tt>eps_diff</tt> is implemented for complex wavefunctions. Note the restart data saved by each theory level
    !% is compatible with all the others.
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
      if(iand(theorylevel, CASIDA_TAMM_DANCOFF) /= 0 .or. iand(theorylevel, CASIDA_PETERSILKA) /= 0) then
        call messages_not_implemented("Tamm-Dancoff and Petersilka theory levels with complex wavefunctions")
      endif
      if((iand(theorylevel, CASIDA_VARIATIONAL) /= 0 &
        .or. iand(theorylevel, CASIDA_CASIDA) /= 0)) then
        message(1) = "Variational and full Casida theory levels do not apply to complex wavefunctions."
        call messages_fatal(1, only_root_writes = .true.)
        ! see section II.D of CV(2) paper regarding this assumption. Would be Eq. 30 with complex wfns.
      endif
    end if

    !%Variable CasidaKohnShamStates
    !%Type string
    !%Section Linear Response::Casida
    !%Default all states
    !%Description
    !% The calculation of the excitation spectrum of a system in the Casida frequency-domain
    !% formulation of linear-response time-dependent density functional theory (TDDFT)
    !% implies the use of a basis set of occupied/unoccupied Kohn-Sham orbitals. This
    !% basis set should, in principle, include all pairs formed by all occupied states,
    !% and an infinite number of unoccupied states. In practice, one has to truncate this
    !% basis set, selecting a number of occupied and unoccupied states that will form the
    !% pairs. These states are specified with this variable. If there are, say, 15 occupied
    !% states, and one sets this variable to the value "10-18", this means that occupied
    !% states from 10 to 15, and unoccupied states from 16 to 18 will be considered.
    !%
    !% This variable is a string in list form, <i>i.e.</i> expressions such as "1,2-5,8-15" are
    !% valid. You should include a non-zero number of unoccupied states and a non-zero number
    !% of occupied states.
    !%End

    write(nst_string,'(i6)') cas%nst
    write(default,'(a,a)') "1-", trim(adjustl(nst_string))
    call parse_string(datasets_check('CasidaKohnShamStates'), default, cas%wfn_list)
    write(message(1),'(a,a)') "Info: States that form the basis: ", trim(cas%wfn_list)
    Call messages_info(1)

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
      do idir = 1, MAX_DIM
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

    cas%restart_dir = trim(tmpdir)//'casida'
    cas%fromScratch = fromScratch
    call io_mkdir(trim(cas%restart_dir))

    if(cas%fromScratch) then ! remove old restart files
      if(cas%triplet) then
        call loct_rm(trim(cas%restart_dir)//'/kernel_triplet')
      else
        call loct_rm(trim(cas%restart_dir)//'/kernel')
      endif

      if(cas%calc_forces) then
        do iatom = 1, sys%geo%natoms
          do idir = 1, sys%gr%sb%dim
            write(restart_filename,'(a,a,i6.6,a,i1)') trim(cas%restart_dir), '/lr_kernel_', iatom, '_', idir
            if(cas%triplet) restart_filename = trim(restart_filename)//'_triplet'
            call loct_rm(trim(restart_filename))
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
      call casida_write(cas, sys)
    endif

    if (sys%st%d%ispin /= SPINORS) then

      if(iand(theorylevel, CASIDA_TAMM_DANCOFF) /= 0) then
        call messages_experimental("Tamm-Dancoff calculation")
        message(1) = "Info: Calculating matrix elements in the Tamm-Dancoff approximation"
        call messages_info(1)
        cas%type = CASIDA_TAMM_DANCOFF
        call casida_work(sys, hm, cas)
        call casida_write(cas, sys)
      endif

      if(iand(theorylevel, CASIDA_VARIATIONAL) /= 0) then
        call messages_experimental("CV(2)-DFT calculation")
        message(1) = "Info: Calculating matrix elements with the CV(2)-DFT theory"
        call messages_info(1)
        cas%type = CASIDA_VARIATIONAL
        call casida_work(sys, hm, cas)
        call casida_write(cas, sys)
      endif

      if(iand(theorylevel, CASIDA_CASIDA) /= 0) then
        message(1) = "Info: Calculating matrix elements with the full Casida method"
        call messages_info(1)
        cas%type = CASIDA_CASIDA
        call casida_work(sys, hm, cas)
        call casida_write(cas, sys)
      endif

      ! Doing this first, if doing the others later, takes longer, because we would use
      ! each Poisson solution for only one matrix element instead of a whole column.
      if(iand(theorylevel, CASIDA_PETERSILKA) /= 0) then
        message(1) = "Info: Calculating resonance energies via the Petersilka approximation"
        call messages_info(1)
        cas%type = CASIDA_PETERSILKA
        call casida_work(sys, hm, cas)
        call casida_write(cas, sys)
      endif

    end if

    call casida_type_end(cas)

    POP_SUB(casida_run)
  end subroutine casida_run

  ! ---------------------------------------------------------
  !> allocates stuff, and constructs the arrays pair_i and pair_j
  subroutine casida_type_init(cas, sys)
    type(casida_t),    intent(inout) :: cas
    type(system_t),    intent(in)    :: sys

    integer :: ist, ast, jpair, ik

    PUSH_SUB(casida_type_init)

    cas%states_are_real = states_are_real(sys%st)

    ! count pairs
    cas%n_pairs = 0
    do ik = 1, cas%nik
      do ast = cas%n_occ(ik) + 1, cas%nst
        if(loct_isinstringlist(ast, cas%wfn_list)) then
          do ist = 1, cas%n_occ(ik)
            if(loct_isinstringlist(ist, cas%wfn_list)) then
              cas%n_pairs = cas%n_pairs + 1
            end if
          end do
        end if
      end do
    end do

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
      SAFE_ALLOCATE(  cas%dtm(1:cas%n_pairs, 1:sys%gr%sb%dim))
    else
      SAFE_ALLOCATE( cas%zmat(1:cas%n_pairs, 1:cas%n_pairs))
      SAFE_ALLOCATE(  cas%ztm(1:cas%n_pairs, 1:sys%gr%sb%dim))
    endif
    SAFE_ALLOCATE(   cas%f(1:cas%n_pairs))
    SAFE_ALLOCATE(   cas%s(1:cas%n_pairs))
    SAFE_ALLOCATE(   cas%w(1:cas%n_pairs))
    SAFE_ALLOCATE(cas%index(1:maxval(cas%n_occ), cas%nst - maxval(cas%n_unocc) + 1:cas%nst, cas%nik))
    SAFE_ALLOCATE( cas%ind(1:cas%n_pairs))

    if(cas%calc_forces) then
      if(cas%states_are_real) then
        SAFE_ALLOCATE(cas%dmat_save(cas%n_pairs, cas%n_pairs))
      else
        SAFE_ALLOCATE(cas%zmat_save(cas%n_pairs, cas%n_pairs))
      endif
      SAFE_ALLOCATE(cas%forces(1:sys%geo%natoms, 1:sys%gr%sb%dim, 1:cas%n_pairs))
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
        if(loct_isinstringlist(ast, cas%wfn_list)) then
          do ist = 1, cas%n_occ(ik)
            if(loct_isinstringlist(ist, cas%wfn_list)) then
              cas%index(ist, ast, ik) = jpair
              cas%pair(jpair)%i = ist
              cas%pair(jpair)%a = ast
              cas%pair(jpair)%sigma = ik
              jpair = jpair + 1
            end if
          end do
        end if
      end do
    end do

    ! now let us take care of initializing the parallel stuff
    cas%parallel_in_eh_pairs = multicomm_strategy_is_parallel(sys%mc, P_STRATEGY_OTHER)
    if(cas%parallel_in_eh_pairs) then
      call mpi_grp_init(cas%mpi_grp, sys%mc%group_comm(P_STRATEGY_OTHER))
    else
      call mpi_grp_init(cas%mpi_grp, -1)
    end if

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

    FLOAT, allocatable :: rho(:, :), rho_spin(:, :)
    FLOAT, target, allocatable :: fxc_spin(:,:,:)
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

    if (cas%type /= CASIDA_EPS_DIFF) then
      ! This is to be allocated here, and is used inside K_term.
      if(.not. cas%triplet) then
        if(cas%states_are_real) then
          SAFE_ALLOCATE(cas%dpot(1:mesh%np))
        else
          SAFE_ALLOCATE(cas%zpot(1:mesh%np))
        endif
      endif
      cas%qi_old = -1
      cas%qa_old = -1
      cas%mu_old = -1
    endif
      
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

    restart_filename = trim(cas%restart_dir)//'/kernel'
    if(cas%triplet) restart_filename = trim(restart_filename)//'_triplet'

    select case(cas%type)
    case(CASIDA_EPS_DIFF)
      call solve_eps_diff()
    case(CASIDA_TAMM_DANCOFF,CASIDA_VARIATIONAL,CASIDA_CASIDA,CASIDA_PETERSILKA)
      if(cas%states_are_real) then
        call dcasida_get_matrix(cas, hm, st, mesh, cas%dmat, cas%fxc, restart_filename)
        cas%dmat = cas%dmat * casida_matrix_factor(cas, sys)
      else
        call zcasida_get_matrix(cas, hm, st, mesh, cas%zmat, cas%fxc, restart_filename)
        cas%zmat = cas%zmat * casida_matrix_factor(cas, sys)
      endif
      call solve_casida()
    end select

    if(cas%calc_forces) call casida_forces_init(cas, sys, mesh, st, hm)

    ! clean up
    if (cas%type /= CASIDA_EPS_DIFF) then
      SAFE_DEALLOCATE_P(cas%fxc)
      if(.not. cas%triplet) then
        if(cas%states_are_real) then
          SAFE_DEALLOCATE_P(cas%dpot)
        else
          SAFE_DEALLOCATE_P(cas%zpot)
        endif
      endif
    end if
    if(cas%type /= CASIDA_EPS_DIFF .or. cas%calc_forces) then
      SAFE_DEALLOCATE_P(cas%rho)
    endif

    POP_SUB(casida_work)

  contains

    ! ---------------------------------------------------------
    subroutine solve_eps_diff

      integer :: ia, actual
      FLOAT, allocatable :: w(:)

      PUSH_SUB(casida_work.solve_eps_diff)

      ! initialize progress bar
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, cas%n_pairs)

      do ia = 1, cas%n_pairs
        cas%w(ia) = st%eigenval(cas%pair(ia)%a, cas%pair(ia)%sigma) - &
                    st%eigenval(cas%pair(ia)%i, cas%pair(ia)%sigma)
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

      if(cas%states_are_real) then
        call doscillator_strengths(cas, mesh, st)
      else
        call zoscillator_strengths(cas, mesh, st)
      endif

      if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

      POP_SUB(casida_work.solve_eps_diff)
    end subroutine solve_eps_diff

    ! ---------------------------------------------------------
    subroutine solve_casida()

      FLOAT :: temp
      integer :: ia, jb
      type(states_pair_t), pointer :: p, q

      PUSH_SUB(casida_work.solve_casida)

      ! all processors with the exception of the first are done
      if (mpi_grp_is_root(cas%mpi_grp)) then

        ! complete the matrix
        do ia = 1, cas%n_pairs
          p => cas%pair(ia)
          temp = st%eigenval(p%a, p%sigma) - st%eigenval(p%i, p%sigma)
            
          do jb = ia, cas%n_pairs
            q => cas%pair(jb)
              
            ! FIXME: need the equivalent of this stuff for forces too.
            if(cas%type == CASIDA_CASIDA) then
              cas%dmat(ia, jb) = M_TWO * sqrt(temp) * cas%dmat(ia, jb) * &
                sqrt(st%eigenval(q%a, q%sigma) - st%eigenval(q%i, q%sigma))
            endif
            if(jb /= ia) cas%dmat(jb, ia) = cas%dmat(ia, jb) ! the matrix is symmetric (FIXME: actually Hermitian)
          end do
          if(cas%type == CASIDA_CASIDA) then
            cas%dmat(ia, ia) = temp**2 + cas%dmat(ia, ia)
          else
            cas%dmat(ia, ia) = cas%dmat(ia, ia) + temp
          endif
        end do

        message(1) = "Info: Diagonalizing matrix for resonance energies."
        call messages_info(1)
        ! now we diagonalize the matrix
        ! for huge matrices, perhaps we should consider ScaLAPACK here...
        call profiling_in(prof, "CASIDA_DIAGONALIZATION")
        if(cas%calc_forces) cas%dmat_save = cas%dmat ! save before gets turned into eigenvectors
        call lalg_eigensolve(cas%n_pairs, cas%dmat, cas%w)
        call profiling_out(prof)

        do ia = 1, cas%n_pairs

          if(cas%type == CASIDA_CASIDA) then
            if(cas%w(ia) < -M_EPSILON) then       
              write(message(1),'(a,i4,a)') 'Casida excitation energy', ia, ' is imaginary.'
              call messages_warning(1)
              cas%w(ia) = -sqrt(-cas%w(ia))
            else
              cas%w(ia) = sqrt(cas%w(ia))
            endif
          else
            if(cas%w(ia) < -M_EPSILON) then
              write(message(1),'(a,i4,a)') 'For whatever reason, excitation energy', ia, ' is negative.'
              write(message(2),'(a)')      'This should not happen.'
              call messages_warning(2)
           endif
          endif

          cas%ind(ia) = ia ! diagonalization returns eigenvalues in order.
        end do

        ! And let us now get the S matrix...
        if(cas%type == CASIDA_CASIDA) then
          do ia = 1, cas%n_pairs
            cas%s(ia) = (M_ONE/cas%el_per_state) / ( st%eigenval(cas%pair(ia)%a, cas%pair(ia)%sigma) - &
                                                     st%eigenval(cas%pair(ia)%i, cas%pair(ia)%sigma) )
          end do
        endif

        if(cas%states_are_real) then
          call doscillator_strengths(cas, mesh, st)
        else
          call zoscillator_strengths(cas, mesh, st)
        endif

      end if

      POP_SUB(casida_work.solve_casida)
    end subroutine solve_casida

  end subroutine casida_work

    ! ---------------------------------------------------------
    subroutine casida_forces_init(cas, sys, mesh, st, hm)
      type(casida_t), intent(inout) :: cas
      type(system_t), intent(inout) :: sys
      type(mesh_t), intent(in) :: mesh
      type(states_t), intent(inout) :: st
      type(hamiltonian_t), intent(inout) :: hm
      
      integer :: ip, iatom, idir, is1, is2, ierr, ik, ia, ist, jst, iunit
      FLOAT, allocatable :: hvar(:,:,:), dlr_hmat1(:,:,:), dl_rho(:,:), kxc(:,:,:,:)
      FLOAT, target, allocatable :: lr_fxc(:,:,:)
      CMPLX, allocatable :: zlr_hmat1(:,:,:)
      type(pert_t) :: ionic_pert
      FLOAT :: factor = CNST(1e6) ! FIXME: allow user to set
      character(len=100) :: restart_filename
      
      PUSH_SUB(casida_forces_init)
      
      if(cas%type == CASIDA_CASIDA) then
        message(1) = "Forces for Casida theory level not implemented"
        call messages_warning(1)
        POP_SUB(casida_forces_init)
        return
      endif

      if(cas%type /= CASIDA_EPS_DIFF) then
        SAFE_ALLOCATE(kxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin, 1:st%d%nspin))
        kxc = M_ZERO
        ! not spin polarized so far
        call xc_get_kxc(sys%ks%xc, mesh, cas%rho, st%d%ispin, kxc(:, :, :, :))
      endif

      message(1) = "Reading vib_modes density for calculating excited-state forces."
      call messages_info(1)
      
      SAFE_ALLOCATE(dl_rho(1:mesh%np, 1:st%d%nspin))
      if (cas%type /= CASIDA_EPS_DIFF) then
        SAFE_ALLOCATE(lr_fxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
      endif

      if(cas%type == CASIDA_EPS_DIFF) then
        cas%dmat_save = M_ZERO
        do ia = 1, cas%n_pairs
          cas%dmat_save(ia, ia) = cas%w(ia)
        enddo
      endif

      SAFE_ALLOCATE(hvar(1:mesh%np, 1:st%d%nspin, 1:1))
      if(cas%states_are_real) then
        SAFE_ALLOCATE(dlr_hmat1(cas%nst, cas%nst, cas%nik))
        SAFE_ALLOCATE(cas%dlr_hmat2(cas%n_pairs, cas%n_pairs))
        SAFE_ALLOCATE(cas%dmat2(cas%n_pairs, cas%n_pairs))
        SAFE_ALLOCATE(cas%dw2(cas%n_pairs))
      else
        SAFE_ALLOCATE(zlr_hmat1(cas%nst, cas%nst, cas%nik))
        SAFE_ALLOCATE(cas%zlr_hmat2(cas%n_pairs, cas%n_pairs))
        SAFE_ALLOCATE(cas%zmat2(cas%n_pairs, cas%n_pairs))
        SAFE_ALLOCATE(cas%zw2(cas%n_pairs))
      endif
      call pert_init(ionic_pert, PERTURBATION_IONIC, sys%gr, sys%geo)

      do iatom = 1, sys%geo%natoms
        call pert_setup_atom(ionic_pert, iatom)

        do idir = 1, mesh%sb%dim

          call pert_setup_dir(ionic_pert, idir)

          call drestart_read_lr_rho(dl_rho, sys%gr, st%d%nspin, &
            VIB_MODES_DIR, phn_rho_tag(iatom, idir), ierr)
          
          if(ierr /= 0) then
            message(1) = "Could not load vib_modes density; previous vib_modes calculation required."
            call messages_fatal(1)
          end if

          call dcalc_hvar(.true., sys, dl_rho, 1, hvar, fxc = cas%fxc)

          if(cas%states_are_real) then
            dlr_hmat1 = M_ZERO
            cas%dlr_hmat2 = M_ZERO
          else
            zlr_hmat1 = M_ZERO
            cas%zlr_hmat2 = M_ZERO
          endif

          do ik = 1, cas%nik
            if(cas%states_are_real) then
              ! occ-occ matrix elements
              call dcasida_lr_hmat1(sys, hm, ionic_pert, hvar, dlr_hmat1, 1, cas%n_occ(ik), ik)
              ! unocc-unocc matrix elements
              call dcasida_lr_hmat1(sys, hm, ionic_pert, hvar, dlr_hmat1, cas%n_occ(ik) + 1, cas%nst, ik)
            else
              ! occ-occ matrix elements
              call zcasida_lr_hmat1(sys, hm, ionic_pert, hvar, zlr_hmat1, 1, cas%n_occ(ik), ik)
              ! unocc-unocc matrix elements
              call zcasida_lr_hmat1(sys, hm, ionic_pert, hvar, zlr_hmat1, cas%n_occ(ik) + 1, cas%nst, ik)
            endif
          enddo

          if(mpi_grp_is_root(mpi_world)) then
            write(restart_filename,'(a,a,i6.6,a,i1)') trim(cas%restart_dir), '/lr_hmat1_', iatom, '_', idir
            iunit = io_open(restart_filename, action = 'write')
            do ik = 1, cas%nik
              do ist = 1, cas%nst
                do jst = ist, cas%nst
                  write(iunit, *) ist, jst, ik, dlr_hmat1(ist, jst, ik)
                enddo
              enddo
            enddo
            call io_close(iunit)
          endif

          ! use them to make two-particle matrix elements (as for eigenvalues)
          do ik = 1, cas%nik
            if(cas%states_are_real) then
              call dcasida_lr_hmat2(cas, st, dlr_hmat1, ik)
            else
              call zcasida_lr_hmat2(cas, st, zlr_hmat1, ik)
            endif
          enddo

          if (cas%type /= CASIDA_EPS_DIFF .and. cas%calc_forces_kernel) then
            forall(ip = 1:mesh%np, is1 = 1:st%d%nspin, is2 = 1:st%d%nspin)
              lr_fxc(ip, is1, is2) = sum(kxc(ip, is1, is2, :) * dl_rho(ip, :))
            end forall

            write(restart_filename,'(a,a,i6.6,a,i1)') trim(cas%restart_dir), '/lr_kernel_', iatom, '_', idir
            if(cas%triplet) restart_filename = trim(restart_filename)//'_triplet'

            if(cas%states_are_real) then
              call dcasida_get_matrix(cas, hm, st, mesh, cas%dmat2, lr_fxc, restart_filename, is_forces = .true.)
              cas%dmat2 = cas%dmat2 * casida_matrix_factor(cas, sys)
            else
              call messages_not_implemented("lr_kernel for complex wfns") ! FIXME
            endif
          else
            if(cas%states_are_real) then
              cas%dmat2 = M_ZERO
            else
              cas%zmat2 = M_ZERO
            endif
          endif

          if(cas%states_are_real) then
            cas%dmat2 = cas%dmat_save * factor + cas%dlr_hmat2 + cas%dmat2
            call lalg_eigensolve(cas%n_pairs, cas%dmat2, cas%dw2)
            do ia = 1, cas%n_pairs
              cas%forces(iatom, idir, cas%ind(ia)) = factor * cas%w(cas%ind(ia)) - cas%dw2(ia)
            enddo
          else
            cas%zmat2 = cas%zmat_save * factor + cas%zlr_hmat2 + cas%zmat2
            call lalg_eigensolve(cas%n_pairs, cas%zmat2, cas%zw2)
            do ia = 1, cas%n_pairs
              cas%forces(iatom, idir, cas%ind(ia)) = factor * cas%w(cas%ind(ia)) - real(cas%zw2(ia), REAL_PRECISION)
            enddo
          endif
        enddo
      enddo

      call pert_end(ionic_pert)
      if(cas%type /= CASIDA_EPS_DIFF) then
        SAFE_DEALLOCATE_A(kxc)
        SAFE_DEALLOCATE_A(lr_fxc)
      endif
      SAFE_DEALLOCATE_A(dl_rho)
      SAFE_DEALLOCATE_A(hvar)
      if(cas%states_are_real) then
        SAFE_DEALLOCATE_A(dlr_hmat1)
        SAFE_DEALLOCATE_P(cas%dmat2)
        SAFE_DEALLOCATE_P(cas%dw2)
      else
        SAFE_DEALLOCATE_A(zlr_hmat1)
        SAFE_DEALLOCATE_P(cas%zmat2)
        SAFE_DEALLOCATE_P(cas%zw2)
      endif
      
      if(cas%calc_forces_scf) then
        call forces_calculate(sys%gr, sys%geo, hm%ep, st)
        do ia = 1, cas%n_pairs
          do iatom = 1, sys%geo%natoms
            do idir = 1, sys%gr%sb%dim
              cas%forces(iatom, idir, ia) = cas%forces(iatom, idir, ia) + sys%geo%atom(iatom)%f(idir)
            enddo
          enddo
        enddo
      endif

      POP_SUB(casida_forces_init)

    end subroutine casida_forces_init

    ! ---------------------------------------------------------
    FLOAT function casida_matrix_factor(cas, sys)
      type(casida_t), intent(in)    :: cas
      type(system_t), intent(in)    :: sys

      FLOAT :: factor
      
      PUSH_SUB(casida_matrix_factors)
      
      factor = M_ONE

      if(cas%type == CASIDA_VARIATIONAL) then
        factor = M_TWO * factor
      endif
      
      if(sys%st%d%ispin == UNPOLARIZED) then
        factor = M_TWO * factor
      endif
      
      POP_SUB(casida_matrix_factor)

    end function casida_matrix_factor

  ! ---------------------------------------------------------
  subroutine qcasida_write(cas)
    type(casida_t), intent(in) :: cas

    integer :: iunit, ia, dim

    if(.not.mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(qcasida_write)

    dim = size(cas%dtm, 2)

    call io_mkdir(CASIDA_DIR)
    iunit = io_open(CASIDA_DIR//'q'//trim(theory_name(cas)), action='write')
    write(iunit, '(a1,a14,1x,a24,1x,a24,1x,a10,3es15.8,a2)') '#','E' , '|<f|exp(iq.r)|i>|^2', &
                                                             '<|<f|exp(iq.r)|i>|^2>','; q = (',cas%qvector(1:dim),')'
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
  subroutine casida_write(cas, sys)
    type(casida_t), intent(in) :: cas
    type(system_t), intent(in) :: sys

    character(len=5) :: str
    character(len=50) :: dir_name
    integer :: iunit, ia, jb, dim, idim
    FLOAT   :: temp

    if(.not.mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(casida_write)

    dim = size(cas%dtm, 2)

    ! output excitation energies and oscillator strengths
    call io_mkdir(CASIDA_DIR)
    iunit = io_open(CASIDA_DIR//trim(theory_name(cas)), action='write')

    if(cas%type == CASIDA_EPS_DIFF) then
      write(iunit, '(2a4)', advance='no') 'From', '  To'
      if(sys%st%d%ispin == SPIN_POLARIZED) then
        write(iunit, '(a5)', advance='no') 'Spin'
      endif
    else
      write(iunit, '(6x)', advance='no')
    endif

    write(iunit, '(1x,a15)', advance='no') 'E [' // trim(units_abbrev(units_out%energy)) // ']' 
    do idim = 1, dim
      write(iunit, '(1x,a15)', advance='no') '<' // index2axis(idim) // '> [' // trim(units_abbrev(units_out%length)) // ']' 
    enddo
    write(iunit, '(1x,a15)') '<f>'

    do ia = 1, cas%n_pairs
      if((cas%type == CASIDA_EPS_DIFF)) then
        write(iunit, '(2i4)', advance='no') cas%pair(cas%ind(ia))%i, cas%pair(cas%ind(ia))%a
        if(sys%st%d%ispin == SPIN_POLARIZED) then
          write(iunit, '(i5)', advance='no') cas%pair(cas%ind(ia))%sigma
        endif
      else
        write(iunit, '(i6)', advance='no') cas%ind(ia)
      end if
      write(iunit, '(99(1x,es15.8))') units_from_atomic(units_out%energy, cas%w(cas%ind(ia))), &
        (units_from_atomic(units_out%length, cas%dtm(cas%ind(ia), idim)), idim=1,dim), cas%f(cas%ind(ia))
    end do
    call io_close(iunit)

    if(cas%qcalc) call qcasida_write(cas)

    if(cas%type /= CASIDA_EPS_DIFF .or. cas%calc_forces) then
      dir_name = CASIDA_DIR//trim(theory_name(cas))//'_excitations'
      call io_mkdir(trim(dir_name))
    endif

    do ia = 1, cas%n_pairs
      
      write(str,'(i5.5)') ia

      ! output eigenvectors
      if(cas%type /= CASIDA_EPS_DIFF) then
        iunit = io_open(trim(dir_name)//'/'//trim(str), action='write')
        ! First, a little header
        write(iunit,'(a,es14.5)') '# Energy ['// trim(units_abbrev(units_out%energy)) // '] = ', &
          units_from_atomic(units_out%energy, cas%w(cas%ind(ia)))
        do idim = 1, dim
          write(iunit,'(a,es14.5)') '# <' // index2axis(idim) // '> ['//trim(units_abbrev(units_out%length))// '] = ', &
            units_from_atomic(units_out%length, cas%dtm(cas%ind(ia), idim))
        enddo

        ! this stuff should go BEFORE calculation of transition matrix elements!
        temp = M_ONE
        ! make the largest component positive, to specify the phase
        if( maxval(cas%dmat(:, cas%ind(ia))) - abs(minval(cas%dmat(:, cas%ind(ia)))) < -M_EPSILON) temp = -temp

        do jb = 1, cas%n_pairs
          write(iunit,*) cas%pair(jb)%i, cas%pair(jb)%a, cas%pair(jb)%sigma, temp * cas%dmat(jb, cas%ind(ia))
        end do

        if(cas%type == CASIDA_TAMM_DANCOFF .or. cas%type == CASIDA_VARIATIONAL .or. cas%type == CASIDA_PETERSILKA) then
          call write_implied_occupations(cas, iunit, cas%ind(ia))
        endif
        call io_close(iunit)
      endif

      if(cas%calc_forces .and. cas%type /= CASIDA_CASIDA) then
        iunit = io_open(trim(dir_name)//'/forces_'//trim(str)//'.xsf', action='write')
        call write_xsf_geometry(iunit, sys%geo, sys%gr%mesh, forces = cas%forces(:, :, cas%ind(ia)))
        call io_close(iunit)
      endif
    end do

    ! Calculate and write the transition densities
    if (cas%states_are_real) then
      call dget_transition_densities(cas, sys)
    else
      call zget_transition_densities(cas, sys)
    end if

    POP_SUB(casida_write)
  end subroutine casida_write

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

  ! ---------------------------------------------------------
  subroutine write_implied_occupations(cas, iunit, ind)
    type(casida_t), intent(in) :: cas
    integer,        intent(in) :: iunit
    integer,        intent(in) :: ind
    
    integer :: ik, ast, ist
    FLOAT :: occ

    PUSH_SUB(write_implied_occupations)

    write(iunit, '(a)')
    write(iunit, '(a)') '%Occupations'

    do ik = 1, cas%nik
      do ist = 1, cas%n_occ(ik)
        occ = M_ONE * cas%el_per_state
        do ast = cas%n_occ(ik) + 1, cas%nst
          if(cas%index(ist, ast, ik) == 0) cycle  ! we were not using this state
          occ = occ - abs(cas%dmat(cas%index(ist, ast, ik), ind))**2
        enddo
        write(iunit, '(f8.6,a)', advance='no') occ, ' | '
      enddo
      do ast = cas%n_occ(ik) + 1, cas%nst
        occ = M_ZERO
        do ist = 1, cas%n_occ(ik)
          if(cas%index(ist, ast, ik) == 0) cycle  ! we were not using this state
          occ = occ + abs(cas%dmat(cas%index(ist, ast, ik), ind))**2
        enddo
        write(iunit, '(f8.6)', advance='no') occ
        if(ast < cas%n_occ(ik) + cas%n_unocc(ik)) write(iunit, '(a)', advance='no') ' | '
      enddo
      write(iunit, '(a)')
    enddo
    
    write(iunit, '(a)') '%'

    POP_SUB(write_implied_occupations)
  end subroutine write_implied_occupations

  logical function isnt_degenerate(cas, st, ia, jb)
    type(casida_t), intent(in) :: cas
    type(states_t), intent(in) :: st
    integer,        intent(in) :: ia
    integer,        intent(in) :: jb

    type(states_pair_t), pointer :: pp, qq

    PUSH_SUB(isnt_degenerate)

    pp => cas%pair(ia)
    qq => cas%pair(jb)

    isnt_degenerate = (abs((st%eigenval(pp%a, pp%sigma) - st%eigenval(pp%i, pp%sigma)) &
      - (st%eigenval(qq%a, qq%sigma) - st%eigenval(qq%i, qq%sigma))) > CNST(1e-8))

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
