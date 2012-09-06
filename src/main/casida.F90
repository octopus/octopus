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

module casida_m
  use calc_mode_m
  use comm_m
  use datasets_m
  use density_m
  use excited_states_m
  use gauss_legendre_m
  use global_m
  use output_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use lalg_adv_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multicomm_m
  use parser_m
  use poisson_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_dim_m
  use system_m
  use unit_m
  use unit_system_m
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

    integer, pointer  :: n_occ(:)       !< number of occupied states
    integer, pointer  :: n_unocc(:)     !< number of unoccupied states
    integer           :: nspin
    integer           :: el_per_state
    character(len=80) :: wfn_list
    character(len=80) :: trandens
    logical           :: triplet        !< use triplet kernel?
    character(len=80) :: restart_file

    integer           :: n_pairs        !< number of pairs to take into account
    type(states_pair_t), pointer :: pair(:)
    integer, pointer  :: index(:,:,:)   !< index(pair(j)%i, pair(j)%a, pair(j)%sigma) = j

    FLOAT,   pointer  :: mat(:,:)       !< general-purpose matrix
    FLOAT,   pointer  :: w(:)           !< The excitation energies.
    FLOAT,   pointer  :: tm(:, :)       !< The transition matrix elements (between the many-particle states)
    FLOAT,   pointer  :: f(:)           !< The (dipole) strengths
    FLOAT,   pointer  :: s(:)           !< The diagonal part of the S-matrix

    ! variables for momentum-transfer-dependent calculation
    logical           :: qcalc
    FLOAT             :: qvector(MAX_DIM)
    FLOAT,   pointer  :: qf(:)
    FLOAT,   pointer  :: qf_avg(:)      !< Directionally averaged intensity
    integer           :: avg_order      !< Quadrature order for directional averaging (Gauss-Legendre scheme) 

    logical           :: parallel_in_eh_pairs
    type(mpi_grp_t)   :: mpi_grp
  end type casida_t

contains

  subroutine casida_run_init()
    
    PUSH_SUB(casida_run_init)
    
    call calc_mode_set_parallelization(P_STRATEGY_OTHER, default = .true.)

    POP_SUB(casida_run_init)
  end subroutine casida_run_init

  ! ---------------------------------------------------------
  subroutine casida_run(sys, hm, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    logical,             intent(inout) :: fromScratch

    type(casida_t) :: cas
    type(block_t) :: blk
    integer :: idir, ik, n_filled, n_partially_filled, n_half_filled, theorylevel
    character(len=80) :: nst_string, default

    PUSH_SUB(casida_run)

    if (simul_box_is_periodic(sys%gr%sb)) then
      message(1) = "Casida oscillator strengths will be incorrect in periodic systems."
      call messages_warning(1)
    end if

    message(1) = 'Info: Starting Casida linear-response calculation.'
    call messages_info(1)

    call restart_look_and_read(sys%st, sys%gr)

    cas%el_per_state = sys%st%smear%el_per_state
    cas%nspin = sys%st%d%nspin

    SAFE_ALLOCATE(  cas%n_occ(1:cas%nspin))
    SAFE_ALLOCATE(cas%n_unocc(1:cas%nspin))

    cas%n_occ(:) = 0
    do ik = 1, cas%nspin
      call occupied_states(sys%st, ik, n_filled, n_partially_filled, n_half_filled)
      cas%n_occ(ik) = n_filled + n_partially_filled + n_half_filled
      cas%n_unocc(ik) = sys%st%nst - cas%n_occ(ik)
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
    !% Only <tt>eps_diff</tt> is available for spinors. Note the restart data saved by each theory level
    !% is compatible with all the others.
    !%Option eps_diff 1
    !% Difference of eigenvalues, <i>i.e.</i> independent-particle approximation.
    !%Option petersilka 2
    !% The Petersilka approximation uses only the diagonal part of the Tamm-Dancoff matrix.
    !% This is acceptable if there is little mixing between single-particle transitions.
    !% Ref: M Petersilka, UJ Gossmann, and EKU Gross, <i>Phys. Rev. Lett.</i> <b>76</b>, 1212 (1996).
    !%Option tamm_dancoff 4
    !% The Tamm-Dancoff approximation uses only occupied-unoccupied transitions and not
    !% unoccupied-occupied transitions.
    !% Ref: S Hirata and M Head-Gordon, <i>Chem. Phys. Lett.</i> <b>314</b>, 291 (1999).
    !%Option variational 8
    !% Second-order constrained variational theory CV(2)-DFT.
    !% Ref: T Ziegler, M Seth, M Krykunov, J Autschbach, and F Wang,
    !% <i>J. Chem. Phys.</i> <b>130</b>, 154102 (2009).
    !%Option lrtddft_casida 16
    !% The full Casida method.
    !% Ref: C Jamorski, ME Casida, and DR Salahub, <i>J. Chem. Phys.</i> <b>104</b>, 5134 (1996)
    !% and ME Casida, "Time-dependent density functional response theory for molecules,"
    !% in <i>Recent Advances in Density Functional Methods</i>, edited by DE Chong, vol. 1
    !% of <i>Recent Advances in Computational Chemistry</i>, pp. 155-192 (World Scientific,
    !% Singapore, 1995).
    !%End

    call parse_integer(datasets_check('CasidaTheoryLevel'), &
      CASIDA_EPS_DIFF + CASIDA_PETERSILKA + CASIDA_CASIDA, theorylevel)

    if (states_are_complex(sys%st) .and. &
      (iand(theorylevel, CASIDA_TAMM_DANCOFF) /= 0 .or. iand(theorylevel, CASIDA_VARIATIONAL) /= 0 &
       .or. iand(theorylevel, CASIDA_CASIDA) /= 0)) then
      message(1) = "Tamm-Dancoff, variational, and full Casida theory levels do not apply to complex wavefunctions."
      call messages_fatal(1)
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

    write(nst_string,'(i6)') sys%st%nst
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
      message(1) = "Info: Calculating IXS/EELS transition rates."
      call messages_info(1)
      cas%qcalc = .true.
    else
      cas%qvector(:) = M_ZERO
      cas%qcalc = .false.
    end if

    !%Variable CasidaQuadratureOrder
    !%Type integer
    !%Section Linear Response::Casida
    !%Default 5
    !%Description
    !% Directionally averaged dynamic structure factor is calculated by
    !% averaging over the results from a set of <i>q</i>-vectors. The vectors
    !% are generated using Gauss-Legendre quadrature scheme [see <i>e.g.</i>
    !% K. Atkinson, <i>J. Austral. Math. Soc.</i> <b>23</b>, 332 (1982)], and this
    !% variable determines the order of the scheme.
    !%End
    call parse_integer(datasets_check('CasidaQuadratureOrder'), 5, cas%avg_order)

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

    ! Initialize structure
    call casida_type_init(cas, sys%gr%sb%dim, sys%mc)

    cas%restart_file = trim(tmpdir)//'casida-restart'
    if(cas%triplet) cas%restart_file = trim(cas%restart_file)//'-triplet'

    if(fromScratch) call loct_rm(trim(cas%restart_file)) ! restart

    ! First, print the differences between KS eigenvalues (first approximation to the
    ! excitation energies, or rather, to the DOS).
    if(iand(theorylevel, CASIDA_EPS_DIFF) /= 0) then
      message(1) = "Info: Approximating resonance energies through KS eigenvalue differences"
      call messages_info(1)
      cas%type = CASIDA_EPS_DIFF
      call casida_work(sys, hm, cas)
      call casida_write(cas, sys)
    endif

    if (sys%st%d%ispin /= SPINORS) then

      if(iand(theorylevel, CASIDA_PETERSILKA) /= 0) then
        message(1) = "Info: Calculating resonance energies via the Petersilka approximation"
        call messages_info(1)
        cas%type = CASIDA_PETERSILKA
        call casida_work(sys, hm, cas)
        call casida_write(cas, sys)
      endif

      if(iand(theorylevel, CASIDA_TAMM_DANCOFF) /= 0) then
        message(1) = "Info: Calculating resonance energies in the Tamm-Dancoff approximation"
        call messages_info(1)
        cas%type = CASIDA_TAMM_DANCOFF
        call casida_work(sys, hm, cas)
        call casida_write(cas, sys)
      endif

      if(iand(theorylevel, CASIDA_VARIATIONAL) /= 0) then
        message(1) = "Info: Calculating resonance energies with the CV(2)-DFT theory"
        call messages_info(1)
        cas%type = CASIDA_VARIATIONAL
        call casida_work(sys, hm, cas)
        call casida_write(cas, sys)
      endif

      if(iand(theorylevel, CASIDA_CASIDA) /= 0) then
        message(1) = "Info: Calculating resonance energies with the full Casida method"
        call messages_info(1)
        cas%type = CASIDA_CASIDA
        call casida_work(sys, hm, cas)
        call casida_write(cas, sys)
      endif

    end if

    call casida_type_end(cas)

    POP_SUB(casida_run)
  end subroutine casida_run

  ! ---------------------------------------------------------
  !> allocates stuff, and constructs the arrays pair_i and pair_j
  subroutine casida_type_init(cas, dim, mc)
    type(casida_t),    intent(inout) :: cas
    integer,           intent(in)    :: dim
    type(multicomm_t), intent(in)    :: mc

    integer :: ist, ast, jpair, ik

    PUSH_SUB(casida_type_init)

    ! count pairs
    cas%n_pairs = 0
    do ik = 1, cas%nspin
      do ast = cas%n_occ(ik) + 1, cas%n_occ(ik) + cas%n_unocc(ik)
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
      call messages_fatal(1)
    end if

    if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

    ! allocate stuff
    SAFE_ALLOCATE(cas%pair(1:cas%n_pairs))
    SAFE_ALLOCATE( cas%mat(1:cas%n_pairs, 1:cas%n_pairs))
    SAFE_ALLOCATE(  cas%tm(1:cas%n_pairs, 1:dim))
    SAFE_ALLOCATE(   cas%f(1:cas%n_pairs))
    SAFE_ALLOCATE(   cas%s(1:cas%n_pairs))
    SAFE_ALLOCATE(   cas%w(1:cas%n_pairs))
    SAFE_ALLOCATE(cas%index(1:maxval(cas%n_occ), minval(cas%n_occ):maxval(cas%n_occ + cas%n_unocc), cas%nspin))

    if(cas%qcalc) then
      SAFE_ALLOCATE( cas%qf    (1:cas%n_pairs))
      SAFE_ALLOCATE( cas%qf_avg(1:cas%n_pairs))
    end if

    cas%index(:,:,:) = 0

    ! create pairs
    jpair = 1
    do ik = 1, cas%nspin
      do ast = cas%n_occ(ik) + 1, cas%n_occ(ik) + cas%n_unocc(ik)
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
    cas%parallel_in_eh_pairs = multicomm_strategy_is_parallel(mc, P_STRATEGY_OTHER)
    if(cas%parallel_in_eh_pairs) then
      call mpi_grp_init(cas%mpi_grp, mc%group_comm(P_STRATEGY_OTHER))
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
    SAFE_DEALLOCATE_P(cas%mat)
    SAFE_DEALLOCATE_P(cas%tm)
    SAFE_DEALLOCATE_P(cas%s)
    SAFE_DEALLOCATE_P(cas%f)
    SAFE_DEALLOCATE_P(cas%w)

    if(cas%qcalc) then
      SAFE_DEALLOCATE_P(cas%qf)
      SAFE_DEALLOCATE_P(cas%qf_avg)
    end if

    SAFE_DEALLOCATE_P(cas%n_occ)
    SAFE_DEALLOCATE_P(cas%n_unocc)

    POP_SUB(casida_type_end)
  end subroutine casida_type_end


  ! ---------------------------------------------------------
  !> this subroutine calculates electronic excitation energies using
  !! the matrix formulation of M. Petersilka, or of M. Casida
  subroutine casida_work(sys, hm, cas)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(in)    :: hm
    type(casida_t),         intent(inout) :: cas

    logical, allocatable :: saved_K(:, :) ! which matrix elements have been loaded
    type(states_t), pointer :: st
    type(mesh_t),   pointer :: mesh

    FLOAT, allocatable :: rho(:, :), rho_spin(:, :), fxc(:,:,:), pot(:)
    integer :: qi_old, qa_old, mu_old

    PUSH_SUB(casida_work)

    ! sanity checks
    ASSERT(cas%type >= CASIDA_EPS_DIFF .and. cas%type <= CASIDA_CASIDA)

    ! some shortcuts
    st => sys%st
    mesh => sys%gr%mesh

    ! initialize stuff
    SAFE_ALLOCATE(saved_K(1:cas%n_pairs, 1:cas%n_pairs))
    cas%mat = M_ZERO
    saved_K = .false.
    cas%tm  = M_ZERO
    cas%f   = M_ZERO
    cas%w   = M_ZERO
    cas%s   = M_ZERO
    if(cas%qcalc) then
      cas%qf     = M_ZERO
      cas%qf_avg = M_ZERO
    end if

    ! load saved matrix elements
    call load_saved()

    if (cas%type /= CASIDA_EPS_DIFF) then
      ! This is to be allocated here, and is used inside K_term.
      if(.not. cas%triplet) then
        SAFE_ALLOCATE(pot(1:mesh%np))
      endif
      qi_old = -1
      qa_old = -1
      mu_old = -1
      
      ! We calculate here the kernel, since it will be needed later.
      SAFE_ALLOCATE(rho(1:mesh%np, 1:st%d%nspin))
      if(cas%triplet) then
        SAFE_ALLOCATE(rho_spin(1:mesh%np, 1:2))
        SAFE_ALLOCATE(fxc(1:mesh%np, 1:2, 1:2))
      else
        SAFE_ALLOCATE(fxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
      endif
      fxc = M_ZERO

      call states_total_density(st, mesh, rho)
      if(cas%triplet) then
        rho_spin(:, 1) = M_HALF * rho(:, 1)
        rho_spin(:, 2) = M_HALF * rho(:, 1)
        call xc_get_fxc(sys%ks%xc, mesh, rho_spin, SPIN_POLARIZED, fxc)
      else
        call xc_get_fxc(sys%ks%xc, mesh, rho, st%d%ispin, fxc)
      endif
    end if

    select case(cas%type)
    case(CASIDA_EPS_DIFF,CASIDA_PETERSILKA)
      call solve_petersilka()
    case(CASIDA_TAMM_DANCOFF,CASIDA_VARIATIONAL,CASIDA_CASIDA)
      call solve_casida()
    end select

    ! clean up
    if (cas%type /= CASIDA_EPS_DIFF) then
      SAFE_DEALLOCATE_A(fxc)
      SAFE_DEALLOCATE_A(rho)
      if(.not. cas%triplet) then
        SAFE_DEALLOCATE_A(pot)
      endif
    end if
    SAFE_DEALLOCATE_A(saved_K)

    POP_SUB(casida_work)

  contains

    ! ---------------------------------------------------------
    subroutine solve_petersilka
      integer :: ia, iunit, idir
      FLOAT   :: ff
      FLOAT, allocatable :: deltav(:), xx(:)

      PUSH_SUB(casida_work.solve_petersilka)

      ! initialize progress bar
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, cas%n_pairs)

      ! file to save matrix elements
      iunit = io_open(trim(cas%restart_file), action='write', &
        position='append', is_tmp=.true.)

      do ia = 1, cas%n_pairs
        cas%w(ia) = st%eigenval(cas%pair(ia)%a, cas%pair(ia)%sigma) - &
                    st%eigenval(cas%pair(ia)%i, cas%pair(ia)%sigma)
        if(cas%w(ia) < -M_EPSILON) then
          message(1) = "There are negative unocc-occ KS eigenvalue differences."
          message(2) = "Probably this indicates an inconsistency in occupations between gs and unocc calculations."
          call messages_fatal(2)
        endif

        if(cas%type == CASIDA_PETERSILKA) then
          if(saved_K(ia, ia)) then
            ff = cas%mat(ia, ia)
          else
            cas%mat(ia, ia) = K_term(cas%pair(ia), cas%pair(ia))
            ff = cas%mat(ia, ia)
            call write_K_term(cas, iunit, ia, ia)
          end if

          cas%w(ia) = cas%w(ia) + cas%el_per_state * ff
          ! note that Petersilka is probably inappropriate for spin-polarized system due to degenerate transitions!
        end if

        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ia, cas%n_pairs)
      end do


      SAFE_ALLOCATE(xx(1:cas%n_pairs))
      SAFE_ALLOCATE(deltav(1:mesh%np))

      do idir = 1, mesh%sb%dim

        deltav(1:mesh%np) = mesh%x(1:mesh%np, idir)
        
        !WARNING: should xx always be real?
        if (states_are_real(st)) then
          xx = dks_matrix_elements(cas, st, mesh, deltav)
        else
          xx = zks_matrix_elements(cas, st, mesh, deltav)
        end if
        
        cas%tm(:, idir) = xx(:)
        cas%tm(:, idir) = sqrt(TOFLOAT(cas%el_per_state)) * cas%tm(:, idir) 

      end do
      SAFE_DEALLOCATE_A(xx)
      SAFE_DEALLOCATE_A(deltav)

      do ia = 1, cas%n_pairs
        cas%f(ia) = (M_TWO / mesh%sb%dim) * cas%w(ia) * sum((abs(cas%tm(ia, :)))**2)
      end do

      if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

      ! close restart file
      call io_close(iunit)
      POP_SUB(casida_work.solve_petersilka)
    end subroutine solve_petersilka


    ! ---------------------------------------------------------
    subroutine solve_casida()
      FLOAT :: temp
      integer :: ip, ia, jb, idir
      integer :: max, actual, iunit, counter
      FLOAT, allocatable :: deltav(:)
      CMPLX, allocatable :: zf(:)

      FLOAT, allocatable :: dx(:)
      CMPLX, allocatable :: zx(:)
      type(states_pair_t), pointer :: p, q

      FLOAT, allocatable :: gaus_leg_points(:), gaus_leg_weights(:)
      integer :: ii, jj
      FLOAT :: theta, phi, qlen
      FLOAT :: qvect(MAX_DIM)

      PUSH_SUB(casida_work.solve_casida)

      max = ceiling((cas%n_pairs*(M_ONE + cas%n_pairs)/M_TWO)/cas%mpi_grp%size)
      counter = 0
      actual = 0
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, max)

      if(.not.mpi_grp_is_root(mpi_world)) cas%mat = M_ZERO

      ! calculate the matrix elements of (v + fxc)
      do jb = 1, cas%n_pairs
        actual = actual + 1
        if(mod(actual, cas%mpi_grp%size) .ne. cas%mpi_grp%rank) cycle
        do ia = jb, cas%n_pairs
          counter = counter + 1
          ! if not loaded, then calculate matrix element
          if(.not.saved_K(ia, jb)) then
            cas%mat(ia, jb) = K_term(cas%pair(ia), cas%pair(jb))
          end if
          if(jb /= ia) cas%mat(jb, ia) = cas%mat(ia, jb) ! the matrix is symmetric
        end do
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(counter, max)
      end do
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(max, max)

      ! sum all matrix elements
      if(cas%parallel_in_eh_pairs) then
        call comm_allreduce(cas%mpi_grp%comm, cas%mat, dim = (/cas%n_pairs, cas%n_pairs/))
      end if
      !if(mpi_grp_is_root(cas%mpi_grp)) print *, "mat =", cas%mat

      ! all processors with the exception of the first are done
      if (mpi_grp_is_root(cas%mpi_grp)) then

        ! complete progress bar
        if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

        ! complete the matrix and output the restart file
        iunit = io_open(trim(cas%restart_file), action='write', &
          position='append', is_tmp=.true.)

        do ia = 1, cas%n_pairs
          p => cas%pair(ia)
          temp = st%eigenval(p%a, p%sigma) - st%eigenval(p%i, p%sigma)
            
          do jb = ia, cas%n_pairs
            q => cas%pair(jb)
            if(.not.saved_K(ia, jb)) call write_K_term(cas, iunit, ia, jb)
              
            if(cas%type == CASIDA_CASIDA) then
              cas%mat(ia, jb) = M_TWO * sqrt(temp) * cas%mat(ia, jb) * &
                sqrt(st%eigenval(q%a, q%sigma) - st%eigenval(q%i, q%sigma))
            else if(cas%type == CASIDA_VARIATIONAL) then
              cas%mat(ia, jb) = M_TWO * cas%mat(ia, jb)
            endif
            if(sys%st%d%ispin == UNPOLARIZED) then
              cas%mat(ia, jb) = M_TWO * cas%mat(ia, jb)
            endif

            if(jb /= ia) cas%mat(jb, ia) = cas%mat(ia, jb) ! the matrix is symmetric
          end do
          if(cas%type == CASIDA_CASIDA) then
            cas%mat(ia, ia) = temp**2 + cas%mat(ia, ia)
          else
            cas%mat(ia, ia) = cas%mat(ia, ia) + temp
          endif
        end do
        call io_close(iunit)

        ! now we diagonalize the matrix
        ! for huge matrices, perhaps we should consider ScaLAPACK here...
        call lalg_eigensolve(cas%n_pairs, cas%mat, cas%w)

        do ia = 1, cas%n_pairs
          if(cas%w(ia) < -M_EPSILON) then
            write(message(1),'(a,i4,a)') 'For whatever reason, excitation energy', ia, ' is negative.'
            write(message(2),'(a)')      'This should not happen.'
            call messages_warning(2)
            cas%w(ia) = M_ZERO
          else
            if(cas%type == CASIDA_CASIDA) cas%w(ia) = sqrt(cas%w(ia))
          end if
        end do

        ! And let us now get the S matrix...
        if(cas%type == CASIDA_CASIDA) then
          do ia = 1, cas%n_pairs
            if(sys%st%d%ispin == UNPOLARIZED) then
              cas%s(ia) = M_HALF / ( st%eigenval(cas%pair(ia)%a, 1) - st%eigenval(cas%pair(ia)%i, 1) )
            elseif(sys%st%d%ispin == SPIN_POLARIZED) then
              cas%s(ia) = M_ONE / ( st%eigenval(cas%pair(ia)%a, cas%pair(ia)%sigma) - &
                                    st%eigenval(cas%pair(ia)%i, cas%pair(ia)%sigma) )
            end if
          end do
        endif

        SAFE_ALLOCATE(deltav(1:mesh%np))
        if (states_are_real(st)) then

          if(cas%qcalc) then
             SAFE_ALLOCATE(zf(1:mesh%np))
             SAFE_ALLOCATE(zx(1:cas%n_pairs))

             ! matrix element
             do ia = 1, cas%n_pairs
               do ip = 1, mesh%np
                 zf(ip) = exp(M_zI * dot_product(cas%qvector(:), mesh%x(ip, :))) * &
                          st%dpsi(ip, 1, cas%pair(ia)%i, cas%pair(ia)%sigma) * &
                          st%dpsi(ip, 1, cas%pair(ia)%a, cas%pair(ia)%sigma)
               end do
               zx(ia) = zmf_integrate(mesh, zf)
             end do

             ! intensity
             do ia = 1, cas%n_pairs
               cas%qf(ia) = abs(ztransition_matrix_element(cas, ia, zx))**2
             end do

             ! do we calculate the average
             if(cas%avg_order .gt. 0) then

               ! use Gauss-Legendre quadrature scheme
               SAFE_ALLOCATE(gaus_leg_points (1:cas%avg_order))
               SAFE_ALLOCATE(gaus_leg_weights(1:cas%avg_order))
               call gauss_legendre_points(cas%avg_order, gaus_leg_points, gaus_leg_weights)

               qlen = sqrt(dot_product(cas%qvector, cas%qvector))
               do ii = 1, cas%avg_order
                 do jj = 1, 2 * cas%avg_order

                   ! construct the q-vector
                   phi   = acos(gaus_leg_points(ii))
                   theta = M_PI * jj / cas%avg_order
                   qvect(1) = qlen * cos(theta) * sin(phi)
                   qvect(2) = qlen * sin(theta) * sin(phi)
                   qvect(3) = qlen * cos(phi)

                   ! matrix elements
                   zx(:) = M_ZERO
                   zf(:) = M_ZERO
                   do ia = 1, cas%n_pairs
                     forall(ip = 1:mesh%np)
                       zf(ip) = exp(M_zI * dot_product(qvect(1:3), mesh%x(ip, 1:3))) * &
                         st%dpsi(ip, 1, cas%pair(ia)%i, cas%pair(ia)%sigma) * &
                         st%dpsi(ip, 1, cas%pair(ia)%a, cas%pair(ia)%sigma)
                     end forall
                     zx(ia) = zmf_integrate(mesh, zf)
                   end do

                   ! intensities
                   do ia = 1, cas%n_pairs
                     cas%qf_avg(ia) = cas%qf_avg(ia) + &
                                      gaus_leg_weights(ii)*abs(ztransition_matrix_element(cas, ia, zx))**2
                   end do

                 end do ! jj (thetas)
               end do ! ii (phis)

               ! normalize: for integral over sphere one would multiply by pi/N, but since
               !            we want the average, the integral must be divided by 4*pi
               forall(ia = 1:cas%n_pairs) cas%qf_avg(ia) = cas%qf_avg(ia) / (4*cas%avg_order)

               ! and finalize
               SAFE_DEALLOCATE_A(gaus_leg_points)
               SAFE_DEALLOCATE_A(gaus_leg_weights)
 
             end if ! averaging

             SAFE_DEALLOCATE_A(zf)
             SAFE_DEALLOCATE_A(zx)

          end if

          SAFE_ALLOCATE(dx(1:cas%n_pairs))
          do idir = 1, mesh%sb%dim
            deltav(1:mesh%np) = mesh%x(1:mesh%np, idir)
            ! let us get now the x vector.
            dx = dks_matrix_elements(cas, st, mesh, deltav)
            ! And now we are able to get the transition matrix elements between many-electron states.
            do ia = 1, cas%n_pairs
              cas%tm(ia, idir) = dtransition_matrix_element(cas, ia, dx)
            end do
          end do
          SAFE_DEALLOCATE_A(dx)
        else
          SAFE_ALLOCATE(zx(1:cas%n_pairs))
          do idir = 1, mesh%sb%dim
            deltav(1:mesh%np) = mesh%x(1:mesh%np, idir)
            ! let us get now the x vector.
            zx = zks_matrix_elements(cas, st, mesh, deltav)
            ! And now we are able to get the transition matrix elements between many-electron states.
            do ia = 1, cas%n_pairs
              cas%tm(ia, idir) = ztransition_matrix_element(cas, ia, zx)
            end do
          end do
          SAFE_DEALLOCATE_A(zx)
        end if
        SAFE_DEALLOCATE_A(deltav)


        ! And the oscillator strengths.
        do ia = 1, cas%n_pairs
          cas%f(ia) = (M_TWO / mesh%sb%dim) * cas%w(ia) * sum( (abs(cas%tm(ia, :)))**2 )
        end do

      end if

#if defined(HAVE_MPI)
      if(cas%parallel_in_eh_pairs) then
        call MPI_Barrier(cas%mpi_grp%comm, mpi_err)
      end if
#endif

      POP_SUB(casida_work.solve_casida)
    end subroutine solve_casida


    ! ---------------------------------------------------------
    ! return the matrix element of <i(p),a(p)|v + fxc|j(q),b(q)>
    FLOAT function K_term(pp, qq)
      type(states_pair_t), intent(in) :: pp, qq

      integer :: pi, qi, sigma, pa, qa, mu
      FLOAT, allocatable :: rho_i(:), rho_j(:)

      PUSH_SUB(casida_work.K_term)

      pi = pp%i
      pa = pp%a
      sigma = pp%sigma

      qi = qq%i
      qa = qq%a
      mu = qq%sigma

      SAFE_ALLOCATE(rho_i(1:mesh%np))
      SAFE_ALLOCATE(rho_j(1:mesh%np))

      if (states_are_real(st)) then
        rho_i(1:mesh%np) =        st%dpsi(1:mesh%np, 1, pi, sigma) *       st%dpsi(1:mesh%np, 1, pa, sigma)
        rho_j(1:mesh%np) =        st%dpsi(1:mesh%np, 1, qi, mu)    *       st%dpsi(1:mesh%np, 1, qa, mu)
      else
        rho_i(1:mesh%np) =        st%zpsi(1:mesh%np, 1, pi, sigma) * conjg(st%zpsi(1:mesh%np, 1, pa, sigma))
        rho_j(1:mesh%np) =  conjg(st%zpsi(1:mesh%np, 1, qi, mu))   *       st%zpsi(1:mesh%np, 1, qa, mu)
      end if

      !  first the Hartree part (only works for real wfs...)
      if(.not. cas%triplet .and. (qi .ne. qi_old  .or.   qa .ne. qa_old   .or.  mu .ne. mu_old)) then
        pot(1:mesh%np) = M_ZERO
        if(hm%theory_level .ne. INDEPENDENT_PARTICLES) &
          call dpoisson_solve(psolver, pot, rho_j, all_nodes=.false.)
      end if

      if(cas%triplet) then
        K_term = M_ZERO
        rho(1:mesh%np, 1) = rho_i(1:mesh%np) * rho_j(1:mesh%np) * M_HALF * (fxc(1:mesh%np, 1, 1) - fxc(1:mesh%np, 1, 2))
      else
        K_term = dmf_dotp(mesh, rho_i(:), pot(:))
        rho(1:mesh%np, 1) = rho_i(1:mesh%np) * rho_j(1:mesh%np) * fxc(1:mesh%np, sigma, mu)
      endif
      K_term = K_term + dmf_integrate(mesh, rho(:, 1))

      qi_old = qi
      qa_old = qa
      mu_old = mu

      SAFE_DEALLOCATE_A(rho_i)
      SAFE_DEALLOCATE_A(rho_j)

      POP_SUB(casida_work.K_term)
    end function K_term

    ! ---------------------------------------------------------
    subroutine load_saved
      integer :: iunit, err
      integer :: ia, jb, ii, aa, is, jj, bb, js
      FLOAT   :: val

      PUSH_SUB(casida_work.load_saved)

      iunit = io_open(trim(cas%restart_file), action='read', &
        status='old', die=.false., is_tmp=.true.)
      if( iunit <= 0) then
        POP_SUB(casida_work.load_saved)
        return
      end if

      do
        read(iunit, fmt=*, iostat=err) ii, aa, is, jj, bb, js, val
        if(err.ne.0) exit

        ia = cas%index(ii, aa, is)
        jb = cas%index(jj, bb, js)

        if(ia > 0 .and. jb > 0) then
          cas%mat(ia, jb) = val
          saved_K(ia, jb) = .true.
          cas%mat(jb, ia) = val
          saved_K(jb, ia) = .true.
        endif
      end do

      if(iunit > 0) call io_close(iunit)
      POP_SUB(casida_work.load_saved)
    end subroutine load_saved

  end subroutine casida_work

  ! ---------------------------------------------------------
  subroutine qcasida_write(cas)
    type(casida_t), intent(in) :: cas

    integer :: iunit, ia, dim
    integer, allocatable :: ind(:)
    FLOAT, allocatable :: w(:)

    if(.not.mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(qcasida_write)

    dim = size(cas%tm, 2)

    SAFE_ALLOCATE(  w(1:cas%n_pairs))
    SAFE_ALLOCATE(ind(1:cas%n_pairs))
    w = cas%w
    call sort(w, ind)

    call io_mkdir(CASIDA_DIR)
    iunit = io_open(CASIDA_DIR//'q'//trim(theory_name(cas)), action='write')
    write(iunit, '(a1,a14,1x,a24,1x,a24,1x,a10,3es15.8,a2)') '#','E' , '|<f|exp(iq.r)|i>|^2', &
                                                             '<|<f|exp(iq.r)|i>|^2>','; q = (',cas%qvector(1:MAX_DIM),')'
    write(iunit, '(a1,a14,1x,a24,1x,a24,1x,10x,a15)')        '#', trim(units_abbrev(units_out%energy)), &
                                                                  trim('-'), &
                                                                  trim('-'), &
                                                                  trim('a.u.')

    if(cas%avg_order.eq.0) then
      do ia = 1, cas%n_pairs
        write(iunit, '(es15.8,es15.8)') units_from_atomic(units_out%energy, cas%w(ind(ia))), cas%qf(ind(ia))
      end do
    else
      do ia = 1, cas%n_pairs
        write(iunit, '(3es15.8)') units_from_atomic(units_out%energy, cas%w(ind(ia))), &
                                  cas%qf    (ind(ia)), &
                                  cas%qf_avg(ind(ia))
      end do
    end if

    call io_close(iunit)

    SAFE_DEALLOCATE_A(w)
    SAFE_DEALLOCATE_A(ind)

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
    integer, allocatable :: ind(:)
    FLOAT, allocatable :: w(:)

    if(.not.mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(casida_write)

    dim = size(cas%tm, 2)

    SAFE_ALLOCATE(  w(1:cas%n_pairs))
    SAFE_ALLOCATE(ind(1:cas%n_pairs))
    w = cas%w
    call sort(w, ind)

    ! output excitation energies and oscillator strengths
    call io_mkdir(CASIDA_DIR)
    iunit = io_open(CASIDA_DIR//trim(theory_name(cas)), action='write')

    if(cas%type == CASIDA_EPS_DIFF .or. (cas%type == CASIDA_PETERSILKA)) then
      write(iunit, '(2a4)', advance='no') 'From', '  To'
      if(cas%nspin == 2) then
        write(iunit, '(a5)', advance='no') 'Spin'
      endif
    else
      write(iunit, '(6x)', advance='no')
    endif

    select case(dim)
    case(1); write(iunit, '(3(1x,a15))') 'E' , '<x>', '<f>'
    case(2); write(iunit, '(4(1x,a15))') 'E' , '<x>', '<y>', '<f>'
    case(3); write(iunit, '(5(1x,a15))') 'E' , '<x>', '<y>', '<z>', '<f>'
    end select
    do ia = 1, cas%n_pairs
      if((cas%type == CASIDA_EPS_DIFF) .or. (cas%type == CASIDA_PETERSILKA)) then
        write(iunit, '(2i4)', advance='no') cas%pair(ind(ia))%i, cas%pair(ind(ia))%a
        if(cas%nspin == 2) then
          write(iunit, '(i5)', advance='no') cas%pair(ind(ia))%sigma
        endif
      else
        write(iunit, '(i6)', advance='no') ind(ia)
      end if
      write(iunit, '(5(1x,es15.8))') units_from_atomic(units_out%energy, cas%w(ind(ia))), &
        (units_from_atomic(units_out%length, cas%tm(ind(ia), idim)), idim=1,dim), cas%f(ind(ia))
    end do
    call io_close(iunit)

    if(cas%qcalc) call qcasida_write(cas)

    if(cas%type == CASIDA_EPS_DIFF .or. cas%type == CASIDA_PETERSILKA) then
      POP_SUB(casida_write)
      return
    end if

    ! output eigenvectors
    dir_name = CASIDA_DIR//trim(theory_name(cas))//'_excitations'
    call io_mkdir(trim(dir_name))
    do ia = 1, cas%n_pairs
      write(str,'(i5.5)') ia
      iunit = io_open(trim(dir_name)//'/'//trim(str), action='write')
      ! First, a little header
      write(iunit,'(a,es14.5)') '# Energy ['// trim(units_abbrev(units_out%energy)) // '] = ', &
                                units_from_atomic(units_out%energy, cas%w(ind(ia)))
        write(iunit,'(a,es14.5)') '# <X> ['//trim(units_abbrev(units_out%length))// '] = ', &
                                  units_from_atomic(units_out%length, cas%tm(ind(ia),1))
      if(dim > 1) &
        write(iunit,'(a,es14.5)') '# <Y> ['//trim(units_abbrev(units_out%length))// '] = ', &
                                  units_from_atomic(units_out%length, cas%tm(ind(ia),2))
      if(dim > 2) &
        write(iunit,'(a,es14.5)') '# <Z> ['//trim(units_abbrev(units_out%length))// '] = ', &
                                  units_from_atomic(units_out%length, cas%tm(ind(ia),3))

      temp = M_ONE
      ! make the largest component positive, to specify the phase
      if( maxval(cas%mat(:, ind(ia))) - abs(minval(cas%mat(:, ind(ia)))) < -M_EPSILON) temp = -temp

      do jb = 1, cas%n_pairs
        write(iunit,*) cas%pair(jb)%i, cas%pair(jb)%a, cas%pair(jb)%sigma, temp * cas%mat(jb, ind(ia))
      end do

      if(cas%type == CASIDA_TAMM_DANCOFF .or. cas%type == CASIDA_VARIATIONAL) then
        call write_implied_occupations(cas, iunit, ind(ia))
      endif
      call io_close(iunit)
    end do

    SAFE_DEALLOCATE_A(w)
    SAFE_DEALLOCATE_A(ind)

    ! Calculate and write the transition densities
    if (states_are_real(sys%st)) then
      call dget_transition_densities(cas, sys)
    else
      call zget_transition_densities(cas, sys)
    end if

    POP_SUB(casida_write)
  end subroutine casida_write

  ! ---------------------------------------------------------
  ! write matrix element to casida_restart file
  subroutine write_K_term(cas, iunit, ia, jb)
    type(casida_t), intent(in) :: cas
    integer,        intent(in) :: iunit
    integer,        intent(in) :: ia
    integer,        intent(in) :: jb

    PUSH_SUB(write_K_term)

    write(iunit,*) cas%pair(ia)%i, cas%pair(ia)%a, cas%pair(ia)%sigma, &
                   cas%pair(jb)%i, cas%pair(jb)%a, cas%pair(jb)%sigma, cas%mat(ia, jb)

    POP_SUB(write_K_term)
  end subroutine write_K_term

  ! ---------------------------------------------------------
  character*80 function theory_name(cas)
    type(casida_t), intent(in) :: cas

    PUSH_SUB(theory_name)

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
        write(message(1),'(a,i6)') 'Unknown Casida theory level ', cas%type
        call messages_fatal(1)
    end select

    POP_SUB(theory_name)
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

    do ik = 1, cas%nspin
      do ist = 1, cas%n_occ(ik)
        occ = M_ONE * cas%el_per_state
        do ast = cas%n_occ(ik) + 1, cas%n_occ(ik) + cas%n_unocc(ik)
          if(cas%index(ist, ast, ik) == 0) cycle  ! we were not using this state
          occ = occ - abs(cas%mat(cas%index(ist, ast, ik), ind))**2
        enddo
        write(iunit, '(f8.6,a)', advance='no') occ, ' | '
      enddo
      do ast = cas%n_occ(ik) + 1, cas%n_occ(ik) + cas%n_unocc(ik)
        occ = M_ZERO
        do ist = 1, cas%n_occ(ik)
          if(cas%index(ist, ast, ik) == 0) cycle  ! we were not using this state
          occ = occ + abs(cas%mat(cas%index(ist, ast, ik), ind))**2
        enddo
        write(iunit, '(f8.6)', advance='no') occ
        if(ast < cas%n_occ(ik) + cas%n_unocc(ik)) write(iunit, '(a)', advance='no') ' | '
      enddo
      write(iunit, '(a)')
    enddo
    
    write(iunit, '(a)') '%'

    POP_SUB(write_implied_occupations)
  end subroutine write_implied_occupations

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
