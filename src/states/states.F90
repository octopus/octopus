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

module states_m
  use blacs_proc_grid_m
  use calc_mode_m
  use comm_m
  use batch_m
  use blas_m
  use datasets_m
  use derivatives_m
  use distributed_m
  use geometry_m
  use global_m
  use grid_m
  use hardware_m
  use io_m
  use kpoints_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use modelmb_particles_m
  use mpi_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2
  use mpi_lib_m
  use multicomm_m
  use parser_m
  use profiling_m
  use simul_box_m
  use smear_m
  use states_dim_m
  use symmetrizer_m
  use types_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m

  implicit none

  private

  public ::                           &
    states_t,                         &
    states_priv_t,                    &
    states_lead_t,                    &
    states_init,                      &
    states_look,                      &
    states_densities_init,            &
    states_allocate_wfns,             &
    states_deallocate_wfns,           &
    states_null,                      &
    states_end,                       &
    states_copy,                      &
    states_generate_random,           &
    states_fermi,                     &
    states_eigenvalues_sum,           &
    states_allocate_free_states,      &
    states_deallocate_free_states,    &
    states_spin_channel,              &
    states_calc_quantities,           &
    state_is_local,                   &
    states_distribute_nodes,          &
    states_wfns_memory,               &
    states_init_block,                &
    states_are_complex,               &
    states_are_real,                  &
    states_set_complex,               &
    states_blacs_blocksize

  type states_lead_t
    CMPLX, pointer     :: intf_psi(:, :, :, :) !< (np, st%d%dim, st%nst, st%d%nik)
    FLOAT, pointer     :: rho(:, :)   !< Density of the lead unit cells.
    CMPLX, pointer     :: self_energy(:, :, :, :, :) !< (np, np, nspin, ncs, nik) self-energy of the leads.
  end type states_lead_t

  type states_priv_t
    private
    type(type_t) :: wfs_type              !< real (TYPE_FLOAT) or complex (TYPE_CMPLX) wavefunctions
  end type states_priv_t

  type states_t
    type(states_dim_t)       :: d
    type(modelmb_particle_t) :: modelmbparticles
    type(states_priv_t)      :: priv                  !< the private components 
    integer                  :: nst                   !< Number of states in each irreducible subspace

    ! pointers to the wavefunctions 
    logical                  :: only_userdef_istates  !< only use user-defined states as initial states in propagation
    FLOAT, pointer           :: dpsi(:,:,:,:)         !< dpsi(sys%gr%mesh%np_part, st%d%dim, st%nst, st%d%nik)
    CMPLX, pointer           :: zpsi(:,:,:,:)         !< zpsi(sys%gr%mesh%np_part, st%d%dim, st%nst, st%d%nik)

    type(batch_t), pointer   :: psib(:, :)            !< A set of wave-functions blocks
    integer                  :: nblocks               !< The number of blocks
    integer, pointer         :: iblock(:, :)          !< A map, that for each state index, returns the index of block that contains it. 
    logical, pointer         :: block_is_local(:, :)  !< It is true if the block is in this node.
    logical                  :: block_initialized     !< For keeping track of the blocks to avoid memory leaks

    logical             :: open_boundaries
    CMPLX, pointer      :: zphi(:, :, :, :)  !< Free states for open-boundary calculations.
    FLOAT, pointer      :: ob_eigenval(:, :) !< Eigenvalues of free states.
    type(states_dim_t)  :: ob_d              !< Dims. of the unscattered systems.
    integer             :: ob_nst            !< nst of the unscattered systems.
    FLOAT, pointer      :: ob_occ(:, :)      !< occupations
    type(states_lead_t) :: ob_lead(2*MAX_DIM)

    !> used for the user-defined wavefunctions (they are stored as formula strings)
    !! (st%d%dim, st%nst, st%d%nik)
    character(len=1024), pointer :: user_def_states(:,:,:)

    !> the densities and currents (after all we are doing DFT :)
    FLOAT, pointer :: rho(:,:)         !< rho(gr%mesh%np_part, st%d%nspin)
    FLOAT, pointer :: current(:, :, :) !<   current(gr%mesh%np_part, gr%sb%dim, st%d%nspin)

    logical        :: nlcc             !< do we have non-linear core corrections
    FLOAT, pointer :: rho_core(:)      !< core charge for nl core corrections
    logical        :: current_in_tau   !< are we using in tau the term which depends on the paramagnetic current?  
    
    !> It may be required to "freeze" the deepest orbitals during the evolution; the density
    !! of these orbitals is kept in frozen_rho. It is different from rho_core.
    FLOAT, pointer :: frozen_rho(:, :)

    FLOAT, pointer :: eigenval(:,:) !< obviously the eigenvalues
    logical        :: fixed_occ     !< should the occupation numbers be fixed?
    FLOAT, pointer :: occ(:,:)      !< the occupation numbers
    logical        :: fixed_spins   !< In spinors mode, the spin direction is set
                                    !< for the initial (random) orbitals.
    FLOAT, pointer :: spin(:, :, :)

    FLOAT          :: qtot          !< (-) The total charge in the system (used in Fermi)
    FLOAT          :: val_charge    !< valence charge

    logical        :: extrastates   ! are there extra states?
    type(smear_t)  :: smear         ! smearing of the electronic occupations

    !> This is stuff needed for the parallelization in states.
    logical                     :: parallel_in_states !< Am I parallel in states?
    type(mpi_grp_t)             :: mpi_grp            !< The MPI group related to the parallelization in states.
    type(mpi_grp_t)             :: dom_st_mpi_grp     !< The MPI group related to the domain-states "plane".
    type(mpi_grp_t)             :: st_kpt_mpi_grp     !< The MPI group related to the states-kpoints "plane".
    type(mpi_grp_t)             :: dom_st_kpt_mpi_grp !< The MPI group related to the domains-states-kpoints "cube".
#ifdef HAVE_SCALAPACK
    type(blacs_proc_grid_t)     :: dom_st_proc_grid   !< The BLACS process grid for the domains states plane
#endif
    integer                     :: lnst               !< Number of states on local node.
    integer                     :: st_start, st_end   !< Range of states processed by local node.
    integer, pointer            :: node(:)            !< To which node belongs each state.
    integer, pointer            :: st_range(:, :)     !< Node r manages states st_range(1, r) to
                                                      !! st_range(2, r) for r = 0, ..., mpi_grp%size-1,
                                                      !! i. e. st_start = st_range(1, r) and
                                                      !! st_end = st_range(2, r) on node r.
    integer, pointer            :: st_num(:)          !< Number of states on node r, i. e.
                                                      !! st_num(r) = st_num(2, r)-st_num(1, r).
    type(multicomm_all_pairs_t) :: ap                 !< All-pairs schedule.

    logical                     :: symmetrize_density
  end type states_t

contains

  ! ---------------------------------------------------------
  subroutine states_null(st)
    type(states_t), intent(inout) :: st

    integer :: il

    PUSH_SUB(states_null)

    nullify(st%dpsi, st%zpsi, st%zphi, st%rho, st%current, st%rho_core, st%frozen_rho, st%eigenval)
    do il=1, NLEADS
      nullify(st%ob_lead(il)%intf_psi, st%ob_lead(il)%rho, st%ob_lead(il)%self_energy)
    end do
    nullify(st%ob_eigenval, st%ob_occ)
    nullify(st%occ, st%spin, st%node, st%user_def_states)
    nullify(st%st_range, st%st_num)
    nullify(st%psib, st%iblock, st%block_is_local)
    nullify(st%ap%schedule)

    st%parallel_in_states = .false.

    ! By default, calculations use real wavefunctions
    st%priv%wfs_type = TYPE_FLOAT

    st%block_initialized = .false.
    call states_dim_null(st%d)
    call states_dim_null(st%ob_d)
    call modelmb_particles_nullify(st%modelmbparticles)
#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_nullify(st%dom_st_proc_grid)
#endif
    st%d%orth_method = 0

    POP_SUB(states_null)
  end subroutine states_null


  ! ---------------------------------------------------------
  subroutine states_init(st, gr, geo)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo

    FLOAT :: excess_charge
    integer :: nempty, ierr, il, ntot
    integer, allocatable :: ob_k(:), ob_st(:), ob_d(:)
    character(len=256)   :: restart_dir

    PUSH_SUB(states_init)

    call states_null(st)

    !%Variable SpinComponents
    !%Type integer
    !%Default unpolarized
    !%Section States
    !%Description
    !% The calculations may be done in three different ways: spin-restricted (TD)DFT (<i>i.e.</i>, doubly
    !% occupied "closed shells"), spin-unrestricted or "spin-polarized" (TD)DFT (<i>i.e.</i> we have two
    !% electronic systems, one with spin up and one with spin down), or making use of two-component
    !% spinors.
    !%Option unpolarized 1
    !% Spin-restricted calculations.
    !%Option polarized 2
    !%Option spin_polarized 2
    !% Spin unrestricted, also known as spin-DFT, SDFT. This mode will double the number of 
    !% wavefunctions necessary for a spin-unpolarized calculation.
    !%Option non_collinear 3
    !%Option spinors 3
    !% The spin-orbitals are two-component spinors. This effectively allows the spin-density to
    !% be oriented non-collinearly: <i>i.e.</i> the magnetization vector is allowed to take different
    !% directions at different points. This vector is always in 3D regardless of <tt>Dimensions</tt>.
    !%End
    call parse_integer(datasets_check('SpinComponents'), UNPOLARIZED, st%d%ispin)
    if(.not.varinfo_valid_option('SpinComponents', st%d%ispin)) call input_error('SpinComponents')
    call messages_print_var_option(stdout, 'SpinComponents', st%d%ispin)
    ! Use of spinors requires complex wavefunctions.
    if (st%d%ispin == SPINORS) st%priv%wfs_type = TYPE_CMPLX

    !%Variable ExcessCharge
    !%Type float
    !%Default 0.0
    !%Section States
    !%Description
    !% The net charge of the system. A negative value means that we are adding 
    !% electrons, while a positive value means we are taking electrons
    !% from the system.
    !%End
    call parse_float(datasets_check('ExcessCharge'), M_ZERO, excess_charge)


    !%Variable TotalStates
    !%Type integer
    !%Default 0
    !%Section States
    !%Description
    !% This variable sets the total number of states that Octopus will
    !% use. This is normally not necessary since by default Octopus
    !% sets the number of states to the minimum necessary to hold the
    !% electrons present in the system. (This default behavior is
    !% obtained by setting <tt>TotalStates</tt>  to 0).
    !%
    !% If you want to add some unoccupied states, probably it is more convenient to use the variable 
    !% <tt>ExtraStates</tt>.
    !%
    !% Note that this number is unrelated to <tt>CalculationMode == unocc</tt>.
    !%End
    call parse_integer(datasets_check('TotalStates'), 0, ntot)
    if (ntot < 0) then
      write(message(1), '(a,i5,a)') "Input: '", ntot, "' is not a valid value for TotalStates."
      call messages_fatal(1)
    end if

    !%Variable ExtraStates
    !%Type integer
    !%Default 0
    !%Section States
    !%Description
    !% The number of states is in principle calculated considering the minimum
    !% numbers of states necessary to hold the electrons present in the system.
    !% The number of electrons is
    !% in turn calculated considering the nature of the species supplied in the
    !% <tt>Species</tt> block, and the value of the <tt>ExcessCharge</tt> variable.
    !% However, one may command <tt>Octopus</tt> to use more states, which is necessary if one wants to
    !% use fractional occupational numbers, either fixed from the beginning through
    !% the <tt>Occupations</tt> block or by prescribing
    !% an electronic temperature with <tt>Smearing</tt>.
    !%
    !% Note that this number is unrelated to <tt>CalculationMode == unocc</tt>.
    !% <tt>ExtraStates</tt> is used for a self-consistent calculation and
    !% the usual convergence criteria on the density do not take into account the
    !% eigenvalues, whereas <tt>unocc</tt> is a non-self-consistent calculation,
    !% and explicitly considers the eigenvalues of the unoccupied states as the
    !% convergence criteria.
    !%End
    call parse_integer(datasets_check('ExtraStates'), 0, nempty)
    if (nempty < 0) then
      write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid value for ExtraStates."
      message(2) = '(0 <= ExtraStates)'
      call messages_fatal(2)
    end if

    st%extrastates = (nempty > 0)

    if(ntot > 0 .and. nempty > 0) then
      message(1) = 'Error: You cannot set TotalStates and ExtraStates at the same time.'
      call messages_fatal(1)
    end if

    ! For non-periodic systems this should just return the Gamma point
    call states_choose_kpoints(st%d, gr%sb, geo)

    call geometry_val_charge(geo, st%val_charge)
    
    if(gr%ob_grid%open_boundaries) then
      ! renormalize charge of central region to match leads (open system, not finite)
      st%val_charge = st%val_charge * (gr%ob_grid%lead(LEFT)%sb%lsize(TRANS_DIR) / gr%sb%lsize(TRANS_DIR))
    end if

    st%qtot = -(st%val_charge + excess_charge)

    do il = 1, NLEADS
      nullify(st%ob_lead(il)%intf_psi)
    end do
    ! When doing open-boundary calculations the number of free states is
    ! determined by the previous periodic calculation.
    st%open_boundaries = gr%ob_grid%open_boundaries
    if(gr%ob_grid%open_boundaries) then
      SAFE_ALLOCATE( ob_k(1:NLEADS))
      SAFE_ALLOCATE(ob_st(1:NLEADS))
      SAFE_ALLOCATE( ob_d(1:NLEADS))
      do il = 1, NLEADS
        restart_dir = trim(trim(gr%ob_grid%lead(il)%info%restart_dir)//'/'// GS_DIR)
        ! first get nst and kpoints of all states
        call states_look(restart_dir, mpi_world, ob_k(il), ob_d(il), ob_st(il), ierr)
        if(ierr.ne.0) then
          message(1) = 'Could not read the number of states of the periodic calculation'
          message(2) = 'from '//restart_dir//'.'
          call messages_fatal(2)
        end if
      end do
      if(NLEADS.gt.1) then
        if(ob_k(LEFT).ne.ob_k(RIGHT).or. &
          ob_st(LEFT).ne.ob_st(LEFT).or. &
          ob_d(LEFT).ne.ob_d(RIGHT)) then
          message(1) = 'The number of states for the left and right leads are not equal.'
          call messages_fatal(1)
        end if
      end if
      st%ob_d%dim = ob_d(LEFT)
      st%ob_nst   = ob_st(LEFT)
      st%ob_d%nik = ob_k(LEFT)
      st%d%nik = st%ob_d%nik
      SAFE_DEALLOCATE_A(ob_d)
      SAFE_DEALLOCATE_A(ob_st)
      SAFE_DEALLOCATE_A(ob_k)
      call distributed_nullify(st%ob_d%kpt, 0)
      if((st%d%ispin.eq.UNPOLARIZED.and.st%ob_d%dim.ne.1) .or.   &
        (st%d%ispin.eq.SPIN_POLARIZED.and.st%ob_d%dim.ne.1) .or. &
        (st%d%ispin.eq.SPINORS.and.st%ob_d%dim.ne.2)) then
        message(1) = 'The spin type of the leads calculation from '&
                     //gr%ob_grid%lead(LEFT)%info%restart_dir
        message(2) = 'and SpinComponents of the current run do not match.'
        call messages_fatal(2)
      end if
      SAFE_DEALLOCATE_P(st%d%kweights)
      SAFE_ALLOCATE(st%d%kweights(1:st%d%nik))
      st%d%kweights = M_ZERO
      st%d%kweights(1) = M_ONE
      SAFE_ALLOCATE(st%ob_d%kweights(1:st%ob_d%nik))
      SAFE_ALLOCATE(st%ob_eigenval(1:st%ob_nst, 1:st%ob_d%nik))
      SAFE_ALLOCATE(st%ob_occ(1:st%ob_nst, 1:st%ob_d%nik))
      st%ob_d%kweights = M_ZERO
      st%ob_eigenval   = huge(st%ob_eigenval)
      st%ob_occ        = M_ZERO
      call read_ob_eigenval_and_occ()
    else
      st%ob_nst   = 0
      st%ob_d%nik = 0
      st%ob_d%dim = 0
    end if

    select case(st%d%ispin)
    case(UNPOLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%d%nspin = 1
      st%d%spin_channels = 1
    case(SPIN_POLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%d%nspin = 2
      st%d%spin_channels = 2
    case(SPINORS)
      st%d%dim = 2
      st%nst = int(st%qtot)
      if(st%nst < st%qtot) st%nst = st%nst + 1
      st%d%nspin = 4
      st%d%spin_channels = 2
    end select

    if(ntot > 0) then
      if(ntot < st%nst) then
        message(1) = 'Error: TotalStates is smaller than the number of states required by the system.'
        call messages_fatal(1)
      end if

      st%extrastates = (ntot > st%nst)
      st%nst = ntot
    end if

    st%nst = st%nst + nempty


    ! FIXME: For now, open-boundary calculations are only possible for
    ! continuum states, i.e. for those states treated by the Lippmann-
    ! Schwinger approach during SCF.
    ! Bound states should be done with extra states, without k-points.
    if(gr%ob_grid%open_boundaries) then
      if(st%nst.ne.st%ob_nst .or. st%d%nik.ne.st%ob_d%nik) then
        message(1) = 'Open-boundary calculations for possibly bound states'
        message(2) = 'are not possible yet. You have to match your number'
        message(3) = 'of states to the number of free states of your previous'
        message(4) = 'periodic run.'
        write(message(5), '(a,i5,a)') 'Your central region contributes ', st%nst, ' states,'
        write(message(6), '(a,i5,a)') 'while your lead calculation had ', st%ob_nst, ' states.'
        write(message(7), '(a,i5,a)') 'Your central region contributes ', st%d%nik, ' k-points,'
        write(message(8), '(a,i5,a)') 'while your lead calculation had ', st%ob_d%nik, ' k-points.'
        call messages_fatal(8)
      end if
    end if

    !%Variable CurrentDFT
    !%Type logical
    !%Default false
    !%Section Hamiltonian
    !%Description
    !% (experimental) If set to yes, Current-DFT will be used. This is the
    !% extension to DFT that should be used when external magnetic fields are
    !% present. The current-dependent part of the XC functional is set using the
    !% <tt>JFunctional</tt> variable. The default is no.
    !%End
    call parse_logical(datasets_check('CurrentDFT'), .false., st%d%cdft)
    if (st%d%cdft) then
      call messages_experimental('Current DFT')

      ! Use of CDFT requires complex wavefunctions
      st%priv%wfs_type = TYPE_CMPLX

      if(st%d%ispin == SPINORS) then
        message(1) = "Sorry, current DFT not working yet for spinors."
        call messages_fatal(1)
      end if
      message(1) = "Info: Using current DFT"
      call messages_info(1)
    end if

    ! Periodic systems require complex wavefunctions
    ! but not if it is Gamma-point only
    if(simul_box_is_periodic(gr%sb)) then
      if(.not. (kpoints_number(gr%sb%kpoints) == 1 .and. kpoints_point_is_gamma(gr%sb%kpoints, 1))) then
        st%priv%wfs_type = TYPE_CMPLX
      endif
    endif

    ! Calculations with open boundaries require complex wavefunctions.
    if(gr%ob_grid%open_boundaries) st%priv%wfs_type = TYPE_CMPLX

    !%Variable OnlyUserDefinedInitialStates
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% If true, then only user-defined states from the block <tt>UserDefinedStates</tt>
    !% will be used as initial states for a time-propagation. No attempt is made
    !% to load ground-state orbitals from a previous ground-state run.
    !%End
    call parse_logical(datasets_check('OnlyUserDefinedInitialStates'), .false., st%only_userdef_istates)

    !%Variable CurrentInTau
    !%Type logical
    !%Default yes
    !%Section States
    !%Description
    !% If true, a term including the (paramagnetic or total) current is included in the calculation ot the kinetic energy density 
    !%End
    call parse_logical(datasets_check('CurrentInTau'), .true., st%current_in_tau)


    ! we now allocate some arrays
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)
    st%occ      = M_ZERO
    ! allocate space for formula strings that define user-defined states
    SAFE_ALLOCATE(st%user_def_states(1:st%d%dim, 1:st%nst, 1:st%d%nik))
    if(st%d%ispin == SPINORS) then
      SAFE_ALLOCATE(st%spin(1:3, 1:st%nst, 1:st%d%nik))
    else
      nullify(st%spin)
    end if

    ! initially we mark all 'formulas' as undefined
    st%user_def_states(1:st%d%dim, 1:st%nst, 1:st%d%nik) = 'undefined'

    call states_read_initial_occs(st, excess_charge)
    call states_read_initial_spins(st)

    nullify(st%zphi)

    st%st_start = 1
    st%st_end = st%nst
    st%lnst = st%nst
    SAFE_ALLOCATE(st%node(1:st%nst))
    st%node(1:st%nst) = 0

    call mpi_grp_init(st%mpi_grp, -1)
    st%parallel_in_states = .false.

    nullify(st%dpsi, st%zpsi)

    call distributed_nullify(st%d%kpt, st%d%nik)

    call modelmb_particles_init (st%modelmbparticles,gr)

    !%Variable SymmetrizeDensity
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% (experimental) If true and the system is periodic, the density will be symmetrized.
    !%End
    call parse_logical(datasets_check('SymmetrizeDensity'), .false., st%symmetrize_density)

    if(st%symmetrize_density) call messages_experimental("Symmetrization of the density")

#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_nullify(st%dom_st_proc_grid)
#endif

    POP_SUB(states_init)

  contains

    subroutine read_ob_eigenval_and_occ()
      integer            :: occs, ist, ik, idim, idir, err
      FLOAT              :: flt, eigenval, occ, kweights
      character          :: char
      character(len=256) :: restart_dir, line, chars

      PUSH_SUB(states_init.read_ob_eigenval_and_occ)

      restart_dir = trim(gr%ob_grid%lead(LEFT)%info%restart_dir)//'/'//GS_DIR

      occs = io_open(trim(restart_dir)//'/occs', action='read', is_tmp=.true., grp=mpi_world)
      if(occs .lt. 0) then
        message(1) = 'Could not read '//trim(restart_dir)//'/occs.'
        call messages_fatal(1)
      end if

      ! Skip two lines.
      call iopar_read(mpi_world, occs, line, err)
      call iopar_read(mpi_world, occs, line, err)

      do
        ! Check for end of file.
        call iopar_read(mpi_world, occs, line, err)

        read(line, '(a)') char
        if(char .eq. '%') exit
        call iopar_backspace(mpi_world, occs)

        ! Extract eigenvalue.
        call iopar_read(mpi_world, occs, line, err)
        ! # occupations | eigenvalue[a.u.] | k-points | k-weights | filename | ik | ist | idim
        read(line, *) occ, char, eigenval, char, (flt, char, idir = 1, gr%sb%dim), kweights, &
           char, chars, char, ik, char, ist, char, idim

        if(st%d%ispin .eq. SPIN_POLARIZED) then
            message(1) = 'Spin-Transport not implemented!'
            call messages_fatal(1)
          if(is_spin_up(ik)) then
            !FIXME
!              st%ob_eigenval(jst, SPIN_UP) = eigenval
!              st%ob_occ(jst, SPIN_UP)      = occ
          else
!              st%ob_eigenval(jst, SPIN_DOWN) = eigenval
!              st%ob_occ(jst, SPIN_DOWN)      = occ
          end if
        else
          st%ob_eigenval(ist, ik) = eigenval
          st%ob_occ(ist, ik)      = occ
          st%ob_d%kweights(ik)    = kweights
        end if
      end do

      call io_close(occs)

      POP_SUB(states_init.read_ob_eigenval_and_occ)
    end subroutine read_ob_eigenval_and_occ
  end subroutine states_init

  ! ---------------------------------------------------------
  !> Reads the state stored in directory "dir", and finds out
  !! the kpoints, dim, and nst contained in it.
  ! ---------------------------------------------------------
  subroutine states_look(dir, mpi_grp, kpoints, dim, nst, ierr)
    character(len=*),  intent(in)    :: dir
    type(mpi_grp_t),   intent(in)    :: mpi_grp
    integer,           intent(out)   :: dim, ierr
    integer,           intent(inout) :: nst, kpoints

    character(len=256) :: line
    character(len=12)  :: filename
    character(len=1)   :: char
    integer :: iunit, iunit2, err, i, ist, idim, ik
    FLOAT :: occ, eigenval

    PUSH_SUB(states_look)

    ierr = 0
    iunit  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(iunit < 0) then
      ierr = -1
    POP_SUB(states_look)
return
    end if
    iunit2 = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(iunit2 < 0) then
      call io_close(iunit, grp = mpi_grp)
      ierr = -1
    POP_SUB(states_look)
return
    end if

    ! Skip two lines.
    call iopar_read(mpi_grp, iunit, line, err)
    call iopar_read(mpi_grp, iunit, line, err)
    call iopar_read(mpi_grp, iunit2, line, err)
    call iopar_read(mpi_grp, iunit2, line, err)

    kpoints = 1
    dim = 1
    nst = 1

    do
      call iopar_read(mpi_grp, iunit, line, i)
      read(line, '(a)') char
      if(i.ne.0.or.char=='%') exit
      read(line, *) ik, char, ist, char, idim, char, filename
      if(idim == 2)    dim     = 2
      call iopar_read(mpi_grp, iunit2, line, err)
      read(line, *) occ, char, eigenval
      if(ik > kpoints) kpoints = ik
      if(ist>nst)      nst     = ist
    end do

    call io_close(iunit, grp = mpi_grp)
    call io_close(iunit2, grp = mpi_grp)

    POP_SUB(states_look)
  end subroutine states_look

  ! ---------------------------------------------------------
  ! Allocate free states.
  subroutine states_allocate_free_states(st, gr)
    type(states_t), intent(inout) :: st
    type(grid_t),   intent(in)    :: gr

    integer :: il

    PUSH_SUB(states_allocate_free_states)

    ! FIXME: spin-polarized free states ignored.
    if(gr%ob_grid%open_boundaries) then
      SAFE_ALLOCATE(st%zphi(1:gr%mesh%np_part, 1:st%ob_d%dim, 1:st%ob_nst, 1:st%ob_d%nik))
      do il = 1, NLEADS
        SAFE_ALLOCATE(st%ob_lead(il)%rho(1:gr%ob_grid%lead(il)%mesh%np, 1:st%d%nspin))
      end do
      st%zphi = M_z0
    else
      nullify(st%zphi)
    end if

    POP_SUB(states_allocate_free_states)
  end subroutine states_allocate_free_states


  ! ---------------------------------------------------------
  ! Deallocate free states.
  subroutine states_deallocate_free_states(st, gr)
    type(states_t), intent(inout) :: st
    type(grid_t),   intent(in)    :: gr

    integer :: il

    PUSH_SUB(states_deallocate_free_states)

    if(gr%ob_grid%open_boundaries) then
      SAFE_DEALLOCATE_P(st%zphi)
      do il = 1, NLEADS
        SAFE_DEALLOCATE_P(st%ob_lead(il)%rho)
      end do
    end if

    POP_SUB(states_deallocate_free_states)
  end subroutine states_deallocate_free_states


  ! ---------------------------------------------------------
  !> Reads from the input file the initial occupations, if the
  !! block "Occupations" is present. Otherwise, it makes an initial
  !! guess for the occupations, maybe using the "Smearing"
  !! variable.
  !! The resulting occupations are placed on the st%occ variable. The
  !! boolean st%fixed_occ is also set to .true., if the occupations are
  !! set by the user through the "Occupations" block; false otherwise.
  subroutine states_read_initial_occs(st, excess_charge)
    type(states_t), intent(inout) :: st
    FLOAT,          intent(in)    :: excess_charge

    integer :: ik, ist, ispin, nspin, ncols, el_per_state
    type(block_t) :: blk
    FLOAT :: rr, charge
    logical :: integral_occs

    PUSH_SUB(states_read_initial_occs)
    !%Variable Occupations
    !%Type block
    !%Section States
    !%Description
    !% The occupation numbers of the orbitals can be fixed through the use of this
    !% variable. For example:
    !%
    !% <tt>%Occupations
    !% <br>&nbsp;&nbsp;2.0 | 2.0 | 2.0 | 2.0 | 2.0
    !% <br>%</tt>
    !%
    !% would fix the occupations of the five states to <i>2.0</i>. There can be
    !% at most as many columns as states in the calculation. If there are fewer columns
    !% than states, then the code will assume that the user is indicating the occupations
    !% of the uppermost states, assigning maximum occupation (i.e. 2 for spin-unpolarized
    !% calculations, 1 otherwise) to the lower states. If <tt>SpinComponents == polarized</tt>
    !% this block should contain two lines, one for each spin channel.
    !% This variable is very useful when dealing with highly symmetric small systems
    !% (like an open-shell atom), for it allows us to fix the occupation numbers
    !% of degenerate states in order to help <tt>octopus</tt> to converge. This is to
    !% be used in conjuction with <tt>ExtraStates</tt>. For example, to calculate the
    !% carbon atom, one would do:
    !%
    !% <tt>ExtraStates = 2
    !% <br>%Occupations
    !% <br>&nbsp;&nbsp;2 | 2/3 | 2/3 | 2/3
    !% <br>%</tt>
    !%
    !% If you want the calculation to be spin-polarized (which makes more sense), you could do:
    !%
    !% <tt>ExtraStates = 2
    !% <br>%Occupations
    !% <br>&nbsp;&nbsp; 2/3 | 2/3 | 2/3
    !% <br>&nbsp;&nbsp; 0   |   0 |   0
    !% <br>%</tt>
    !%
    !% Note that in this case the first state is absent, the code will calculate four states
    !% (two because there are four electrons, plus two because <tt>ExtraStates</tt> = 2), and since
    !% it finds only three columns, it will occupy the first state with one electron for each
    !% of the spin options.
    !%
    !% If the sum of occupations is not equal to the total charge set by <tt>ExcessCharge</tt>,
    !% an error message is printed. If <tt>FromScratch = no</tt>, this block is ignored when restart information
    !% is read, and the previous occupations are used.
    !%End

    integral_occs = .true.

    if(st%open_boundaries) then
      st%fixed_occ = .true.
      st%occ  = st%ob_occ
      st%d%kweights = st%ob_d%kweights
      st%qtot = M_ZERO
      do ist = 1, st%nst
        st%qtot = st%qtot + sum(st%occ(ist, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
      end do

    else
      occ_fix: if(parse_block(datasets_check('Occupations'), blk)==0) then
        ! read in occupations
        st%fixed_occ = .true.

        ! Reads the number of columns in the first row. This assumes that all rows
        ! have the same column number; otherwise the code will stop with an error.
        ncols = parse_block_cols(blk, 0)
        if(ncols > st%nst) then
          call input_error("Occupations")
        end if
        ! Now we fill all the "missing" states with the maximum occupation.
        if(st%d%ispin == UNPOLARIZED) then
          el_per_state = M_TWO
        else
          el_per_state = M_ONE
        endif

        do ik = 1, st%d%nik
          do ist = 1, st%nst - ncols
            st%occ(ist, ik) = el_per_state
          end do
        end do
        do ik = 1, st%d%nik
          do ist = st%nst - ncols + 1, st%nst 
            call parse_block_float(blk, ik-1, ist-1-(st%nst-ncols), st%occ(ist, ik))
            integral_occs = integral_occs .and. &
              abs((st%occ(ist, ik) - el_per_state) * st%occ(ist, ik)) .le. M_EPSILON
          end do
        end do
        call parse_block_end(blk)

      else
        st%fixed_occ = .false.
        integral_occs = .false.

        ! first guess for occupation...paramagnetic configuration
        rr = M_ONE
        if(st%d%ispin == UNPOLARIZED) rr = M_TWO

        st%occ  = M_ZERO
        st%qtot = -(st%val_charge + excess_charge)

        nspin = 1
        if(st%d%nspin == 2) nspin = 2

        do ik = 1, st%d%nik, nspin
          charge = M_ZERO
          do ispin = ik, ik + nspin - 1
            do ist = 1, st%nst
              st%occ(ist, ispin) = min(rr, -(st%val_charge + excess_charge) - charge)
              charge = charge + st%occ(ist, ispin)
            end do
          end do
        end do

      end if occ_fix
    end if

    call smear_init(st%smear, st%d%ispin, st%fixed_occ, integral_occs)

    if(.not. smear_is_semiconducting(st%smear) .and. .not. st%smear%method == SMEAR_FIXED_OCC .and. st%nst * 2 .le. st%qtot) then
      message(1) = "Smearing needs unoccupied states (via ExtraStates) to be useful."
      call messages_warning(1)
    endif

    ! sanity check
    charge = M_ZERO
    do ist = 1, st%nst
      charge = charge + sum(st%occ(ist, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
    end do
    if(abs(charge - st%qtot) > CNST(1e-6)) then
      message(1) = "Initial occupations do not integrate to total charge."
      write(message(2), '(6x,f12.6,a,f12.6)') charge, ' != ', st%qtot
      call messages_fatal(2)
    end if

    POP_SUB(states_read_initial_occs)
  end subroutine states_read_initial_occs


  ! ---------------------------------------------------------
  !> Reads, if present, the "InitialSpins" block. This is only
  !! done in spinors mode; otherwise the routine does nothing. The
  !! resulting spins are placed onto the st%spin pointer. The boolean
  !! st%fixed_spins is set to true if (and only if) the InitialSpins
  !! block is present.
  subroutine states_read_initial_spins(st)
    type(states_t), intent(inout) :: st

    integer :: i, j
    type(block_t) :: blk

    PUSH_SUB(states_read_initial_spins)

    st%fixed_spins = .false.
    if(st%d%ispin .ne. SPINORS) then
      POP_SUB(states_read_initial_spins)
      return
    end if

    !%Variable InitialSpins
    !%Type block
    !%Section States
    !%Description
    !% The spin character of the initial random guesses for the spinors can
    !% be fixed by making use of this block. Note that this will not "fix" the
    !% the spins during the calculation (this cannot be done in spinors mode, in
    !% being able to change the spins is why the spinors mode exists in the first
    !% place).
    !%
    !% This block is meaningless and ignored if the run is not in spinors mode
    !% (<tt>SpinComponents = spinors</tt>). 
    !%
    !% The structure of the block is very simple: each column contains the desired
    !% &lt;<i>S_x</i>&gt;, &lt;<i>S_y</i>&gt;, &lt;<i>S_z</i>&gt; for each spinor.
    !% If the calculation is for a periodic system
    !% and there is more than one <i>k</i>-point, the spins of all the <i>k</i>-points are
    !% the same.
    !%
    !% For example, if we have two spinors, and we want one in the <i>Sx</i> "down" state,
    !% and another one in the <i>Sx</i> "up" state:
    !%
    !% <tt>%InitialSpins
    !% <br>&nbsp;&nbsp;  0.5 | 0.0 | 0.0
    !% <br>&nbsp;&nbsp; -0.5 | 0.0 | 0.0
    !% <br>%</tt>
    !%
    !% WARNING: if the calculation is for a system described by pseudopotentials (as
    !% opposed to user-defined potentials or model systems), this option is
    !% meaningless since the random spinors are overwritten by the atomic orbitals.
    !%
    !% There are a couple of physical constraints that have to be fulfilled:
    !%
    !% (A) | &lt;<i>S_i</i>&gt; | &lt;= 1/2
    !%
    !% (B) &lt;<i>S_x</i>&gt;^2 + &lt;<i>S_y</i>&gt;^2 + &lt;<i>S_z</i>&gt;^2 = 1/4
    !%
    !%End
    spin_fix: if(parse_block(datasets_check('InitialSpins'), blk)==0) then
      do i = 1, st%nst
        do j = 1, 3
          call parse_block_float(blk, i-1, j-1, st%spin(j, i, 1))
        end do
        ! This checks (B).
        if( abs(sum(st%spin(1:3, i, 1)**2) - M_FOURTH) > CNST(1.0e-6)) call input_error('InitialSpins')
      end do
      call parse_block_end(blk)
      ! This checks (A). In fact (A) follows from (B), so maybe this is not necessary...
      if(any(abs(st%spin(:, :, :)) > M_HALF)) then
        call input_error('InitialSpins')
      end if
      st%fixed_spins = .true.
      do i = 2, st%d%nik
        st%spin(:, :, i) = st%spin(:, :, 1)     
      end do
    end if spin_fix

    POP_SUB(states_read_initial_spins)
  end subroutine states_read_initial_spins


  ! ---------------------------------------------------------
  ! Allocates the KS wavefunctions defined within a states_t structure.
  subroutine states_allocate_wfns(st, mesh, wfs_type, ob_mesh)
    type(states_t),         intent(inout)   :: st
    type(mesh_t),           intent(in)      :: mesh
    type(type_t), optional, intent(in)      :: wfs_type
    type(mesh_t), optional, intent(in)      :: ob_mesh(:)

    integer :: ip, ik, ist, idim, st1, st2, k1, k2, size, il
    logical :: force

    PUSH_SUB(states_allocate_wfns)

    if(associated(st%dpsi).or.associated(st%zpsi)) then
      message(1) = "Trying to allocate wavefunctions that are already allocated."
      call messages_fatal(1)
    end if
    
    if (present(wfs_type)) then
      ASSERT(wfs_type == TYPE_FLOAT .or. wfs_type == TYPE_CMPLX)
      st%priv%wfs_type = wfs_type
    end if

    !%Variable ForceComplex
    !%Type logical
    !%Default no
    !%Section Execution::Debug
    !%Description
    !% Normally <tt>Octopus</tt> determines automatically the type necessary
    !% for the wavefunctions. When set to yes this variable will
    !% force the use of complex wavefunctions. 
    !%
    !% Warning: This variable is designed for testing and
    !% benchmarking and normal users need not use it.
    !%
    !%End
    call parse_logical(datasets_check('ForceComplex'), .false., force)

    if(force) call states_set_complex(st)

    st1 = st%st_start
    st2 = st%st_end
    k1 = st%d%kpt%start
    k2 = st%d%kpt%end
    size = mesh%np_part

    if (states_are_real(st)) then
      SAFE_ALLOCATE(st%dpsi(1:size, 1:st%d%dim, st1:st2, k1:k2))

      forall(ik=k1:k2, ist=st1:st2, idim=1:st%d%dim, ip=1:size)
        st%dpsi(ip, idim, ist, ik) = M_ZERO
      end forall

    else
      SAFE_ALLOCATE(st%zpsi(1:size, 1:st%d%dim, st1:st2, k1:k2))

      forall(ik=k1:k2, ist=st1:st2, idim=1:st%d%dim, ip=1:size)
        st%zpsi(ip, idim, ist, ik) = M_Z0
      end forall
    end if

    if(present(ob_mesh)) then
      ASSERT(st%open_boundaries)

      do il = 1, NLEADS
        SAFE_ALLOCATE(st%ob_lead(il)%intf_psi(1:ob_mesh(il)%np, 1:st%d%dim, st1:st2, k1:k2))
        st%ob_lead(il)%intf_psi = M_z0
      end do
    end if

    call states_init_block(st)

    POP_SUB(states_allocate_wfns)
  end subroutine states_allocate_wfns
  
  ! -----------------------------------------------------
  
  subroutine states_init_block(st)
    type(states_t),    intent(inout)   :: st

    integer :: ib, iq, ist
    logical :: same_node
    integer, allocatable :: bstart(:), bend(:)

    PUSH_SUB(states_init_block)

    SAFE_ALLOCATE(bstart(1:st%nst))
    SAFE_ALLOCATE(bend(1:st%nst))
    SAFE_ALLOCATE(st%iblock(1:st%nst, 1:st%d%nik))
    st%iblock = 0

    ! count and assign blocks
    ib = 0
    st%nblocks = 0
    bstart(1) = 1
    do ist = 1, st%nst
      INCR(ib, 1)

      st%iblock(ist, st%d%kpt%start:st%d%kpt%end) = st%nblocks + 1

      same_node = .true.
      if(st%parallel_in_states .and. ist /= st%nst) then
        ! We have to avoid that states that are in different nodes end
        ! up in the same block
        same_node = (st%node(ist + 1) == st%node(ist))
      end if

      if(ib == st%d%block_size .or. ist == st%nst .or. .not. same_node) then
        ib = 0
        INCR(st%nblocks, 1)
        bend(st%nblocks) = ist
        if(ist /= st%nst) bstart(st%nblocks + 1) = ist + 1
      end if
    end do

!!$    ! some debug output that I will keep here for the moment
!!$    if(mpi_grp_is_root(mpi_world)) then
!!$      print*, "NST       ", st%nst
!!$      print*, "BLOCKSIZE ", st%d%block_size
!!$      print*, "NBLOCKS   ", st%nblocks
!!$        
!!$      print*, "===============" 
!!$      do ist = 1, st%nst
!!$        print*, st%node(ist), ist, st%iblock(ist, 1)
!!$      end do
!!$      print*, "===============" 
!!$      
!!$      do ib = 1, st%nblocks
!!$        print*, ib, bstart(ib), bend(ib)
!!$      end do
!!$      
!!$    end if
    
    SAFE_ALLOCATE(st%psib(1:st%nblocks, 1:st%d%nik))
    SAFE_ALLOCATE(st%block_is_local(1:st%nblocks, 1:st%d%nik))
    st%block_is_local = .false.

    do ib = 1, st%nblocks
      if(bstart(ib) >= st%st_start .and. bend(ib) <= st%st_end) then
        do iq = st%d%kpt%start, st%d%kpt%end
          st%block_is_local(ib, iq) = .true.

          if (states_are_real(st)) then
            ASSERT(associated(st%dpsi))
            call batch_init(st%psib(ib, iq), st%d%dim, bstart(ib), bend(ib), st%dpsi(:, :, bstart(ib):bend(ib), iq))
          else
            ASSERT(associated(st%zpsi))
            call batch_init(st%psib(ib, iq), st%d%dim, bstart(ib), bend(ib), st%zpsi(:, :, bstart(ib):bend(ib), iq))
          end if
          
        end do
      end if
    end do

    st%block_initialized = .true.

    SAFE_DEALLOCATE_A(bstart)
    SAFE_DEALLOCATE_A(bend)
    POP_SUB(states_init_block)
  end subroutine states_init_block


  ! ---------------------------------------------------------
  !> Deallocates the KS wavefunctions defined within a states_t structure.
  subroutine states_deallocate_wfns(st)
    type(states_t), intent(inout) :: st

    integer :: il, ib, iq

    PUSH_SUB(states_deallocate_wfns)

    if (states_are_real(st)) then
      SAFE_DEALLOCATE_P(st%dpsi)
    else
      SAFE_DEALLOCATE_P(st%zpsi)
    end if

    if (st%block_initialized) then
       do ib = 1, st%nblocks
          do iq = st%d%kpt%start, st%d%kpt%end
             if(st%block_is_local(ib, iq)) call batch_end(st%psib(ib, iq))
          end do
       end do

       SAFE_DEALLOCATE_P(st%psib)
       SAFE_DEALLOCATE_P(st%block_is_local)
       SAFE_DEALLOCATE_P(st%iblock)
       st%block_initialized = .false.
    end if

    if(st%open_boundaries) then
      do il = 1, NLEADS
        SAFE_DEALLOCATE_P(st%ob_lead(il)%intf_psi)
      end do
    end if

    POP_SUB(states_deallocate_wfns)
  end subroutine states_deallocate_wfns


  ! ---------------------------------------------------------
  subroutine states_densities_init(st, gr, geo, mc)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo
    type(multicomm_t), intent(in)    :: mc

    PUSH_SUB(states_densities_init)

    ! allocate arrays for charge and current densities
    SAFE_ALLOCATE(st%rho(1:gr%fine%mesh%np_part, 1:st%d%nspin))
    st%rho  = M_ZERO
    if(st%d%cdft) then
      SAFE_ALLOCATE(st%current(1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:st%d%nspin))
      st%current = M_ZERO
    end if
    st%nlcc = geo%nlcc
    if(st%nlcc) then
      SAFE_ALLOCATE(st%rho_core(1:gr%fine%mesh%np))
      st%rho_core(:) = M_ZERO
    end if

    ! This has to be here as it requires mc%nthreads that is not available in states_init
    call states_exec_init()

    POP_SUB(states_densities_init)
    
  contains

    subroutine states_exec_init()
      integer :: default

      PUSH_SUB(states_densities_init.states_exec_init)

      !%Variable StatesBlockSize
      !%Type integer
      !%Default max(4, 2*nthreads)
      !%Section Execution::Optimization
      !%Description
      !% Some routines work over blocks of eigenfunctions, which
      !% generally improves performance at the expense of increased
      !% memory consumption. This variable selects the size of the
      !% blocks to be used.
      !%End

      call parse_integer(datasets_check('StatesBlockSize'), max(4, 2*mc%nthreads), st%d%block_size)
      if(st%d%block_size < 1) then
        message(1) = "Error: The variable 'StatesBlockSize' must be greater than 0."
        call messages_fatal(1)
      end if

      st%d%block_size = min(st%d%block_size, st%nst)

      !%Variable StatesOrthogonalization
      !%Type integer
      !%Default gram_schmidt
      !%Section Execution::Optimization
      !%Description
      !% The full orthogonalization method used by some
      !% eigensolvers. The default is gram_schmidt. When state
      !% parallelization the default is par_gram_schmidt.
      !%Option gram_schmidt 1
      !% The standard Gram-Schmidt orthogonalization implemented using
      !% Blas/Lapack. Can be used with domain parallelization but not
      !% state parallelization.
      !%Option par_gram_schmidt 2
      !% The standard Gram-Schmidt orthogonalization implemented using
      !% Scalapack. Compatible with states parallelization.
      !%Option mgs 3
      !% Modified Gram-Schmidt orthogonalization.
      !%Option qr 4
      !% (Experimental) Orthogonalization is performed based on a QR
      !% decomposition based on Lapack routines _getrf and _orgqr.
      !% Compatible with states parallelization. 
      !%Option old_gram_schmidt 5
      !% Old Gram-Schmidt implementation, compatible with states
      !% parallelization.
      !%End

      if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then
#ifdef HAVE_SCALAPACK
        default = ORTH_PAR_GS
#else        
        default = ORTH_OLDGS
#endif
      else
        default = ORTH_GS
      end if
      
      call parse_integer(datasets_check('StatesOrthogonalization'), default, st%d%orth_method)

      if(.not.varinfo_valid_option('StatesOrthogonalization', st%d%orth_method)) call input_error('StatesOrthogonalization')
      call messages_print_var_option(stdout, 'StatesOrthogonalization', st%d%orth_method)


      if(st%d%orth_method == ORTH_QR) call messages_experimental("QR Orthogonalization")

#ifdef HAVE_MPI
      if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then
        select case(st%d%orth_method)
        case(ORTH_OLDGS, ORTH_PAR_GS)
        case(ORTH_QR)
#ifndef HAVE_SCALAPACK
          message(1) = 'The QR orthogonalizer requires scalapack to work with state-parallelization.'
          call messages_fatal(1)
#endif
        case default
          message(1) = 'The selected orthogonalization method cannot work with state-parallelization.'
          call messages_fatal(1)
        end select
      end if
#endif

      ! some checks for ingenious users
      if(st%d%orth_method == ORTH_PAR_GS) then
#ifndef HAVE_MPI
        message(1) = 'The parallel gram schmidt orthogonalizer can only be used in parallel.'
        call messages_fatal(1)
#else
#ifndef HAVE_SCALAPACK
        message(1) = 'The parallel gram schmidt orthogonalizer requires scalapack.'
        call messages_fatal(1)
#endif
        if(st%dom_st_mpi_grp%size == 1) then
          message(1) = 'The parallel gram schmidt orthogonalizer is designed to be used with domain or state parallelization.'
          call messages_warning(1)
        end if
#endif
      end if

      POP_SUB(states_densities_init.states_exec_init)
    end subroutine states_exec_init
  end subroutine states_densities_init


  ! ---------------------------------------------------------
  subroutine states_copy(stout, stin)
    type(states_t), intent(inout) :: stout
    type(states_t), intent(in)    :: stin

    PUSH_SUB(states_copy)

    call states_null(stout)

    stout%priv%wfs_type   = stin%priv%wfs_type
    call states_dim_copy(stout%d, stin%d)
    stout%nst        = stin%nst
    stout%qtot       = stin%qtot
    stout%val_charge = stin%val_charge
    call smear_copy(stout%smear, stin%smear)
    stout%parallel_in_states = stin%parallel_in_states
    stout%lnst       = stin%lnst
    stout%st_start   = stin%st_start
    stout%st_end     = stin%st_end
    call loct_pointer_copy(stout%dpsi, stin%dpsi)
    call loct_pointer_copy(stout%zpsi, stin%zpsi)
    call loct_pointer_copy(stout%user_def_states, stin%user_def_states)

    call loct_pointer_copy(stout%rho, stin%rho)
    call loct_pointer_copy(stout%current, stin%current)
    stout%nlcc = stin%nlcc
    call loct_pointer_copy(stout%rho_core, stin%rho_core)
    call loct_pointer_copy(stout%frozen_rho, stin%frozen_rho)
    call loct_pointer_copy(stout%eigenval, stin%eigenval)
    stout%fixed_occ = stin%fixed_occ
    call loct_pointer_copy(stout%occ, stin%occ)
    stout%fixed_spins = stin%fixed_spins
    call loct_pointer_copy(stout%spin, stin%spin)
    call loct_pointer_copy(stout%node, stin%node)
    call mpi_grp_copy(stout%mpi_grp, stin%mpi_grp)
    call loct_pointer_copy(stout%st_range, stin%st_range)
    call loct_pointer_copy(stout%st_num, stin%st_num)

    call modelmb_particles_copy(stout%modelmbparticles, stin%modelmbparticles)

    if(stin%parallel_in_states) call multicomm_all_pairs_copy(stout%ap, stin%ap)

    stout%open_boundaries = stin%open_boundaries
    ! Some of the "open boundaries" variables are not copied.

    stout%symmetrize_density = stin%symmetrize_density

    if(associated(stout%dpsi) .or. associated(stout%zpsi)) then
      call states_init_block(stout)
    end if

    stout%dom_st_kpt_mpi_grp = stin%dom_st_kpt_mpi_grp
    stout%st_kpt_mpi_grp = stin%st_kpt_mpi_grp
    stout%st_kpt_mpi_grp = stin%st_kpt_mpi_grp

#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_copy(stin%dom_st_proc_grid, stout%dom_st_proc_grid)
#endif

    POP_SUB(states_copy)
  end subroutine states_copy


  ! ---------------------------------------------------------
  subroutine states_end(st)
    type(states_t), intent(inout) :: st

    integer :: il

    PUSH_SUB(states_end)
    
#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_end(st%dom_st_proc_grid)
#endif

    call states_deallocate_wfns(st)

    SAFE_DEALLOCATE_P(st%user_def_states)

    SAFE_DEALLOCATE_P(st%rho)
    SAFE_DEALLOCATE_P(st%current)
    SAFE_DEALLOCATE_P(st%rho_core)
    SAFE_DEALLOCATE_P(st%frozen_rho)
    SAFE_DEALLOCATE_P(st%eigenval)

    SAFE_DEALLOCATE_P(st%occ)
    SAFE_DEALLOCATE_P(st%spin)

    SAFE_DEALLOCATE_P(st%node)
    SAFE_DEALLOCATE_P(st%st_range)
    SAFE_DEALLOCATE_P(st%st_num)

    if(st%parallel_in_states) then
      SAFE_DEALLOCATE_P(st%ap%schedule)
    end if

    SAFE_DEALLOCATE_P(st%zphi)
    call states_dim_end(st%d)

    call states_dim_end(st%ob_d)
    SAFE_DEALLOCATE_P(st%ob_eigenval)
    SAFE_DEALLOCATE_P(st%ob_occ)
    do il = 1, NLEADS
      SAFE_DEALLOCATE_P(st%ob_lead(il)%self_energy)
    end do

    call modelmb_particles_end(st%modelmbparticles)

    POP_SUB(states_end)
  end subroutine states_end

  ! ---------------------------------------------------------
  !> generate a hydrogen s-wavefunction around a random point
  subroutine states_generate_random(st, mesh, ist_start_, ist_end_)
    type(states_t),    intent(inout) :: st
    type(mesh_t),      intent(in)    :: mesh
    integer, optional, intent(in)    :: ist_start_, ist_end_

    integer :: ist, ik, id, ist_start, ist_end, jst, seed
    CMPLX   :: alpha, beta

    PUSH_SUB(states_generate_random)

    ist_start = st%st_start
    if(present(ist_start_)) ist_start = max(ist_start, ist_start_)
    ist_end = st%st_end
    if(present(ist_end_)) ist_end = min(ist_end, ist_end_)

    if(st%parallel_in_states) then
      seed = st%mpi_grp%rank
    else
      seed = 0
    end if
    if(st%d%kpt%parallel) then
      seed = st%d%kpt%mpi_grp%rank
    else
      seed = 0
    end if

    select case(st%d%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = ist_start, ist_end
          if (states_are_real(st)) then
            call dmf_random(mesh, st%dpsi(:, 1, ist, ik), seed)
          else
            call zmf_random(mesh, st%zpsi(:, 1, ist, ik), seed)
          end if
          st%eigenval(ist, ik) = M_ZERO
        end do
      end do

    case(SPINORS)

      ASSERT(states_are_complex(st))

      if(st%fixed_spins) then
        
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = ist_start, ist_end
            call zmf_random(mesh, st%zpsi(:, 1, ist, ik))
            ! In this case, the spinors are made of a spatial part times a vector [alpha beta]^T in 
            ! spin space (i.e., same spatial part for each spin component). So (alpha, beta)
            ! determines the spin values. The values of (alpha, beta) can be be obtained
            ! with simple formulae from <Sx>, <Sy>, <Sz>.
            !
            ! Note that here we orthonormalize the orbital part. This ensures that the spinors
            ! are untouched later in the general orthonormalization, and therefore the spin values
            ! of each spinor remain the same.
            do jst = ist_start, ist - 1
              st%zpsi(:, 1, ist, ik) = st%zpsi(:, 1, ist, ik) - &
                                       zmf_dotp(mesh, st%zpsi(:, 1, ist, ik), st%zpsi(:, 1, jst, ik)) * &
                                       st%zpsi(:, 1, jst, ik)
            end do
            st%zpsi(:, 1, ist, ik) = st%zpsi(:, 1, ist, ik) / zmf_nrm2(mesh, st%zpsi(:, 1, ist, ik))
            st%zpsi(:, 2, ist, ik) = st%zpsi(:, 1, ist, ik)
            alpha = cmplx(sqrt(M_HALF + st%spin(3, ist, ik)), M_ZERO, REAL_PRECISION)
            beta  = cmplx(sqrt(M_ONE - abs(alpha)**2), M_ZERO, REAL_PRECISION)
            if(abs(alpha) > M_ZERO) then
              beta = cmplx(st%spin(1, ist, ik) / abs(alpha), st%spin(2, ist, ik) / abs(alpha), REAL_PRECISION)
            end if
            st%zpsi(:, 1, ist, ik) = alpha * st%zpsi(:, 1, ist, ik)
            st%zpsi(:, 2, ist, ik) = beta * st%zpsi(:, 2, ist, ik)
            st%eigenval(ist, ik) = M_ZERO
          end do
        end do
      else
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = ist_start, ist_end
            do id = 1, st%d%dim
              call zmf_random(mesh, st%zpsi(:, id, ist, ik))
            end do
            st%eigenval(ist, ik) = M_HUGE
          end do
        end do
      end if

    end select

    POP_SUB(states_generate_random)
  end subroutine states_generate_random

  ! ---------------------------------------------------------
  subroutine states_fermi(st, mesh)
    type(states_t), intent(inout) :: st
    type(mesh_t),   intent(in)    :: mesh

    ! Local variables.
    integer            :: ist, ik
    FLOAT              :: charge
#if defined(HAVE_MPI)
    integer            :: jj
    integer            :: tmp
    FLOAT, allocatable :: lspin(:, :) ! To exchange spin.
#endif

    PUSH_SUB(states_fermi)

    call smear_find_fermi_energy(st%smear, st%eigenval, st%occ, st%qtot, &
      st%d%nik, st%nst, st%d%kweights)

    call smear_fill_occupations(st%smear, st%eigenval, st%occ, &
      st%d%nik, st%nst)

    ! check if everything is OK
    charge = M_ZERO
    do ist = 1, st%nst
      charge = charge + sum(st%occ(ist, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
    end do
    if(abs(charge-st%qtot) > CNST(1e-6)) then
      message(1) = 'Occupations do not integrate to total charge.'
      write(message(2), '(6x,f12.8,a,f12.8)') charge, ' != ', st%qtot
      call messages_warning(2)
    end if

    if(st%d%ispin == SPINORS) then
      ASSERT(states_are_complex(st))

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          st%spin(1:3, ist, ik) = state_spin(mesh, st%zpsi(:, :, ist, ik))
        end do
#if defined(HAVE_MPI)
        if(st%parallel_in_states) then
          SAFE_ALLOCATE(lspin(1:3, 1:st%lnst))
          lspin = st%spin(1:3, st%st_start:st%st_end, ik)
          do jj = 1, 3
            call lmpi_gen_allgatherv(st%lnst, lspin(jj, :), tmp, st%spin(jj, :, ik), st%mpi_grp)
          end do
          SAFE_DEALLOCATE_A(lspin)
        end if
#endif
      end do
    end if

    POP_SUB(states_fermi)
  end subroutine states_fermi


  ! ---------------------------------------------------------
  !> function to calculate the eigenvalues sum using occupations as weights
  function states_eigenvalues_sum(st, alt_eig) result(tot)
    type(states_t), intent(in)  :: st
    FLOAT, optional, intent(in) :: alt_eig(st%st_start:st%st_end, 1:st%d%nik)
    FLOAT                       :: tot

    integer :: ik

    PUSH_SUB(states_eigenvalues_sum)

    tot = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      if(present(alt_eig)) then
        tot = tot + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik) * &
          alt_eig(st%st_start:st%st_end, ik))
      else
        tot = tot + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik) * &
          st%eigenval(st%st_start:st%st_end, ik))
      end if
    end do

    if(st%parallel_in_states .or. st%d%kpt%parallel) call comm_allreduce(st%st_kpt_mpi_grp%comm, tot)

    POP_SUB(states_eigenvalues_sum)
  end function states_eigenvalues_sum

  ! -------------------------------------------------------
  integer pure function states_spin_channel(ispin, ik, dim)
    integer, intent(in) :: ispin, ik, dim

    select case(ispin)
    case(1); states_spin_channel = 1
    case(2); states_spin_channel = mod(ik+1, 2)+1
    case(3); states_spin_channel = dim
    case default; states_spin_channel = -1
    end select

  end function states_spin_channel


  ! ---------------------------------------------------------
  subroutine states_distribute_nodes(st, mc)
    type(states_t),    intent(inout) :: st
    type(multicomm_t), intent(in)    :: mc

#ifdef HAVE_MPI
    integer :: inode, ist
#endif

    PUSH_SUB(states_distribute_nodes)

    ! Defaults.
    st%node(:)            = 0
    st%st_start           = 1
    st%st_end             = st%nst
    st%lnst               = st%nst
    st%parallel_in_states = .false.
    call mpi_grp_init(st%mpi_grp, mc%group_comm(P_STRATEGY_STATES))
    call mpi_grp_init(st%dom_st_kpt_mpi_grp, mc%dom_st_kpt_comm)
    call mpi_grp_init(st%dom_st_mpi_grp, mc%dom_st_comm)
    call mpi_grp_init(st%st_kpt_mpi_grp, mc%st_kpt_comm)

#ifdef HAVE_SCALAPACK
    if(calc_mode_scalapack_compat() .and. .not. st%d%kpt%parallel) then
      call blacs_proc_grid_init(st%dom_st_proc_grid, st%dom_st_mpi_grp)
    end if
#endif

#if defined(HAVE_MPI)
    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then
      st%parallel_in_states = .true.
      
      call multicomm_create_all_pairs(st%mpi_grp, st%ap)

     if(st%nst < st%mpi_grp%size) then
       message(1) = "Have more processors than necessary"
       write(message(2),'(i4,a,i4,a)') st%mpi_grp%size, " processors and ", st%nst, " states."
       call messages_fatal(2)
     end if

     SAFE_ALLOCATE(st%st_range(1:2, 0:st%mpi_grp%size-1))
     SAFE_ALLOCATE(st%st_num(0:st%mpi_grp%size-1))

     call multicomm_divide_range(st%nst, st%mpi_grp%size, st%st_range(1, :), st%st_range(2, :), &
       lsize = st%st_num, scalapack_compat = calc_mode_scalapack_compat())

     message(1) = "Info: Parallelization in states"
     call messages_info(1)

     do inode = 0, st%mpi_grp%size - 1
       write(message(1),'(a,i4,a,i5,a,i6,a,i6)') &
            'Info: Nodes in states-group ', inode, ' will manage ', st%st_num(inode), ' states:', &
            st%st_range(1, inode), " - ", st%st_range(2, inode)
       call messages_info(1)

       do ist = st%st_range(1, inode), st%st_range(2, inode)
         st%node(ist) = inode
       end do
     end do

     st%st_start = st%st_range(1, st%mpi_grp%rank)
     st%st_end   = st%st_range(2, st%mpi_grp%rank)
     st%lnst     = st%st_num(st%mpi_grp%rank)

   end if
#endif

    POP_SUB(states_distribute_nodes)
  end subroutine states_distribute_nodes


  ! ---------------------------------------------------------
  subroutine states_set_complex(st)
    type(states_t),    intent(inout) :: st

    PUSH_SUB(states_set_complex)
    st%priv%wfs_type = TYPE_CMPLX

    POP_SUB(states_set_complex)
  end subroutine states_set_complex

  ! ---------------------------------------------------------
  pure logical function states_are_complex(st) result (wac)
    type(states_t),    intent(in) :: st

    wac = (st%priv%wfs_type == TYPE_CMPLX)

  end function states_are_complex


  ! ---------------------------------------------------------
  pure logical function states_are_real(st) result (war)
    type(states_t),    intent(in) :: st

    war = (st%priv%wfs_type == TYPE_FLOAT)

  end function states_are_real

  ! --------------------------------------------------------- 
  !
  !> This function can calculate several quantities that depend on
  !! derivatives of the orbitals from the states and the density.
  !! The quantities to be calculated depend on the arguments passed.
  subroutine states_calc_quantities(der, st, &
    kinetic_energy_density, paramagnetic_current, density_gradient, density_laplacian, gi_kinetic_energy_density)
    type(derivatives_t),     intent(in)    :: der
    type(states_t),          intent(inout) :: st
    FLOAT, optional, target, intent(out)   :: kinetic_energy_density(:,:)       !< The kinetic energy density.
    FLOAT, optional, target, intent(out)   :: paramagnetic_current(:,:,:)       !< The paramagnetic current.
    FLOAT, optional,         intent(out)   :: density_gradient(:,:,:)           !< The gradient of the density.
    FLOAT, optional,         intent(out)   :: density_laplacian(:,:)            !< The Laplacian of the density.
    FLOAT, optional,         intent(out)   :: gi_kinetic_energy_density(:,:)    !< The gauge-invariant kinetic energy density.

    FLOAT, pointer :: jp(:, :, :)
    FLOAT, pointer :: tau(:, :)
    CMPLX, allocatable :: wf_psi(:,:), gwf_psi(:,:,:), lwf_psi(:,:)
    CMPLX   :: c_tmp
    integer :: sp, is, ik, ik_tmp, ist, i_dim, st_dim, ii
    FLOAT   :: ww, kpoint(1:MAX_DIM)
    logical :: something_to_do

    PUSH_SUB(states_calc_quantities)

    something_to_do = present(kinetic_energy_density) .or. present(paramagnetic_current) .or. &
      present(density_gradient) .or. present(density_laplacian)
    ASSERT(something_to_do)

    SAFE_ALLOCATE( wf_psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gwf_psi(1:der%mesh%np, 1:der%mesh%sb%dim, 1:st%d%dim))
    if(present(density_laplacian)) SAFE_ALLOCATE(lwf_psi(1:der%mesh%np, 1:st%d%dim))

    sp = 1
    if(st%d%ispin == SPIN_POLARIZED) sp = 2


    nullify(tau)
    if(present(kinetic_energy_density)) tau => kinetic_energy_density

    nullify(jp)
    if(present(paramagnetic_current)) jp => paramagnetic_current

    ! for the gauge-invariant kinetic energy density we need the
    ! current and the kinetic energy density
    if(present(gi_kinetic_energy_density)) then 
      if(.not. present(paramagnetic_current) .and. states_are_complex(st)) then
        SAFE_ALLOCATE(jp(1:der%mesh%np, 1:der%mesh%sb%dim, 1:st%d%nspin))
      end if
      if(.not. present(kinetic_energy_density)) then
        SAFE_ALLOCATE(tau(1:der%mesh%np, 1:st%d%nspin))
      end if
    end if

    if(associated(tau)) tau = M_ZERO
    if(associated(jp)) jp = M_ZERO
    if(present(density_gradient)) density_gradient(:,:,:) = M_ZERO
    if(present(density_laplacian)) density_laplacian(:,:) = M_ZERO
    if(present(gi_kinetic_energy_density)) gi_kinetic_energy_density = M_ZERO

    do is = 1, sp
      do ik_tmp = st%d%kpt%start, st%d%kpt%end, sp
        ik = ik_tmp + is - 1

        kpoint(1:der%mesh%sb%dim) = kpoints_get_point(der%mesh%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))

        do ist = st%st_start, st%st_end

          ! all calculations will be done with complex wavefunctions
          if (states_are_real(st)) then
            wf_psi(:,:) = cmplx(st%dpsi(:,:, ist, ik), KIND=REAL_PRECISION)
          else
            wf_psi(:,:) = st%zpsi(:,:, ist, ik)
          end if

          ! calculate gradient of the wavefunction
          do st_dim = 1, st%d%dim
            call zderivatives_grad(der, wf_psi(:,st_dim), gwf_psi(:,:,st_dim))
          end do

          ! calculate the Laplacian of the wavefunction
          if (present(density_laplacian)) then
            do st_dim = 1, st%d%dim
              call zderivatives_lapl(der, wf_psi(:,st_dim), lwf_psi(:,st_dim))
            end do
          end if

          ww = st%d%kweights(ik)*st%occ(ist, ik)

          if(present(density_laplacian)) then
            density_laplacian(1:der%mesh%np, is) = density_laplacian(1:der%mesh%np, is) + &
              ww*M_TWO*real(conjg(wf_psi(1:der%mesh%np, 1))*lwf_psi(1:der%mesh%np, 1))
            if(st%d%ispin == SPINORS) then
              density_laplacian(1:der%mesh%np, 2) = density_laplacian(1:der%mesh%np, 2) + &
                ww*M_TWO*real(conjg(wf_psi(1:der%mesh%np, 2))*lwf_psi(1:der%mesh%np, 2))
              density_laplacian(1:der%mesh%np, 3) = density_laplacian(1:der%mesh%np, 3) + &
                ww*real (lwf_psi(1:der%mesh%np, 1)*conjg(wf_psi(1:der%mesh%np, 2)) + &
                wf_psi(1:der%mesh%np, 1)*conjg(lwf_psi(1:der%mesh%np, 2)))
              density_laplacian(1:der%mesh%np, 4) = density_laplacian(1:der%mesh%np, 4) + &
                ww*aimag(lwf_psi(1:der%mesh%np, 1)*conjg(wf_psi(1:der%mesh%np, 2)) + &
                wf_psi(1:der%mesh%np, 1)*conjg(lwf_psi(1:der%mesh%np, 2)))
            end if
          end if

          do i_dim = 1, der%mesh%sb%dim
            if(present(density_gradient)) &
              density_gradient(1:der%mesh%np, i_dim, is) = density_gradient(1:der%mesh%np, i_dim, is) + &
              ww*M_TWO*real(conjg(wf_psi(1:der%mesh%np, 1))*gwf_psi(1:der%mesh%np, i_dim, 1))
            if(present(density_laplacian)) &
              density_laplacian(1:der%mesh%np, is) = density_laplacian(1:der%mesh%np, is)         + &
              ww*M_TWO*real(conjg(gwf_psi(1:der%mesh%np, i_dim, 1))*gwf_psi(1:der%mesh%np, i_dim, 1))

            if(associated(jp)) then
              if (.not.(states_are_real(st))) then
                jp(1:der%mesh%np, i_dim, is) = jp(1:der%mesh%np, i_dim, is) + &
                  ww*aimag(conjg(wf_psi(1:der%mesh%np, 1))*gwf_psi(1:der%mesh%np, i_dim, 1) - &
                  M_zI*(wf_psi(1:der%mesh%np, 1))**2*kpoint(i_dim ) )
              else
                jp(1:der%mesh%np, i_dim, is) = M_ZERO
              end if
            end if

            if (associated(tau)) then 
              tau (1:der%mesh%np, is)   = tau (1:der%mesh%np, is)        + &
                ww*abs(gwf_psi(1:der%mesh%np, i_dim, 1))**2  &
                + ww*abs(kpoint(i_dim))**2*abs(wf_psi(1:der%mesh%np, 1))**2  &
                - ww*M_TWO*aimag(conjg(wf_psi(1:der%mesh%np, 1))*kpoint(i_dim)*gwf_psi(1:der%mesh%np, i_dim, 1) )
            end if

            if(present(gi_kinetic_energy_density)) then
              ASSERT(associated(tau))
              if(states_are_complex(st) .and. st%current_in_tau) then
                ASSERT(associated(jp))
                gi_kinetic_energy_density(1:der%mesh%np, is) = tau(1:der%mesh%np, is) - &
                  jp(1:der%mesh%np, i_dim, 1)**2/st%rho(1:der%mesh%np, 1)
              else
                gi_kinetic_energy_density(1:der%mesh%np, is) = tau(1:der%mesh%np, is)
              end if
            end if

            if(st%d%ispin == SPINORS) then
              if(present(density_gradient)) then
                density_gradient(1:der%mesh%np, i_dim, 2) = density_gradient(1:der%mesh%np, i_dim, 2) + &
                  ww*M_TWO*real(conjg(wf_psi(1:der%mesh%np, 2))*gwf_psi(1:der%mesh%np, i_dim, 2))
                density_gradient(1:der%mesh%np, i_dim, 3) = density_gradient(1:der%mesh%np, i_dim, 3) + ww* &
                  real (gwf_psi(1:der%mesh%np, i_dim, 1)*conjg(wf_psi(1:der%mesh%np, 2)) + &
                  wf_psi(1:der%mesh%np, 1)*conjg(gwf_psi(1:der%mesh%np, i_dim, 2)))
                density_gradient(1:der%mesh%np, i_dim, 4) = density_gradient(1:der%mesh%np, i_dim, 4) + ww* &
                  aimag(gwf_psi(1:der%mesh%np, i_dim, 1)*conjg(wf_psi(1:der%mesh%np, 2)) + &
                  wf_psi(1:der%mesh%np, 1)*conjg(gwf_psi(1:der%mesh%np, i_dim, 2)))
              end if

              if(present(density_laplacian)) then
                density_laplacian(1:der%mesh%np, 2) = density_laplacian(1:der%mesh%np, 2)         + &
                  ww*M_TWO*real(conjg(gwf_psi(1:der%mesh%np, i_dim, 2))*gwf_psi(1:der%mesh%np, i_dim, 2))
                density_laplacian(1:der%mesh%np, 3) = density_laplacian(1:der%mesh%np, 3)         + &
                  ww*M_TWO*real (gwf_psi(1:der%mesh%np, i_dim, 1)*conjg(gwf_psi(1:der%mesh%np, i_dim, 2)))
                density_laplacian(1:der%mesh%np, 4) = density_laplacian(1:der%mesh%np, 4)         + &
                  ww*M_TWO*aimag(gwf_psi(1:der%mesh%np, i_dim, 1)*conjg(gwf_psi(1:der%mesh%np, i_dim, 2)))
              end if

              ! the expression for the paramagnetic current with spinors is
              !     j = ( jp(1)             jp(3) + i jp(4) ) 
              !         (-jp(3) + i jp(4)   jp(2)           )
              if(associated(jp)) then
                jp(1:der%mesh%np, i_dim, 2) = jp(1:der%mesh%np, i_dim, 2) + &
                  ww*aimag(conjg(wf_psi(1:der%mesh%np, 2))*gwf_psi(1:der%mesh%np, i_dim, 2))
                do ii = 1, der%mesh%np
                  c_tmp = conjg(wf_psi(ii, 1))*gwf_psi(ii, i_dim, 2) - wf_psi(ii, 2)*conjg(gwf_psi(ii, i_dim, 1))
                  jp(ii, i_dim, 3) = jp(ii, i_dim, 3) + ww* real(c_tmp)
                  jp(ii, i_dim, 4) = jp(ii, i_dim, 4) + ww*aimag(c_tmp)
                end do
              end if

              ! the expression for the paramagnetic current with spinors is
              !     t = ( tau(1)              tau(3) + i tau(4) ) 
              !         ( tau(3) - i tau(4)   tau(2)            )
              if(associated(tau)) then
                tau (1:der%mesh%np, 2) = tau (1:der%mesh%np, 2) + ww*abs(gwf_psi(1:der%mesh%np, i_dim, 2))**2
                do ii = 1, der%mesh%np
                  c_tmp = conjg(gwf_psi(ii, i_dim, 1))*gwf_psi(ii, i_dim, 2)
                  tau(ii, 3) = tau(ii, 3) + ww* real(c_tmp)
                  tau(ii, 4) = tau(ii, 4) + ww*aimag(c_tmp)
                end do
              end if

              ASSERT(.not. present(gi_kinetic_energy_density))
              
            end if !SPINORS

          end do

        end do
      end do
    end do

    SAFE_DEALLOCATE_A(wf_psi)
    SAFE_DEALLOCATE_A(gwf_psi)
    SAFE_DEALLOCATE_A(lwf_psi)

    if(.not. present(paramagnetic_current)) then
      SAFE_DEALLOCATE_P(jp)
    end if

    if(.not. present(kinetic_energy_density)) then
      SAFE_DEALLOCATE_P(tau)
    end if

    if(st%parallel_in_states .or. st%d%kpt%parallel) call reduce_all(st%st_kpt_mpi_grp)

    POP_SUB(states_calc_quantities)

  contains 

    subroutine reduce_all(grp)
      type(mpi_grp_t), intent(in)  :: grp

      PUSH_SUB(states_calc_quantities.reduce_all)

      if(associated(tau)) call comm_allreduce(grp%comm, tau, dim = (/der%mesh%np, st%d%nspin/))

      if(present(gi_kinetic_energy_density)) &
        call comm_allreduce(grp%comm, gi_kinetic_energy_density, dim = (/der%mesh%np, st%d%nspin/))

      if (present(density_laplacian)) call comm_allreduce(grp%comm, density_laplacian, dim = (/der%mesh%np, st%d%nspin/))
      
      do is = 1, st%d%nspin
        if(associated(jp)) call comm_allreduce(grp%comm, jp(:, :, is), dim = (/der%mesh%np, der%mesh%sb%dim/))

        if(present(density_gradient)) &
          call comm_allreduce(grp%comm, density_gradient(:, :, is), dim = (/der%mesh%np, der%mesh%sb%dim/))
      end do
      
      POP_SUB(states_calc_quantities.reduce_all)
    end subroutine reduce_all

  end subroutine states_calc_quantities


  ! ---------------------------------------------------------
  function state_spin(mesh, f1) result(spin)
    type(mesh_t), intent(in) :: mesh
    CMPLX,        intent(in) :: f1(:, :)
    CMPLX                    :: spin(1:3)

    CMPLX :: z

    PUSH_SUB(state_spin)

    z = zmf_dotp(mesh, f1(:, 1) , f1(:, 2))

    spin(1) = M_TWO*z
    spin(2) = M_TWO*aimag(z)
    spin(3) = zmf_dotp(mesh, f1(:, 1), f1(:, 1)) - zmf_dotp(mesh, f1(:, 2), f1(:, 2))
    spin = M_HALF*spin ! spin is half the sigma matrix.

    POP_SUB(state_spin)
  end function state_spin

  ! ---------------------------------------------------------
  logical function state_is_local(st, ist)
    type(states_t), intent(in) :: st
    integer,        intent(in) :: ist

    PUSH_SUB(state_is_local)

    state_is_local = ist.ge.st%st_start.and.ist.le.st%st_end

    POP_SUB(state_is_local)
  end function state_is_local


  ! ---------------------------------------------------------

  real(8) function states_wfns_memory(st, mesh) result(memory)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: mesh

    PUSH_SUB(states_wfns_memory)
    memory = 0.0_8

    ! orbitals
    memory = memory + REAL_PRECISION*dble(mesh%np_part_global)*st%d%dim*dble(st%nst)*st%d%kpt%nglobal

    POP_SUB(states_wfns_memory)
  end function states_wfns_memory

  ! ---------------------------------------------------------

  subroutine states_blacs_blocksize(st, mesh, blocksize, total_np)
    type(states_t),  intent(in)    :: st
    type(mesh_t),    intent(in)    :: mesh
    integer,         intent(out)   :: blocksize(2)
    integer,         intent(out)   :: total_np

#ifdef HAVE_SCALAPACK
    ! We need to select the block size of the decomposition. This is
    ! tricky, since not all processors have the same number of
    ! points.
    !
    ! What we do for now is to use the maximum of the number of
    ! points and we set to zero the remaining points.

    if (mesh%parallel_in_domains) then
      blocksize(1) = maxval(mesh%vp%np_local) + &
        (st%d%dim - 1)*maxval(mesh%vp%np_local + mesh%vp%np_bndry + mesh%vp%np_ghost)
    else 
      blocksize(1) = mesh%np + (st%d%dim - 1)*mesh%np_part
    end if

    if (st%parallel_in_states) then
      blocksize(2) = maxval(st%st_num)
    else
      blocksize(2) = st%nst
    end if
    
    total_np = blocksize(1)*st%dom_st_proc_grid%nprow
  

    ASSERT(st%d%dim*mesh%np_part >= blocksize(1))
#endif

  end subroutine states_blacs_blocksize

end module states_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
