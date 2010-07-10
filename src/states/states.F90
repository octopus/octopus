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
  use batch_m
  use blas_m
  use calc_mode_m
  use datasets_m
  use derivatives_m
  use distributed_m
  use geometry_m
  use global_m
  use grid_m
  use hardware_m
  use io_m
  use io_function_m
  use kpoints_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use modelmb_particles_m
  use mpi_m
  use mpi_lib_m
  use multicomm_m
  use multigrid_m
  use ob_green_m
  use ob_interface_m
  use profiling_m
  use simul_box_m
  use smear_m
  use states_dim_m
  use symmetrizer_m
  use unit_m
  use unit_system_m
  use types_m
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
    states_write_eigenvalues,         &
    states_write_dos,                 &
    states_write_tpa,                 &
    states_write_bands,               &
    states_write_fermi_energy,        &
    states_spin_channel,              &
    states_dens_accumulate,           &
    states_dens_accumulate_batch,     &
    states_dens_reduce,               &
    states_calc_dens,                 &
    states_calc_quantities,           &
    state_is_local,                   &
    states_dump,                      &
    states_distribute_nodes,          &
    states_wfns_memory,               &
    states_freeze_orbitals,           &
    states_total_density,             &
    states_init_self_energy,          &
    states_write_proj_lead_wf,        &
    states_read_proj_lead_wf

  public ::                           &
    states_are_complex,               &
    states_are_real,                  &
    states_set_complex

  type states_lead_t
    CMPLX, pointer     :: intf_psi(:, :, :, :) !< (np, st%d%dim, st%nst, st%d%nik)
    FLOAT, pointer     :: rho(:, :)   !< Density of the lead unit cells.
    CMPLX, pointer     :: self_energy(:, :, :, :, :) !< (np, np, nspin, ncs, nik) self-energy of the leads.
  end type states_lead_t

  type states_priv_t
    private
    integer :: wfs_type              !< real (TYPE_FLOAT) or complex (TYPE_CMPLX) wavefunctions
  end type states_priv_t

  type states_t
    type(states_dim_t)       :: d
    type(modelmb_particle_t) :: modelmbparticles
    type(states_priv_t)      :: priv !< the private components 
    integer :: nst                   !< Number of states in each irreducible subspace

    ! pointers to the wavefunctions 
    logical :: only_userdef_istates !< only use user-defined states as initial states in propagation
    FLOAT, pointer :: dpsi(:,:,:,:) !< dpsi(sys%gr%mesh%np_part, st%d%dim, st%nst, st%d%nik)
    CMPLX, pointer :: zpsi(:,:,:,:) !< zpsi(sys%gr%mesh%np_part, st%d%dim, st%nst, st%d%nik)

    logical             :: open_boundaries
    CMPLX, pointer      :: zphi(:, :, :, :)  !< Free states for open-boundary calculations.
    FLOAT, pointer      :: ob_eigenval(:, :) !< Eigenvalues of free states.
    type(states_dim_t)  :: ob_d              !< Dims. of the unscattered systems.
    integer             :: ob_nst            !< nst of the unscattered systems.
    integer             :: ob_ncs            !< No. of continuum states of open system.
                                             !< ob_ncs = ob_nst*st%ob_d%nik / st%d%nik
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
    type(mpi_grp_t)             :: dom_st             !< The MPI group related to the domain-states "plane".
    integer                     :: lnst               !< Number of states on local node.
    integer                     :: st_start, st_end   !< Range of states processed by local node.
    integer, pointer            :: node(:)            !< To which node belongs each state.
    integer, pointer            :: st_range(:, :)     !< Node r manages states st_range(1, r) to
                                                      !< st_range(2, r) for r = 0, ..., mpi_grp%size-1,
                                                      !< i. e. st_start = st_range(1, r) and
                                                      !< st_end = st_range(2, r) on node r.
    integer, pointer            :: st_num(:)          !< Number of states on node r, i. e.
                                                      !< st_num(r) = st_num(2, r)-st_num(1, r).
    type(multicomm_all_pairs_t) :: ap                 !< All-pairs schedule.

    logical                     :: symmetrize_density
  end type states_t

contains

  ! ---------------------------------------------------------
  subroutine states_null(st)
    type(states_t), intent(inout) :: st

    integer :: il

    call push_sub('states.states_null')

    nullify(st%dpsi, st%zpsi, st%zphi, st%rho, st%current, st%rho_core, st%frozen_rho, st%eigenval)
    do il=1, NLEADS
      nullify(st%ob_lead(il)%intf_psi, st%ob_lead(il)%rho, st%ob_lead(il)%self_energy)
    end do
    nullify(st%ob_eigenval, st%ob_occ)
    nullify(st%occ, st%spin, st%node, st%user_def_states)
    nullify(st%d%kweights)
    nullify(st%st_range, st%st_num)

    ! By default, calculations use real wavefunctions
    st%priv%wfs_type = TYPE_FLOAT

    call modelmb_particles_nullify(st%modelmbparticles)

    call pop_sub('states.states_null')
  end subroutine states_null


  ! ---------------------------------------------------------
  subroutine states_init(st, gr, geo)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo

    FLOAT :: excess_charge
    integer :: nempty, ierr, il
    integer, allocatable :: ob_k(:), ob_st(:), ob_d(:)
    character(len=256)   :: restart_dir

    call push_sub('states.states_init')

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
    !% be oriented non-collinearly - <i>i.e.</i> the magnetization vector is allowed to take different
    !% directions at different points.
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
      call write_fatal(2)
    end if

    st%extrastates = (nempty > 0)

    ! For non-periodic systems this should just return the Gamma point
    call states_choose_kpoints(st%d, gr%sb, geo)

    call geometry_val_charge(geo, st%val_charge)
    
    if(gr%ob_grid%open_boundaries) excess_charge = -st%val_charge

    st%qtot = -(st%val_charge + excess_charge)

    do il = 1, NLEADS
      nullify(st%ob_lead(il)%intf_psi)
    end do
    st%open_boundaries = .false.
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
          ! then count the occupied ones
        call states_look(restart_dir, mpi_world, ob_k(il), ob_d(il), ob_st(il), ierr, .true.)
        if(ierr.ne.0) then
          message(1) = 'Could not read the number of states of the periodic calculation'
          message(2) = 'from '//restart_dir//'.'
          call write_fatal(2)
        end if
      end do
      if(NLEADS.gt.1) then
        if(ob_k(LEFT).ne.ob_k(RIGHT).or. &
          ob_st(LEFT).ne.ob_st(LEFT).or. &
          ob_d(LEFT).ne.ob_d(RIGHT)) then
          message(1) = 'The number of states for the left and right leads are not equal.'
          call write_fatal(1)
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
        call write_fatal(2)
      end if
      ! If the system is spin-polarized one half of the free states
      ! goes to the spin-up k-index, the other half to the spin-down
      ! k-index; we therefore divide by st%d%nik.
      st%ob_ncs = st%ob_d%nik*st%ob_nst / st%d%nik
      st%ob_ncs = 1
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
      st%ob_ncs   = 0
      st%ob_d%nik = 0
      st%ob_d%dim = 0
    end if

    select case(st%d%ispin)
    case(UNPOLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty + st%ob_ncs
      st%d%nspin = 1
      st%d%spin_channels = 1
    case(SPIN_POLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty + st%ob_ncs
      ! st%d%nik = st%d%nik*2
      st%d%nspin = 2
      st%d%spin_channels = 2
    case(SPINORS)
      st%d%dim = 2
      st%nst = int(st%qtot)
      if(st%nst < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty + st%ob_ncs
      st%d%nspin = 4
      st%d%spin_channels = 2
    end select

    ! FIXME: For now, open-boundary calculations are only possible for
    ! continuum states, i.e. for those states treated by the Lippmann-
    ! Schwinger approach during SCF.
    if(gr%ob_grid%open_boundaries) then
      if(st%nst.ne.st%ob_nst .or. st%d%nik.ne.st%ob_d%nik) then
        message(1) = 'Open-boundary calculations for possibly bound states'
        message(2) = 'are not possible yet. You have to match your number'
        message(3) = 'of states to the number of free states of your previous'
        message(4) = 'periodic run.'
        write(message(5), '(a,i5,a)') 'Your finite system contributes ', st%nst-st%ob_nst, ' states,'
        write(message(6), '(a,i5,a)') 'while your periodic calculation had ', st%ob_nst, ' states.'
        write(message(7), '(a,i5,a)') 'Your finite system contributes ', st%d%nik-st%ob_d%nik, ' k-points,'
        write(message(8), '(a,i5,a)') 'while your periodic calculation had ', st%ob_d%nik, ' k-points.'
        call write_fatal(8)
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
      call messages_devel_version('Current DFT')

      ! Use of CDFT requires complex wavefunctions
      st%priv%wfs_type = TYPE_CMPLX

      if(st%d%ispin == SPINORS) then
        message(1) = "Sorry, current DFT not working yet for spinors."
        call write_fatal(1)
      end if
      message(1) = "Info: Using current DFT"
      call write_info(1)
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

    if(st%symmetrize_density) call messages_devel_version("Symmetrization of the density")

    call pop_sub('states.states_init')

  contains

    subroutine read_ob_eigenval_and_occ()
      integer            :: occs, jk(1:st%ob_nst), ist, ik, idim, err
      FLOAT              :: flt, eigenval, occ, kweights
      character          :: char
      character(len=256) :: restart_dir, line, chars

      call push_sub('states.states_init.read_ob_eigenval_and_occ')

      restart_dir = trim(gr%ob_grid%lead(LEFT)%info%restart_dir)//'/'//GS_DIR

      occs = io_open(trim(restart_dir)//'/occs', action='read', is_tmp=.true., grp=mpi_world)
      if(occs.lt.0) then
        message(1) = 'Could not read '//trim(restart_dir)//'/occs.'
        call write_fatal(1)
      end if

      ! Skip two lines.
      call iopar_read(mpi_world, occs, line, err)
      call iopar_read(mpi_world, occs, line, err)

      jk(:) = 1
      do
        ! Check for end of file.
        call iopar_read(mpi_world, occs, line, err)

        read(line, '(a)') char
        if(char.eq.'%') exit
        call iopar_backspace(mpi_world, occs)

        ! Extract eigenvalue.
        call iopar_read(mpi_world, occs, line, err)
        read(line, *) occ, char, eigenval, char, flt, char, flt, char, flt, char, &
          kweights, char, chars, char, ik, char, ist, char, idim
        if(occ > M_EPSILON) then
          if(st%d%ispin.eq.SPIN_POLARIZED) then
              message(1) = 'Spin-Transport not implemented!'
              call write_fatal(1)
            if(is_spin_up(ik)) then
              !FIXME
!              st%ob_eigenval(jst, SPIN_UP) = eigenval
!              st%ob_occ(jst, SPIN_UP)      = occ
            else
!              st%ob_eigenval(jst, SPIN_DOWN) = eigenval
!              st%ob_occ(jst, SPIN_DOWN)      = occ
            end if
          else
            st%ob_eigenval(ist, jk(ist)) = eigenval
            st%ob_occ(ist, jk(ist))      = occ
            st%ob_d%kweights(jk(ist))    = kweights
          end if
          jk(ist) = jk(ist) + 1
        end if
      end do

      call io_close(occs)

      call pop_sub('states.states_init.read_ob_eigenval_and_occ')
    end subroutine read_ob_eigenval_and_occ
  end subroutine states_init


  ! ---------------------------------------------------------
  ! Allocate free states.
  subroutine states_allocate_free_states(st, gr)
    type(states_t), intent(inout) :: st
    type(grid_t),   intent(in)    :: gr

    integer :: il

    call push_sub('states.states_allocate_free_states')

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

    call pop_sub('states.states_allocate_free_states')
  end subroutine states_allocate_free_states


  ! ---------------------------------------------------------
  ! Deallocate free states.
  subroutine states_deallocate_free_states(st, gr)
    type(states_t), intent(inout) :: st
    type(grid_t),   intent(in)    :: gr

    integer :: il

    call push_sub('states.states_deallocate_free_states')

    if(gr%ob_grid%open_boundaries) then
      SAFE_DEALLOCATE_P(st%zphi)
      do il = 1, NLEADS
        SAFE_DEALLOCATE_P(st%ob_lead(il)%rho)
      end do
    end if

    call pop_sub('states.states_deallocate_free_states')
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

    integer :: ik, ist, ispin, nspin, ncols
    type(block_t) :: blk
    FLOAT :: rr, charge

    call push_sub('states.states_read_initial_occs')
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
        do ik = 1, st%d%nik
          do ist = 1, st%nst - ncols
            if(st%d%ispin == UNPOLARIZED) then
              st%occ(ist, ik) = M_TWO
            else
              st%occ(ist, ik) = M_ONE
            end if
          end do
        end do
        do ik = 1, st%d%nik
          do ist = st%nst - ncols + 1, st%nst 
            call parse_block_float(blk, ik-1, ist-1-(st%nst-ncols), st%occ(ist, ik))
          end do
        end do
        call parse_block_end(blk)

      else
        st%fixed_occ = .false.

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

    call smear_init(st%smear, st%d%ispin, st%fixed_occ)

    if(.not. smear_is_semiconducting(st%smear) .and. .not. st%extrastates) then
      message(1) = "Non-semiconductor smearing needs ExtraStates > 0 to be useful."
      call write_warning(1)
    endif

    ! sanity check
    charge = M_ZERO
    do ist = 1, st%nst
      charge = charge + sum(st%occ(ist, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
    end do
    if(abs(charge - st%qtot) > CNST(1e-6)) then
      message(1) = "Initial occupations do not integrate to total charge."
      write(message(2), '(6x,f12.6,a,f12.6)') charge, ' != ', st%qtot
      call write_fatal(2)
    end if

    call pop_sub('states.states_read_initial_occs')
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

    call push_sub('states.states_read_initial_spins')

    st%fixed_spins = .false.
    if(st%d%ispin .ne. SPINORS) then
      call pop_sub('states.states_read_initial_spins')
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

    call pop_sub('states.states_read_initial_spins')
  end subroutine states_read_initial_spins


  ! ---------------------------------------------------------
  ! Allocates the KS wavefunctions defined within a states_t structure.
  subroutine states_allocate_wfns(st, mesh, wfs_type, ob_mesh)
    type(states_t),    intent(inout)   :: st
    type(mesh_t),      intent(in)      :: mesh
    integer, optional, intent(in)      :: wfs_type
    type(mesh_t), optional, intent(in) :: ob_mesh(:)

    integer :: ip, ik, ist, idim, st1, st2, k1, k2, size, il
    logical :: force

    call push_sub('states.states_allocate_wfns')

    if(associated(st%dpsi).or.associated(st%zpsi)) then
      message(1) = "Trying to allocate wavefunctions that are already allocated."
      call write_fatal(1)
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

    if(present(ob_mesh).and.calc_mode_is(CM_TD).and.st%open_boundaries) then
      do il = 1, NLEADS
        SAFE_ALLOCATE(st%ob_lead(il)%intf_psi(1:ob_mesh(il)%np, 1:st%d%dim, st1:st2, k1:k2))
        st%ob_lead(il)%intf_psi = M_z0
      end do

    end if

    call pop_sub('states.states_allocate_wfns')
  end subroutine states_allocate_wfns


  ! ---------------------------------------------------------
  !> Deallocates the KS wavefunctions defined within a states_t structure.
  subroutine states_deallocate_wfns(st)
    type(states_t), intent(inout) :: st

    integer :: il

    call push_sub('states.states_deallocate_wfns')

    if (states_are_real(st)) then
      SAFE_DEALLOCATE_P(st%dpsi)
    else
      SAFE_DEALLOCATE_P(st%zpsi)
    end if

    if(st%open_boundaries .and. calc_mode_is(CM_TD)) then
      do il = 1, NLEADS
        SAFE_DEALLOCATE_P(st%ob_lead(il)%intf_psi)
      end do
    end if

    call pop_sub('states.states_deallocate_wfns')
  end subroutine states_deallocate_wfns


  ! ---------------------------------------------------------
  subroutine states_densities_init(st, gr, geo, mc)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo
    type(multicomm_t), intent(in)    :: mc

    call push_sub('states.states_densities_init')

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

    call pop_sub('states.states_densities_init')
    
  contains

    subroutine states_exec_init()

      !%Variable StatesBlockSize
      !%Type integer
      !%Default max(4, 2*nthreads)
      !%Section Execution::Optimization
      !%Description
      !% Some routines work over blocks of eigenfunctions, this
      !% generally improves performance at the expense of increased
      !% memory consumption. This variable selects the size of the
      !% blocks to be used.
      !%End

      call parse_integer(datasets_check('StatesBlockSize'), max(4, 2*mc%nthreads), st%d%block_size)
      if(st%d%block_size < 1) then
        message(1) = "Error: The variable 'StatesBlockSize' must be greater than 0."
        call write_fatal(1)
      end if

      st%d%block_size = min(st%d%block_size, st%nst)

      !%Variable StatesWindowSize
      !%Type integer
      !%Default nst
      !%Section Execution::Optimization
      !%Description
      !% (experimental) The orthogonalization of the wavefunctions can
      !% be done by windows; this value selects the size of the
      !% windows. The default size is total number of states, which
      !% disables window orthogonalization.
      !%End

      call parse_integer(datasets_check('StatesWindowSize'), st%nst, st%d%window_size)
      if(st%d%block_size < 1) then
        message(1) = "Error: The variable 'StatesBlockSize' must be greater than 0."
        call write_fatal(1)
      end if

      st%d%window_size = min(st%d%window_size, st%nst)

      !%Variable StatesOrthogonalization
      !%Type integer
      !%Default gs
      !%Section Execution::Optimization
      !%Description
      !% Select the full orthogonalization method used by some
      !% eigensolvers. The default is gram_schmidt.
      !%Option gram_schmidt 1
      !% The standard Gram-Schmidt orthogonalization implemented using
      !% Blas level 3 routines.
      !%Option mgs 2
      !% Modified Gram-Schmidt orthogonalization.
      !%Option qr 3
      !% (Experimental) Orthogonalization is performed based on a QR
      !% decomposition based on Lapack routines _getrf and _orgqr.
      !%End

      call parse_integer(datasets_check('StatesOrthogonalization'), ORTH_GS, st%d%orth_method)
      if(st%d%orth_method == ORTH_QR) call messages_devel_version("QR Orthogonalization")

    end subroutine states_exec_init
  end subroutine states_densities_init


  ! ---------------------------------------------------------
  subroutine states_copy(stout, stin)
    type(states_t), intent(inout) :: stout
    type(states_t), intent(in)    :: stin

    call push_sub('states.states_copy')

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

    call pop_sub('states.states_copy')
  end subroutine states_copy


  ! ---------------------------------------------------------
  subroutine states_end(st)
    type(states_t), intent(inout) :: st

    integer :: il

    call push_sub('states.states_end')
    
    SAFE_DEALLOCATE_P(st%dpsi)
    SAFE_DEALLOCATE_P(st%zpsi)
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
    if(st%open_boundaries) then
      call states_dim_end(st%ob_d)
      SAFE_DEALLOCATE_P(st%ob_eigenval)
      SAFE_DEALLOCATE_P(st%ob_occ)
      do il = 1, NLEADS
        SAFE_DEALLOCATE_P(st%ob_lead(il)%self_energy)
      end do
    end if

    SAFE_DEALLOCATE_P(st%user_def_states)

    call modelmb_particles_end(st%modelmbparticles)

    call pop_sub('states.states_end')
  end subroutine states_end


  ! ---------------------------------------------------------
  !> Calculates the new density out the wavefunctions and
  !! occupations...
  subroutine states_dens_accumulate(st, gr, ist, ik, rho)
    type(states_t), intent(inout) :: st
    type(grid_t),   intent(in)    :: gr
    integer,        intent(in)    :: ist
    integer,        intent(in)    :: ik
    FLOAT,          intent(inout) :: rho(:,:)
    
    type(batch_t) :: psib

    call push_sub('states.states_dens_accumulate')

    if(states_are_real(st)) then
      call batch_init(psib, st%d%dim, ist, ist, st%dpsi(:, :, ist:ist, ik))
    else
      call batch_init(psib, st%d%dim, ist, ist, st%zpsi(:, :, ist:ist, ik))
    end if

    call states_dens_accumulate_batch(st, gr, ik, psib, rho)

    call batch_end(psib)

    call pop_sub('states.states_dens_accumulate')
  end subroutine states_dens_accumulate

  ! ---------------------------------------------------

  subroutine states_dens_accumulate_batch(st, gr, ik, psib, rho)
    type(states_t), intent(in)    :: st
    type(grid_t),   intent(in)    :: gr
    integer,        intent(in)    :: ik
    type(batch_t),  intent(inout) :: psib
    FLOAT, target,  intent(inout) :: rho(:,:)
    
    integer :: ist, ist2, ip, ispin
    CMPLX   :: term, psi1, psi2
    FLOAT, pointer :: dpsi(:, :)
    CMPLX, pointer :: zpsi(:, :)
    FLOAT, pointer :: crho(:)
    FLOAT, allocatable :: frho(:)
    type(profile_t), save :: prof

    call push_sub('states.states_dens_accumulate_batch')
    call profiling_in(prof, "CALC_DENSITY")

    ASSERT(ubound(rho, dim = 1) == gr%fine%mesh%np .or. ubound(rho, dim = 1) == gr%fine%mesh%np_part)

    ispin = states_dim_get_spin_index(st%d, ik)

    if(gr%have_fine_mesh) then
      SAFE_ALLOCATE(crho(1:gr%mesh%np_part))
      crho = M_ZERO
    else
      crho => rho(:, ispin)
    end if

    if(states_are_real(st)) then
      do ist = 1, psib%nst
        ist2 = psib%states(ist)%ist
        dpsi => psib%states(ist)%dpsi

        forall(ip = 1:gr%mesh%np)
          crho(ip) = crho(ip) + st%d%kweights(ik) * st%occ(ist2, ik) * dpsi(ip, 1)**2
        end forall
      end do
    else
      do ist = 1, psib%nst
        ist2 = psib%states(ist)%ist
        zpsi => psib%states(ist)%zpsi

        forall(ip = 1:gr%mesh%np)
          crho(ip) = crho(ip) + st%d%kweights(ik) * st%occ(ist2, ik) * &
            (real(zpsi(ip, 1), REAL_PRECISION)**2 + aimag(zpsi(ip, 1))**2)
        end forall
      end do
    end if
    
    if(gr%have_fine_mesh) then
      SAFE_ALLOCATE(frho(1:gr%fine%mesh%np))
      call dmultigrid_coarse2fine(gr%fine%tt, gr%der, gr%fine%mesh, crho, frho, order = 2)
      ! some debugging output that I will keep here for the moment, XA
      !      call doutput_function(1, "./", "n_fine", gr%fine%mesh, frho, unit_one, ierr)
      !      call doutput_function(1, "./", "n_coarse", gr%mesh, crho, unit_one, ierr)
      forall(ip = 1:gr%fine%mesh%np) rho(ip, ispin) = rho(ip, ispin) + frho(ip)
      SAFE_DEALLOCATE_P(crho)
      SAFE_DEALLOCATE_A(frho)
    end if

    if(st%d%ispin == SPINORS) then ! in this case wavefunctions are always complex

      ASSERT(.not. gr%have_fine_mesh)

     do ist = 1, psib%nst
        ist2 = psib%states(ist)%ist
        zpsi => psib%states(ist)%zpsi

        do ip = 1, gr%fine%mesh%np
          
          psi1 = zpsi(ip, 1)
          psi2 = zpsi(ip, 2)
          
          rho(ip, 2) = rho(ip, 2) + &
            st%d%kweights(ik) * st%occ(ist2, ik) * (real(psi2, REAL_PRECISION)**2 + aimag(psi2)**2)
        
          term = st%d%kweights(ik) * st%occ(ist2, ik) * psi1 * conjg(psi2)
          rho(ip, 3) = rho(ip, 3) + real(term, REAL_PRECISION)
          rho(ip, 4) = rho(ip, 4) + aimag(term)

        end do
      end do
    end if
    
    call profiling_out(prof)

    call pop_sub('states.states_dens_accumulate_batch')
  end subroutine states_dens_accumulate_batch

  ! ---------------------------------------------------

  subroutine states_dens_reduce(st, gr, rho)
    type(states_t), intent(in)    :: st
    type(grid_t),   intent(in)    :: gr
    FLOAT,          intent(inout) :: rho(:,:)

    type(symmetrizer_t) :: symmetrizer
    FLOAT,  allocatable :: symmrho(:)
    integer :: ispin, np
#ifdef HAVE_MPI
    FLOAT,  allocatable :: reduce_rho(:)
    type(profile_t), save :: reduce_prof
#endif

    np = gr%fine%mesh%np

    call push_sub('states.states_dens_reduce')

#ifdef HAVE_MPI
    ! reduce over states
    if(st%parallel_in_states) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
#ifndef HAVE_MPI2
      SAFE_ALLOCATE(reduce_rho(1:np))
#endif
      do ispin = 1, st%d%nspin
#ifndef HAVE_MPI2
        call blas_copy(np, rho(1, ispin), 1, reduce_rho(1), 1)
#endif
        call MPI_Allreduce(MPI_IN_PLACE_OR(reduce_rho(1)), rho(1, ispin), np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      end do
      SAFE_DEALLOCATE_A(reduce_rho)
      call profiling_out(reduce_prof)
    end if

    ! reduce over k-points
    if(st%d%kpt%parallel) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
      SAFE_ALLOCATE(reduce_rho(1:np))
      do ispin = 1, st%d%nspin
        call MPI_Allreduce(rho(1, ispin), reduce_rho(1), np, MPI_FLOAT, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
        call blas_copy(np, reduce_rho(1), 1, rho(1, ispin), 1)
      end do
      SAFE_DEALLOCATE_A(reduce_rho)
      call profiling_out(reduce_prof)
    end if
#endif

    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmrho(1:np))
      call symmetrizer_init(symmetrizer, gr%fine%mesh)

      do ispin = 1, st%d%nspin
        call dsymmetrizer_apply(symmetrizer, rho(:, ispin), symmrho)
        rho(1:np, ispin) = symmrho(1:np)
      end do

      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmrho)
    end if

    call pop_sub('states.states_dens_reduce')
  end subroutine states_dens_reduce


  ! ---------------------------------------------------------
  !> Computes the density from the orbitals in st. If rho is
  !! present, the density is placed there; if it is not present,
  !! the density is placed in st%rho.
  ! ---------------------------------------------------------
  subroutine states_calc_dens(st, gr, rho)
    type(states_t),          intent(inout)  :: st
    type(grid_t),            intent(in)     :: gr
    FLOAT, optional, target, intent(out)    :: rho(:,:)

    integer :: ik
    FLOAT, pointer :: dens(:, :)
    type(batch_t)  :: psib

    call push_sub('states.states_calc_dens')

    if(present(rho)) then
      dens => rho
    else
      dens => st%rho
    end if

    ASSERT(ubound(dens, dim = 1) == gr%fine%mesh%np .or. ubound(dens, dim = 1) == gr%fine%mesh%np_part)

    dens(1:gr%fine%mesh%np, 1:st%d%nspin) = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      if(states_are_real(st)) then
        call batch_init(psib, st%d%dim, st%st_start, st%st_end, st%dpsi(:, :, st%st_start:, ik))
      else
        call batch_init(psib, st%d%dim, st%st_start, st%st_end, st%zpsi(:, :, st%st_start:, ik))
      end if

      call states_dens_accumulate_batch(st, gr, ik, psib, dens)
      
      call batch_end(psib)
    end do

    call states_dens_reduce(st, gr, dens)

    nullify(dens)
    call pop_sub('states.states_calc_dens')
  end subroutine states_calc_dens

  ! ---------------------------------------------------------
  !> generate a hydrogen s-wavefunction around a random point
  subroutine states_generate_random(st, mesh, ist_start_, ist_end_)
    type(states_t),    intent(inout) :: st
    type(mesh_t),      intent(in)    :: mesh
    integer, optional, intent(in)    :: ist_start_, ist_end_

    integer :: ist, ik, id, ist_start, ist_end, jst, seed
    CMPLX   :: alpha, beta

    call push_sub('states.states_generate_random')

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

      ASSERT(st%priv%wfs_type == TYPE_CMPLX)

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

    call pop_sub('states.states_generate_random')
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

    call push_sub('states.states_fermi')

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
      call write_warning(2)
    end if

    if(st%d%ispin == SPINORS) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          if (st%priv%wfs_type == TYPE_FLOAT) then
            write(message(1),'(a)') 'Internal error in states_fermi.'
            call write_fatal(1)
          else
            st%spin(1:3, ist, ik) = state_spin(mesh, st%zpsi(:, :, ist, ik))
          end if
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

    call pop_sub('states.states_fermi')
  end subroutine states_fermi


  ! ---------------------------------------------------------
  !> function to calculate the eigenvalues sum using occupations as weights
  function states_eigenvalues_sum(st, alt_eig) result(tot)
    type(states_t), intent(in)  :: st
    FLOAT, optional, intent(in) :: alt_eig(st%st_start:st%st_end, 1:st%d%nik)
    FLOAT                       :: tot

    integer :: ik
#ifdef HAVE_MPI
    FLOAT :: tot_temp
#endif

    call push_sub('states.states_eigenvalues_sum')

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

#ifdef HAVE_MPI
    if(st%parallel_in_states) then
      call MPI_Allreduce(tot, tot_temp, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      tot = tot_temp
    end if

    if(st%d%kpt%parallel) then
      call MPI_Allreduce(tot, tot_temp, 1, MPI_FLOAT, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
      tot = tot_temp
    end if
#endif

    call pop_sub('states.states_eigenvalues_sum')
  end function states_eigenvalues_sum


  ! ---------------------------------------------------------
  subroutine states_write_eigenvalues(iunit, nst, st, sb, error)
    integer,           intent(in) :: iunit, nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb
    FLOAT, optional,   intent(in) :: error(nst, st%d%nik)

    integer ik, ist, ns, is
    FLOAT :: occ, kpoint(1:3)
    character(len=80) tmp_str(MAX_DIM), cspin

    call push_sub('states.states_write_eigenvalues')

    ns = 1
    if(st%d%nspin == 2) ns = 2

    message(1) = 'Eigenvalues [' // trim(units_abbrev(units_out%energy)) // ']'
    call write_info(1, iunit)
    if (st%d%nik > ns) then
      message(1) = 'k-points [' // trim(units_abbrev(unit_one/units_out%length)) //']'
      call write_info(1, iunit)
    end if

    if(.not. mpi_grp_is_root(mpi_world)) then
      call pop_sub('states.states_write_eigenvalues')
      return
    end if

    if(present(error)) then
      if(st%d%ispin .eq. SPINORS) then
        write(message(1), '(a4,1x,a5,1x,a12,1x,a12,2x,a4,4x,a4,4x,a4,5x,a5)')   &
          '#st',' Spin',' Eigenvalue', 'Occupation ', '<Sx>', '<Sy>', '<Sz>', 'Error'
      else
        write(message(1), '(a4,1x,a5,1x,a12,4x,a12,1x,a10)')   &
          '#st',' Spin',' Eigenvalue', 'Occupation ', 'Error'
      end if
    else
      if(st%d%ispin .eq. SPINORS) then
        write(message(1), '(a4,1x,a5,1x,a12,1x,a12,2x,a4,4x,a4,4x,a4)')   &
          '#st',' Spin',' Eigenvalue', 'Occupation ', '<Sx>', '<Sy>', '<Sz>'
      else
        write(message(1), '(a4,1x,a5,1x,a12,4x,a12,1x)')       &
          '#st',' Spin',' Eigenvalue', 'Occupation '
      end if
    end if
    call write_info(1, iunit)

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then
        kpoint = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik))
        kpoint = units_from_atomic(unit_one/units_out%length, kpoint)
        write(message(1), '(a,i4,3(a,f12.6),a)') '#k =', ik, ', k = (',  kpoint(1), ',', kpoint(2), ',', kpoint(3), ')'
        call write_info(1, iunit)
      end if

      do ist = 1, nst
        do is = 0, ns-1
          if(ist > st%nst) then
            occ = M_ZERO
          else
            occ = st%occ(ist, ik+is)
          end if

          if(is .eq. 0) cspin = 'up'
          if(is .eq. 1) cspin = 'dn'
          if(st%d%ispin .eq. UNPOLARIZED .or. st%d%ispin .eq. SPINORS) cspin = '--'

          write(tmp_str(1), '(i4,3x,a2)') ist, trim(cspin)
          if(simul_box_is_periodic(sb)) then
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,3x,4f5.2)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik)), occ, st%spin(1:3, ist, ik)
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(ist, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik+is)), occ
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(ist, ik), ')'
            end if
          else
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,5x,f5.2,3x,3f8.4)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik)), occ, st%spin(1:3, ist, ik)
              if(present(error)) write(tmp_str(3), '(a3,es7.1,a1)')'  (', error(ist, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik+is)), occ
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(ist, ik), ')'
            end if
          end if
          if(present(error)) then
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))//trim(tmp_str(3))
          else
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))
          end if
          call write_info(1, iunit)
        end do
      end do
    end do

    call pop_sub('states.states_write_eigenvalues')
  end subroutine states_write_eigenvalues


  ! ---------------------------------------------------------
  subroutine states_write_bands(dir, nst, st, sb)
    character(len=*),  intent(in) :: dir    
    integer,           intent(in) :: nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb

    integer :: idir, ist, ik, ns, is
    integer, allocatable :: iunit(:)
    FLOAT   :: factor(MAX_DIM), kpoint(1:MAX_DIM)
    logical :: grace_mode, gnuplot_mode
    character(len=80) :: filename    

    if(.not.mpi_grp_is_root(mpi_world)) return

    call push_sub('states.states_write_bands')

    !%Variable OutputBandsGnuplotMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in Gnuplot-friendly format.
    !%End
    call parse_logical(datasets_check('OutputBandsGnuplotMode'), .true., gnuplot_mode)

    !%Variable OutputBandsGraceMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in Grace-friendly format.
    !%End
    call parse_logical(datasets_check('OutputBandsGraceMode'), .false., grace_mode)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    SAFE_ALLOCATE(iunit(0:ns-1))

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do idir = 1, MAX_DIM
      factor(idir) = M_ONE
      if (sb%klattice(idir, idir) /= M_ZERO) factor(idir) = sb%klattice(idir, idir)
    end do

    if (gnuplot_mode) then
      do is = 0, ns-1
        if (ns .gt. 1) then
          write(filename, '(a,i1.1,a)') 'bands-gp-', is+1,'.dat'
        else
          write(filename, '(a)') 'bands-gp.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a, i6)') '# kx ky kz (unscaled), kx ky kz (scaled), bands:', nst
      end do

      ! output bands in gnuplot format
      do ist = 1, nst
        do ik = 1, st%d%nik, ns
          do is = 0, ns - 1
            kpoint = M_ZERO
            kpoint = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik + is))
            write(iunit(is), '(1x,6f14.8,3x,f14.8)')            &
              kpoint(1:sb%dim),                                 & ! unscaled
              kpoint(1:sb%dim)/factor(1:sb%dim),                & ! scaled
              units_from_atomic(units_out%energy, st%eigenval(ist, ik + is))
          end do
        end do
        do is = 0, ns-1
          write(iunit(is), '(a)') ''
        end do
      end do
      do is = 0, ns-1
        call io_close(iunit(is))
      end do
    end if

    if (grace_mode) then
      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i1.1,a)') 'bands-grace-', is+1,'.dat'
        else
          write(filename, '(a)') 'bands-grace.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a, i6)') '# kx ky kz (unscaled), kx ky kz (scaled), bands:', nst
      end do

      ! output bands in xmgrace format, i.e.:
      ! k_x, k_y, k_z, e_1, e_2, ..., e_n
      do ik = 1, st%d%nik, ns
        do is = 0, ns-1
          kpoint = M_ZERO
          kpoint = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik + is))
          write(iunit(is), '(1x,6f14.8,3x,16384f14.8)')         &
            kpoint(1:MAX_DIM),                                  & ! unscaled
            kpoint(1:MAX_DIM)/factor(1:MAX_DIM),                & ! scaled
            (units_from_atomic(units_out%energy, st%eigenval(ist, ik+is)), ist = 1, nst)
        end do
      end do
      do is = 0, ns-1
        call io_close(iunit(is))
      end do        
    end if

    SAFE_DEALLOCATE_A(iunit)

    call pop_sub('states.states_write_bands')
  end subroutine states_write_bands

  ! ---------------------------------------------------------
  subroutine states_write_tpa(dir, gr, st)
    character(len=*), intent(in) :: dir
    type(grid_t),     intent(in) :: gr
    type(states_t),   intent(in) :: st

    type(block_t) :: blk
    integer       :: ncols, icoord, ist, ik, tpa_initialst, tpa_initialk
    integer       :: iunit

    FLOAT, allocatable  :: ff(:)
    CMPLX, allocatable  :: cff(:)
    FLOAT, allocatable  :: osc(:)
    FLOAT               :: transition_energy, osc_strength, dsf

    FLOAT, parameter    :: M_THRESHOLD = CNST(1.0e-6)
    logical             :: use_qvector = .false.
    FLOAT, allocatable  :: qvector(:)

    call push_sub('states.states_write_tpa')

    ! find the orbital with half-occupation
    tpa_initialst = -1
    do ist = 1, st%nst
      do ik = 1, st%d%nik
        if (abs(st%occ(ist,ik)-0.5) .lt. M_THRESHOLD) then
          tpa_initialst = ist
          tpa_initialk  = ik
        end if
      end do
    end do

    ! make sure that half-occupancy was found
    if(tpa_initialst.eq.-1) then
      if(mpi_grp_is_root(mpi_world)) then

        message(1) = 'No orbital with half-occupancy found. TPA output is not written.'
        call write_warning(1)

    call pop_sub('states.states_write_tpa')
return

      end if
    end if

    !%Variable MomentumTransfer
    !%Type block
    !%Section States
    !%Description
    !% Momentum-transfer vector <i>q</i> to be used when calculating matrix elements
    !% &lt;f|exp(iq.r)|i&gt;. This enables the calculation of the dynamical structure factor,
    !% which is closely related to generalized oscillator strengths.
    !% If the vector is not given, but TPA output is requested (<tt>Output = TPA</tt>),
    !% only the oscillator strengths are written in the output file.
    !% For example, to use <i>q</i> = (0.1, 0.2, 0.3), set
    !%
    !% <tt>%MomentumTransfer
    !% <br>&nbsp;&nbsp; 0.1 | 0.2 | 0.3
    !% <br>%</tt>
    !%End
    if(parse_block(datasets_check('MomentumTransfer'),blk)==0) then

      ! check if input makes sense
      ncols = parse_block_cols(blk, 0)

      if(ncols .ne. gr%mesh%sb%dim ) then ! wrong size

        if(mpi_grp_is_root(mpi_world)) then
          message(1) = 'Inconsistent size of momentum-transfer vector. It will not be used in the TPA calculation.'
          call write_warning(1)
        end if

      else ! correct size

        use_qvector = .true.
        SAFE_ALLOCATE(qvector(1:gr%mesh%sb%dim))

        do icoord = 1,gr%mesh%sb%dim    !for x,y,z
          call parse_block_float(blk, 0, icoord-1, qvector(icoord))
          qvector(icoord) = units_to_atomic(unit_one / units_inp%length, qvector(icoord))
        end do

      end if

    end if

    ! calculate the matrix elements

    SAFE_ALLOCATE(ff(1:gr%mesh%np))
    if(use_qvector) then
      SAFE_ALLOCATE(cff(1:gr%mesh%np))
    end if
    SAFE_ALLOCATE(osc(1:gr%mesh%sb%dim))

    ! root writes output to file

    if(mpi_grp_is_root(mpi_world)) then

      iunit = io_open(trim(dir)//'/'//trim('tpa_xas'), action='write')    

      ! header
      if(use_qvector) then
        write (message(1),'(a1,a30,3(es14.5,1x),a1)') '#', ' momentum-transfer vector : (', &
          (units_from_atomic(unit_one / units_out%length, qvector(icoord)), icoord=1, gr%mesh%sb%dim),')'
        select case(gr%mesh%sb%dim)
          case(1); write(message(2), '(a1,4(a15,1x))') '#', 'E' , '<x>', '<f>', 'S(q,omega)'
          case(2); write(message(2), '(a1,5(a15,1x))') '#', 'E' , '<x>', '<y>', '<f>', 'S(q,omega)'
          case(3); write(message(2), '(a1,6(a15,1x))') '#', 'E' , '<x>', '<y>', '<z>', '<f>', 'S(q,omega)'
        end select
        call write_info(2,iunit)
      else
        select case(gr%mesh%sb%dim)
          case(1); write(message(1), '(a1,3(a15,1x))') '#', 'E' , '<x>', '<f>'
          case(2); write(message(1), '(a1,4(a15,1x))') '#', 'E' , '<x>', '<y>', '<f>'
          case(3); write(message(1), '(a1,5(a15,1x))') '#', 'E' , '<x>', '<y>', '<z>', '<f>'
        end select
        call write_info(1,iunit)
      end if

    end if

    ! loop through every state
    do ist = 1,st%nst

      ! final states are the unoccupied ones
      if (abs(st%occ(ist,tpa_initialk)) .lt. M_THRESHOLD) then

        osc_strength=M_ZERO
        transition_energy=st%eigenval(ist,tpa_initialk)-st%eigenval(tpa_initialst,tpa_initialk)

        ! dipole matrix elements <f|x|i> etc. -> oscillator strengths
        do icoord=1,gr%mesh%sb%dim    ! for x,y,z

          ff(1:gr%mesh%np) = st%dpsi(1:gr%mesh%np,1,tpa_initialst,tpa_initialk) * &
                       &  gr%mesh%x(1:gr%mesh%np,icoord)                        * &
                       &  st%dpsi(1:gr%mesh%np,1,ist,tpa_initialk)
          osc(icoord)  = dmf_integrate(gr%mesh, ff)
          osc_strength = osc_strength + 2.0/real(gr%mesh%sb%dim)*transition_energy*abs(osc(icoord))**2.0

        end do

        ! matrix elements <f|exp(iq.r)|i> -> dynamic structure factor
        if (use_qvector) then

          cff(1:gr%mesh%np) = cmplx(st%dpsi(1:gr%mesh%np,1,tpa_initialst,tpa_initialk)) * &
                       &   cmplx(st%dpsi(1:gr%mesh%np,1,ist,tpa_initialk))
          do icoord=1,gr%mesh%sb%dim    ! for x,y,z
            cff(1:gr%mesh%np) = cff(1:gr%mesh%np) * exp(M_zI*gr%mesh%x(1:gr%mesh%np,icoord)*qvector(icoord))
          end do

          dsf = abs(zmf_integrate(gr%mesh, cff))**2.0
        end if

        ! write oscillator strengths (+ dynamic structure factor if qvector if given) into file
        if(mpi_grp_is_root(mpi_world)) then

          if(use_qvector) then
            write(message(1), '(1x,6(es15.8,1x))') units_from_atomic(units_out%energy, transition_energy), osc(:), osc_strength, &
                                                   units_from_atomic(unit_one/units_out%energy, dsf)
          else
            write(message(1), '(1x,6(es15.8,1x))') units_from_atomic(units_out%energy, transition_energy), osc(:), osc_strength
          endif

          call write_info(1,iunit)

        end if

      end if

    end do

    ! finally close the file
    if(mpi_grp_is_root(mpi_world)) then
      call io_close(iunit)
    end if

    SAFE_DEALLOCATE_A(ff)
    if(use_qvector) then
      SAFE_DEALLOCATE_A(cff)
    end if
    SAFE_DEALLOCATE_A(osc)
    if (use_qvector) then
      SAFE_DEALLOCATE_A(qvector)
    end if
  
    call pop_sub('states.states_write_tpa')
 
  end subroutine states_write_tpa

  ! ---------------------------------------------------------
  subroutine states_write_dos(dir, st)
    character(len=*), intent(in) :: dir
    type(states_t),   intent(in) :: st

    integer :: ie, ik, ist, epoints, is, ns
    integer, allocatable :: iunit(:)
    FLOAT   :: emin, emax, de, gamma, energy
    FLOAT   :: evalmax, evalmin, tdos, eextend
    FLOAT, allocatable :: dos(:,:,:)
    character(len=64)  :: filename

    call push_sub('states.states_write_dos')

    evalmin = minval(st%eigenval)
    evalmax = maxval(st%eigenval)
    ! we extend the energy mesh by this amount
    eextend  = (evalmax - evalmin) / M_FOUR

    !%Variable DOSEnergyMin
    !%Type float
    !%Section Output
    !%Description
    !% Lower bound for the energy mesh of the DOS.
    !% The default is the lowest eigenvalue, minus a quarter of the total range of eigenvalues.
    !%End
    call parse_float(datasets_check('DOSEnergyMin'), units_from_atomic(units_inp%energy, evalmin - eextend), emin)
    emin = units_to_atomic(units_inp%energy, emin)

    !%Variable DOSEnergyMax
    !%Type float
    !%Section Output
    !%Description
    !% Upper bound for the energy mesh of the DOS.
    !% The default is the highest eigenvalue, plus a quarter of the total range of eigenvalues.
    !%End
    call parse_float(datasets_check('DOSEnergyMax'), units_from_atomic(units_inp%energy, evalmax + eextend), emax)
    emax = units_to_atomic(units_inp%energy, emax)

    !%Variable DOSEnergyPoints
    !%Type integer
    !%Default 500
    !%Section Output
    !%Description
    !% Determines how many energy points <tt>Octopus</tt> should use for 
    !% the DOS energy grid.
    !%End
    call parse_integer(datasets_check('DOSEnergyPoints'), 500, epoints)

    !%Variable DOSGamma
    !%Type float
    !%Default 0.008 Ha
    !%Section Output
    !%Description
    !% Determines the width of the Lorentzian which is used for the DOS sum.
    !%End
    call parse_float(datasets_check('DOSGamma'), &
      units_from_atomic(units_inp%energy, CNST(0.008)), gamma)
    gamma = units_to_atomic(units_inp%energy, gamma)

    ! spacing for energy mesh
    de = (emax - emin) / (epoints - 1)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    ! space for state-dependent DOS
    SAFE_ALLOCATE(dos(1:epoints, 1:st%nst, 0:ns-1))
    SAFE_ALLOCATE(iunit(0:ns-1))    

    ! compute band/spin-resolved density of states
    do ist = 1, st%nst

      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i4.4,a,i1.1,a)') 'dos-', ist, '-', is+1,'.dat'
        else
          write(filename, '(a,i4.4,a)') 'dos-', ist, '.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a)') '# energy, band resolved DOS'
      end do

      do ie = 1, epoints
        energy = emin + (ie - 1) * de
        dos(ie, ist, :) = M_ZERO
        ! sum up Lorentzians
        do ik = 1, st%d%nik, ns
          do is = 0, ns-1
            dos(ie, ist, is) = dos(ie, ist, is) + st%d%kweights(ik+is) * M_ONE/M_Pi * &
              gamma / ( (energy - st%eigenval(ist, ik+is))**2 + gamma**2 )
          end do
        end do
        do is = 0, ns-1
          write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                        units_from_atomic(unit_one / units_out%energy, dos(ie, ist, is))
          call write_info(1, iunit(is))
        end do
      end do

      do is = 0, ns-1
        call io_close(iunit(is))
      end do
    end do

    ! for spin-polarized calculations also output spin-resolved tDOS
    if(st%d%nspin .gt. 1) then    
      do is = 0, ns-1
        write(filename, '(a,i1.1,a)') 'total-dos-', is+1,'.dat'
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a)') '# energy, total DOS (spin-resolved)'

        do ie = 1, epoints
          energy = emin + (ie - 1) * de
          tdos = M_ZERO
          do ist = 1, st%nst
            tdos = tdos + dos(ie, ist, is)
          end do
          write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                        units_from_atomic(unit_one / units_out%energy, tdos)
          call write_info(1, iunit(is))
        end do

        call io_close(iunit(is))
      end do
    end if


    iunit(0) = io_open(trim(dir)//'/'//'total-dos.dat', action='write')    

    ! compute total density of states
    do ie = 1, epoints
      energy = emin + (ie - 1) * de
      tdos = M_ZERO
      do ist = 1, st%nst
        do is = 0, ns-1
          tdos = tdos + dos(ie, ist, is)
        end do
      end do
      write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                    units_from_atomic(unit_one / units_out%energy, tdos)
      call write_info(1, iunit(0))
    end do

    call io_close(iunit(0))

    SAFE_DEALLOCATE_A(iunit)
    SAFE_DEALLOCATE_A(dos)

    call pop_sub('states.states_write_dos')
  end subroutine states_write_dos


  ! ---------------------------------------------------------
  subroutine states_write_fermi_energy(dir, st, mesh, sb)
    character(len=*),  intent(in) :: dir
    type(states_t), intent(inout) :: st
    type(mesh_t),      intent(in) :: mesh
    type(simul_box_t), intent(in) :: sb

    integer :: iunit, idir
    FLOAT :: maxdos
    FLOAT :: factor(MAX_DIM)

    call push_sub('states.states_write_fermi_energy')

    call states_fermi(st, mesh)

    iunit = io_open(trim(dir)//'/'//'bands-efermi.dat', action='write')    

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do idir = 1, MAX_DIM
      factor(idir) = M_ONE
      if (sb%klattice(idir, idir) /= M_ZERO) factor(idir) = sb%klattice(idir, idir)
    end do

    ! write Fermi energy in a format that can be used together 
    ! with bands.dat
    write(message(1), '(a)') '# Fermi energy in a format compatible with bands-gp.dat'

    write(message(2), '(7f12.6)')          &
      minval(sb%kpoints%reduced%point(1,:)),           &
      minval(sb%kpoints%reduced%point(2,:)),           &
      minval(sb%kpoints%reduced%point(3,:)),           &
      minval(sb%kpoints%reduced%point(1,:)/factor(1)), &
      minval(sb%kpoints%reduced%point(2,:)/factor(2)), &
      minval(sb%kpoints%reduced%point(3,:)/factor(3)), &
      units_from_atomic(units_out%energy, st%smear%e_fermi)

    ! Gamma point
    write(message(3), '(7f12.6)')          &
      (M_ZERO, idir = 1, 6),               &
      units_from_atomic(units_out%energy, st%smear%e_fermi)

    write(message(4), '(7f12.6)')          &
      maxval(sb%kpoints%reduced%point(1,:)),           &
      maxval(sb%kpoints%reduced%point(2,:)),           &
      maxval(sb%kpoints%reduced%point(3,:)),           &
      maxval(sb%kpoints%reduced%point(1,:)/factor(1)), &
      maxval(sb%kpoints%reduced%point(2,:)/factor(2)), &
      maxval(sb%kpoints%reduced%point(3,:)/factor(3)), &
      units_from_atomic(units_out%energy, st%smear%e_fermi)

    call write_info(4, iunit)
    call io_close(iunit)

    ! now we write the same information so that it can be used 
    ! together with total-dos.dat
    iunit = io_open(trim(dir)//'/'//'total-dos-efermi.dat', action='write')    

    write(message(1), '(a)') '# Fermi energy in a format compatible with total-dos.dat'    

    ! this is the maximum that tdos can reach
    maxdos = sum(st%d%kweights) * st%nst

    write(message(2), '(4f12.6)') units_from_atomic(units_out%energy, st%smear%e_fermi), M_ZERO
    write(message(3), '(4f12.6)') units_from_atomic(units_out%energy, st%smear%e_fermi), maxdos

    call write_info(3, iunit)
    call io_close(iunit)

    call pop_sub('states.states_write_fermi_energy')
  end subroutine states_write_fermi_energy
  ! ---------------------------------------------------------


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
    integer :: sn, sn1, r, j, k
#endif

    call push_sub('states.states_distribute_nodes')

    ! Defaults.
    st%node(:)            = 0
    st%st_start           = 1
    st%st_end             = st%nst
    st%lnst               = st%nst
    st%parallel_in_states = .false.
    call mpi_grp_init(st%mpi_grp, -1)

#if defined(HAVE_MPI)
    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then
      st%parallel_in_states = .true.
      call mpi_grp_init(st%mpi_grp, mc%group_comm(P_STRATEGY_STATES))
      call mpi_grp_init(st%dom_st, mc%dom_st_comm)
      call multicomm_create_all_pairs(st%mpi_grp, st%ap)

     if(st%nst < st%mpi_grp%size) then
       message(1) = "Have more processors than necessary"
       write(message(2),'(i4,a,i4,a)') st%mpi_grp%size, " processors and ", st%nst, " states."
       call write_fatal(2)
     end if

     SAFE_ALLOCATE(st%st_range(1:2, 0:st%mpi_grp%size-1))
     SAFE_ALLOCATE(st%st_num(0:st%mpi_grp%size-1))

     call multicomm_divide_range(st%nst, st%mpi_grp%size, st%st_range(1, :), st%st_range(2, :), st%st_num)

     do k = 0, st%mpi_grp%size - 1
       write(message(1),'(a,i4,a,i7,a,i7)') &
            'Info: Nodes in states-group ', k, ' will manage states', st%st_range(1, k), " - ", st%st_range(2, k)
       call write_info(1)
       if(st%mpi_grp%rank .eq. k) then
         st%st_start = st%st_range(1, k)
         st%st_end   = st%st_range(2, k)
         st%lnst     = st%st_num(k)
       endif
     end do

     sn  = st%nst/st%mpi_grp%size
     sn1 = sn + 1
     r  = mod(st%nst, st%mpi_grp%size)
     do j = 1, r
       st%node((j-1)*sn1+1:j*sn1) = j - 1
     end do
     k = sn1*r
     call MPI_Barrier(st%mpi_grp%comm, mpi_err)
     do j = 1, st%mpi_grp%size - r
       st%node(k+(j-1)*sn+1:k+j*sn) = r + j - 1
     end do
   end if
#endif

    call pop_sub('states.states_distribute_nodes')
  end subroutine states_distribute_nodes


  ! ---------------------------------------------------------
  subroutine states_set_complex(st)
    type(states_t),    intent(inout) :: st

    st%priv%wfs_type = TYPE_CMPLX

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
  subroutine states_dump(st, iunit)
    type(states_t), intent(in) :: st
     integer,       intent(in) :: iunit

     call push_sub('states.states_dump')
     
     write(iunit, '(a20,1i10)')  'nst=                ', st%nst
     write(iunit, '(a20,1i10)')  'dim=                ', st%d%dim
     write(iunit, '(a20,1i10)')  'nik=                ', st%d%nik

     call pop_sub('states.states_dump')
  end subroutine states_dump


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
#if defined(HAVE_MPI)
    FLOAT, allocatable :: tmp_reduce(:)
    integer :: mpi_err
#endif

    call push_sub('states.states_calc_quantities')

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

        kpoint = M_ZERO
        kpoint = kpoints_get_point(der%mesh%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))

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
    if (present(density_laplacian)) then
      SAFE_DEALLOCATE_A(lwf_psi)
    end if

    if(.not. present(paramagnetic_current)) then
      SAFE_DEALLOCATE_P(jp)
    end if

    if(.not. present(kinetic_energy_density)) then
      SAFE_DEALLOCATE_P(tau)
    end if

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      call reduce_all(st%mpi_grp)
    end if
    if(st%d%kpt%parallel) then
      call reduce_all(st%d%kpt%mpi_grp)
    end if
#endif

    call pop_sub('states.states_calc_quantities')

#if defined(HAVE_MPI)
  contains 

    subroutine reduce_all(grp)
      type(mpi_grp_t), intent(in)  :: grp

      call push_sub('states.states_calc_quantities.reduce_all')

      SAFE_ALLOCATE(tmp_reduce(1:der%mesh%np))

      do is = 1, st%d%nspin
        if(associated(tau)) then
          call MPI_Allreduce(tau(1, is), tmp_reduce(1), der%mesh%np, MPI_FLOAT, MPI_SUM, grp%comm, mpi_err)
          tau(1:der%mesh%np, is) = tmp_reduce(1:der%mesh%np)       
        end if

        if(present(gi_kinetic_energy_density)) then
          call MPI_Allreduce(gi_kinetic_energy_density(1, is), tmp_reduce(1), der%mesh%np, MPI_FLOAT, MPI_SUM, grp%comm, mpi_err)
          gi_kinetic_energy_density(1:der%mesh%np, is) = tmp_reduce(1:der%mesh%np)       
        end if

        if (present(density_laplacian)) then
          call MPI_Allreduce(density_laplacian(1, is), tmp_reduce(1), der%mesh%np, MPI_FLOAT, MPI_SUM, grp%comm, mpi_err)
          density_laplacian(1:der%mesh%np, is) = tmp_reduce(1:der%mesh%np)       
        end if

        do i_dim = 1, der%mesh%sb%dim
          if(associated(jp)) then
            call MPI_Allreduce(jp(1, i_dim, is), tmp_reduce(1), der%mesh%np, MPI_FLOAT, MPI_SUM, grp%comm, mpi_err)
            jp(1:der%mesh%np, i_dim, is) = tmp_reduce(1:der%mesh%np)
          end if

          if(present(density_gradient)) then
            call MPI_Allreduce(density_gradient(1, i_dim, is), tmp_reduce(1), der%mesh%np, MPI_FLOAT, MPI_SUM, grp%comm, mpi_err)
            density_gradient(1:der%mesh%np, i_dim, is) = tmp_reduce(1:der%mesh%np)
          end if
        end do

      end do
      SAFE_DEALLOCATE_A(tmp_reduce)

      call pop_sub('states.states_calc_quantities.reduce_all')
    end subroutine reduce_all

#endif            
  end subroutine states_calc_quantities


  ! ---------------------------------------------------------
  function state_spin(m, f1) result(s)
    FLOAT, dimension(3) :: s
    type(mesh_t), intent(in) :: m
    CMPLX,  intent(in) :: f1(:, :)

    CMPLX :: z

    call push_sub('states.state_spin')

    z = zmf_dotp(m, f1(:, 1) , f1(:, 2))

    s(1) = M_TWO * z
    s(2) = M_TWO * aimag(z)
    s(3) = zmf_dotp(m, f1(:, 1), f1(:, 1)) - zmf_dotp(m, f1(:, 2), f1(:, 2))
    s = s * M_HALF ! spin is half the sigma matrix.

    call pop_sub('states.state_spin')
  end function state_spin


  ! ---------------------------------------------------------
  !> Reads the state stored in directory "dir", and finds out
  !! the kpoints, dim, and nst contained in it.
  ! ---------------------------------------------------------
  subroutine states_look(dir, mpi_grp, kpoints, dim, nst, ierr, only_occupied)
    character(len=*),  intent(in)    :: dir
    type(mpi_grp_t),   intent(in)    :: mpi_grp
    integer,           intent(out)   :: dim, ierr
    integer,           intent(inout) :: nst, kpoints
    logical, optional, intent(in)    :: only_occupied

    character(len=256) :: line
    character(len=12)  :: filename
    character(len=1)   :: char
    integer :: iunit, iunit2, err, i, ist, idim, ik
    integer, allocatable :: counter_kpoints(:), counter_nst(:)
    FLOAT :: occ, eigenval
    logical :: only_occupied_

    call push_sub('states.states_look')

    only_occupied_ = .false.
    if(present(only_occupied)) only_occupied_ = only_occupied
    if(only_occupied_) then
      ASSERT(nst>0 .and. kpoints>0)
      SAFE_ALLOCATE(counter_kpoints(1:nst))
      SAFE_ALLOCATE(counter_nst(1:kpoints))
    end if
    ierr = 0
    iunit  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(iunit < 0) then
      ierr = -1
    call pop_sub('states.states_look')
return
    end if
    iunit2 = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(iunit2 < 0) then
      call io_close(iunit, grp = mpi_grp)
      ierr = -1
    call pop_sub('states.states_look')
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
    if(only_occupied_) then
      counter_kpoints(:) = 0
      counter_nst(:) = 0
    end if

    do
      call iopar_read(mpi_grp, iunit, line, i)
      read(line, '(a)') char
      if(i.ne.0.or.char=='%') exit
      read(line, *) ik, char, ist, char, idim, char, filename
      if(idim == 2)    dim     = 2
      call iopar_read(mpi_grp, iunit2, line, err)
      read(line, *) occ, char, eigenval
      if(only_occupied_) then
        if(occ > M_EPSILON) then
          counter_kpoints(ist) = counter_kpoints(ist) + 1
          counter_nst(ik) = counter_nst(ik) + 1
          kpoints = maxval(counter_kpoints(:))
          nst = maxval(counter_nst(:))
        end if
      else
        if(ik > kpoints) kpoints = ik
        if(ist>nst)      nst     = ist
      end if
    end do

    call io_close(iunit, grp = mpi_grp)
    call io_close(iunit2, grp = mpi_grp)
    if(only_occupied_) then
      SAFE_DEALLOCATE_A(counter_kpoints)
      SAFE_DEALLOCATE_A(counter_nst)
    end if

    call pop_sub('states.states_look')
  end subroutine states_look


  ! ---------------------------------------------------------
  logical function state_is_local(st, ist)
    type(states_t), intent(in) :: st
    integer,        intent(in) :: ist

    call push_sub('states.state_is_local')

    state_is_local = ist.ge.st%st_start.and.ist.le.st%st_end

    call pop_sub('states.state_is_local')
  end function state_is_local


  ! ---------------------------------------------------------
  subroutine states_freeze_orbitals(st, gr, mc, n)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(multicomm_t), intent(in)    :: mc
    integer,           intent(in)    :: n

    integer :: ist, ik
    type(states_t) :: staux

    call push_sub('states.states_freeze_orbitals')

    if(n >= st%nst) then
      write(message(1),'(a)') 'Attempting to freeze a number of orbitals which is larger or equal to'
      write(message(2),'(a)') 'the total number. The program has to stop.'
      call write_fatal(2)
    end if

    if(.not.associated(st%frozen_rho)) then
      SAFE_ALLOCATE(st%frozen_rho(1:gr%mesh%np, 1:st%d%dim))
      st%frozen_rho = M_ZERO
    end if

    st%frozen_rho = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        if(ist > n) cycle
        call states_dens_accumulate(st, gr, ist, ik, st%frozen_rho)
      end do
    end do
    call states_dens_reduce(st, gr, st%frozen_rho)

    call states_copy(staux, st)

    st%nst = st%nst - n

    call states_deallocate_wfns(st)
    call states_distribute_nodes(st, mc)
    call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX)

#if defined(HAVE_MPI) 

    if(staux%parallel_in_states) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = staux%st_start, staux%st_end
          if(ist <= n) cycle
          if(.not.state_is_local(st, ist-n)) then
            call mpi_send(staux%zpsi(1, 1, ist, ik), gr%mesh%np_part*st%d%dim, MPI_CMPLX, staux%node(ist), &
              ist, st%mpi_grp%comm, mpi_err)

            call mpi_recv(st%zpsi(1, 1, ist-n, ik), gr%mesh%np_part*st%d%dim, MPI_CMPLX, st%node(ist-n), &
              ist, st%mpi_grp%comm, mpi_err)
          else
            st%zpsi(:, :, ist-n, ik) = staux%zpsi(:, :, ist, ik)
          end if
   
        end do
      end do
   else
     do ik = st%d%kpt%start, st%d%kpt%end
       do ist = staux%st_start, staux%st_end
         if(ist <= n) cycle
         st%zpsi(:, :, ist-n, ik) = staux%zpsi(:, :, ist, ik)
       end do
     end do
   end if

#else

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%zpsi(:, :, ist, ik) = staux%zpsi(:, :, n + ist, ik)
      end do
    end do

#endif

    SAFE_DEALLOCATE_P(st%occ)
    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)
    st%occ      = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%occ(ist, ik) = staux%occ(n+ist, ik)
        st%eigenval(ist, ik) = staux%eigenval(n+ist, ik)
      end do
    end do

    call states_end(staux)
    call pop_sub('states.states_freeze_orbitals')
  end subroutine states_freeze_orbitals


  ! ---------------------------------------------------------
  !> this routine calculates the total electronic density,
  !! which is the sum of the part coming from the orbitals, the
  !! non-linear core corrections and the frozen orbitals
  subroutine states_total_density(st, mesh, rho)
    type(states_t), intent(in)  :: st
    type(mesh_t),   intent(in)  :: mesh
    FLOAT,          intent(out) :: rho(:,:)

    integer :: is, ip

    call push_sub('states.states_total_density')

    forall(is = 1:st%d%nspin, ip = 1:mesh%np)
      rho(ip, is) = st%rho(ip, is)
    end forall

    if(associated(st%rho_core)) then
      forall(is = 1:st%d%spin_channels, ip = 1:mesh%np)
        rho(ip, is) = rho(ip, is) + st%rho_core(ip)/st%d%nspin
      end forall
    end if

    ! Add, if it exists, the frozen density from the inner orbitals.
    if(associated(st%frozen_rho)) then
      forall(is = 1:st%d%spin_channels, ip = 1:mesh%np)
        rho(ip, is) = rho(ip, is) + st%frozen_rho(ip, is)
      end forall
    end if

    call pop_sub('states.states_total_density')
  end subroutine states_total_density


  ! ---------------------------------------------------------
  real(8) function states_wfns_memory(st, mesh) result(memory)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: mesh

    call push_sub('states.states_wfns_memory')
    memory = 0.0_8

    ! orbitals
    memory = memory + REAL_PRECISION*dble(mesh%np_part_global)*st%d%dim*dble(st%nst)*st%d%kpt%nglobal

    call pop_sub('states.states_wfns_memory')
  end function states_wfns_memory


  ! ---------------------------------------------------------
  !> initialize the self-energy of the leads
  subroutine states_init_self_energy(st, gr, nspin, d_ispin, lead)
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: nspin
    integer,             intent(in)    :: d_ispin
    type(lead_t),        intent(in)    :: lead(:) ! Diagonal and off-diagonal block of the lead Hamiltonian.

    character(len=2)      :: spin
    character(len=256)    :: fmt, fname, fname_real, fname_imag
    FLOAT                 :: energy
    integer  :: np, ik, ist, il, ispin, s1, s2, k1, k2
    integer  :: green_real, green_imag, irow

    call push_sub('states.states_init_self_energy')

    ! Calculate self-energy of the leads.
    ! FIXME: For spinors, this calculation is almost certainly wrong.
    ASSERT(st%ob_nst == st%nst)
    ASSERT(st%ob_d%nik == st%d%nik)

    s1 = st%st_start
    s2 = st%st_end
    k1 = st%d%kpt%start
    k2 = st%d%kpt%end

    do il = 1, NLEADS
      np = gr%intf(il)%np_intf
      SAFE_ALLOCATE(st%ob_lead(il)%self_energy(1:np, 1:np, 1:nspin, s1:s2, k1:k2))
    end do
    call messages_print_stress(stdout, "Lead self-energy")
    message(1) = ' st#     k#  Spin      Lead     Energy'
    call write_info(1)
#ifdef HAVE_MPI 
    ! wait for all processors to finish 
    if(st%d%kpt%parallel) then 
      call MPI_Barrier(st%d%kpt%mpi_grp%comm, mpi_err) 
    end if
#endif
    do ik = k1, k2
      do ist = s1, s2
        energy = st%ob_eigenval(ist, ik)
        do il = 1, NLEADS
          np = gr%intf(il)%np_intf
          do ispin = 1, nspin
            select case(d_ispin)
            case(UNPOLARIZED)
              spin = '--'
            case(SPIN_POLARIZED)
              if(is_spin_up(ik)) then
                spin = 'up'
              else
                spin = 'dn'
              end if
              ! This is nonsense, but at least all indices are present.
            case(SPINORS)
              if(ispin.eq.1) then
                spin = 'up'
              else
                spin = 'dn'
              end if
            end select
            write(message(1), '(i4,3x,i4,3x,a2,5x,a6,1x,f12.6)') ist, ik, &
              trim(spin), trim(LEAD_NAME(il)), energy
            call write_info(1)

            call lead_self_energy(energy, lead(il)%h_diag(:, :, ispin), lead(il)%h_offdiag(:, :), &
              gr%intf(il), st%ob_lead(il)%self_energy(:, :, ispin, ist, ik), .true.)

            ! Write the entire self-energy to a file.
            if(in_debug_mode) then
              call io_mkdir('debug/open_boundaries')
              write(fname_real, '(3a,i4.4,a,i3.3,a,i1.1,a)') 'debug/open_boundaries/self-energy-', &
                trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', ispin, '.real'
              write(fname_imag, '(3a,i4.4,a,i3.3,a,i1.1,a)') 'debug/open_boundaries/self-energy-', &
                trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', ispin, '.imag'
              green_real = io_open(fname_real, action='write', grp=st%d%kpt%mpi_grp, is_tmp=.false.)
              green_imag = io_open(fname_imag, action='write', grp=st%d%kpt%mpi_grp, is_tmp=.false.)

              write(fmt, '(a,i6,a)') '(', np, 'e24.16)'
              do irow = 1, np
                write(green_real, fmt) real(st%ob_lead(il)%self_energy(irow, :, ispin, ist, ik))
                write(green_imag, fmt) aimag(st%ob_lead(il)%self_energy(irow, :, ispin, ist, ik))
              end do
              call io_close(green_real)
              call io_close(green_imag)
            end if
          end do
        end do
      end do
    end do
    call messages_print_stress(stdout)

    call pop_sub('states.states_init_self_energy')
  end subroutine states_init_self_energy


  ! ---------------------------------------------------------
  ! write H_(C,apha)*Psi_(alpha) without using Psi_(alpha)
  subroutine states_write_proj_lead_wf(sb, dir, intf, st)
    type(simul_box_t), intent(in) :: sb
    character(len=*),  intent(in) :: dir   ! directory
    type(interface_t), intent(in) :: intf(:)
    type(states_t),    intent(in) :: st

    integer :: ik, ist, idim, il, np, ip, iunit
    CMPLX, allocatable :: psi(:, :), phi(:, :), hpsi(:, :), self_energy(:, :)
    character(len=256) :: fname
    FLOAT :: kpoint(1:3)

    call push_sub('states.write_proj_lead_wf')

    np = maxval(intf(1:NLEADS)%np_intf)

    SAFE_ALLOCATE(psi(1:np, 1:st%d%dim))
    SAFE_ALLOCATE(phi(1:np, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:np, 1:st%d%dim))
    SAFE_ALLOCATE(self_energy(1:np, 1:np))

    hpsi(:, :) = M_z0

    if(mpi_grp_is_root(mpi_world)) call io_mkdir(dir, is_tmp=.true.)
#ifdef HAVE_MPI
    if(st%parallel_in_states.or.st%d%kpt%parallel) call MPI_Barrier(st%d%kpt%mpi_grp%comm, mpi_err)
#endif

    do ik = st%d%kpt%start, st%d%kpt%end

      kpoint = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik))

      do ist = st%st_start, st%st_end
        do il = 1, NLEADS
          np = intf(il)%np_intf
          forall(ip = 1:np)
            psi(ip, :) = st%zpsi(intf(il)%index(ip), :, ist, ik)
            phi(ip, :) = st%zphi(intf(il)%index(ip), :, ist, ik)
          end forall
          do idim = 1, st%d%dim
            self_energy(1:np, 1:np) = st%ob_lead(il)%self_energy(1:np, 1:np, idim, ist, ik)
            call lalg_gemv(np, np, M_z1, self_energy(1:np, 1:np), psi(1:np, idim), M_z0, hpsi(1:np, idim))
            if((il.eq.LEFT).and.(-kpoint(1).gt.M_ZERO) .or. (il.eq.RIGHT).and.(-kpoint(1).lt.M_ZERO)) then
              ! add the reflecting part
              self_energy(1:np, 1:np) = transpose(conjg(self_energy(1:np, 1:np))) - self_energy(1:np, 1:np)
              call lalg_gemv(np, np, M_z1, self_energy(1:np, 1:np), phi(1:np, idim), M_z1, hpsi(1:np, idim))
            end if
            write(fname, '(3a,i4.4,a,i3.3,a,i1.1)') 'src0-', trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', idim
            iunit = io_open(trim(dir)//trim(fname), action='write', form='unformatted', is_tmp=.true.)
            if(iunit.lt.0) then
              message(1) = 'Cannot write term for source term to file.'
              call write_warning(1)
              call io_close(iunit)
    call pop_sub('states.write_proj_lead_wf')
return
            end if
            ! Write parameters.
            write(iunit) np
            ! Write matrix.
            write(iunit) hpsi(1:np, idim)
            call io_close(iunit)
          end do
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(phi)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(self_energy)

    call pop_sub('states.write_proj_lead_wf')
  end subroutine states_write_proj_lead_wf

  ! ---------------------------------------------------------
  subroutine states_read_proj_lead_wf(dir, intf, st, src0)
    character(len=*),  intent(in)    :: dir   ! directory
    type(interface_t), intent(in)    :: intf
    type(states_t),    intent(in)    :: st
    CMPLX,             intent(inout) :: src0(:, : ,:, :)

    integer :: ik, ist, idim, np, ip, iunit
    character(len=256) :: fname

    call push_sub('states.states_read_proj_lead_wf')

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          ! Try to open file.
          write(fname, '(3a,i4.4,a,i3.3,a,i1.1)') 'src0-', trim(LEAD_NAME(intf%il)), '-', ist, '-', ik, '-', idim
          iunit = io_open(trim(dir)//trim(fname), action='read', status='old', die=.false., is_tmp=.true., form='unformatted')
          if(iunit.lt.0) then ! no file found
            message(1) = 'Cannot read src(0) from file.'
            call write_fatal(1)
          end if

          ! Now read the data.
          read(iunit) np

          if(np.ne.size(src0, 1)) then
            message(1) = 'Size mismatch! Cannot read src(0) from file.'
            call write_fatal(1)
          end if

          ! because we use a sliced array we have to remap the index
          read(iunit) src0(1:np, idim, ist-st%st_start+1, ik-st%d%kpt%start+1)

          call io_close(iunit)
        end do
      end do
    end do

    call pop_sub('states.states_read_proj_lead_wf')
  end subroutine states_read_proj_lead_wf


end module states_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
