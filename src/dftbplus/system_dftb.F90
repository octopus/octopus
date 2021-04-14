!! Copyright (C) 2020 F. Bonafé, H. Appel
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

module system_dftb_oct_m
  use algorithm_oct_m
  use clock_oct_m
#ifdef HAVE_DFTBPLUS
  use dftbplus
#endif
  use geometry_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use iso_c_binding
  use lasers_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m
  use quantity_oct_m
  use space_oct_m
  use species_oct_m
  use system_oct_m
  use tdfunction_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::           &
    system_dftb_t,    &
    system_dftb_init

   type, extends(system_t) :: system_dftb_t
    integer :: n_atom
    FLOAT, allocatable :: coords(:,:), gradients(:,:)
    FLOAT, allocatable :: acc(:,:)
    FLOAT, allocatable :: tot_force(:,:)
    FLOAT, allocatable :: vel(:,:)
    FLOAT, allocatable :: prev_tot_force(:,:) !< Used for the SCF convergence criterium
    integer, allocatable :: species(:)
    integer              :: dynamics
    FLOAT, allocatable :: mass(:)
    FLOAT, allocatable :: atom_charges(:,:) !< shape is (n_atoms, n_spin)
    character(len=LABEL_LEN), allocatable  :: labels(:)
    FLOAT, allocatable :: prev_acc(:,:,:) !< A storage of the prior times.
    FLOAT :: scc_tolerance
    type(geometry_t), pointer :: geo => NULL()
    type(c_ptr) :: output_handle(2)
    type(ion_dynamics_t) :: ions
    class(lasers_t), pointer :: ext_lasers => null()
    logical :: laser_field
    FLOAT :: field(3)
    FLOAT :: energy
#ifdef HAVE_DFTBPLUS
    type(TDftbPlus) :: dftbp
#endif
  contains
    procedure :: init_interaction => system_dftb_init_interaction
    procedure :: initial_conditions => system_dftb_initial_conditions
    procedure :: do_td_operation => system_dftb_do_td
    procedure :: iteration_info => system_dftb_iteration_info
    procedure :: output_start => system_dftb_output_start
    procedure :: output_write => system_dftb_output_write
    procedure :: output_finish => system_dftb_output_finish
    procedure :: is_tolerance_reached => system_dftb_is_tolerance_reached
    procedure :: update_quantity => system_dftb_update_quantity
    procedure :: update_exposed_quantity => system_dftb_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => system_dftb_copy_quantities_to_interaction
    procedure :: update_interactions_start => system_dftb_update_interactions_start
    procedure :: update_interactions_finish => system_dftb_update_interactions_finish
    final :: system_dftb_finalize
  end type system_dftb_t

  interface system_dftb_t
    procedure system_dftb_constructor
  end interface system_dftb_t

  !> Parameters.
  integer, parameter :: &
    EHRENFEST = 1,   &
    BO        = 2

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function system_dftb_constructor(namespace) result(sys)
    class(system_dftb_t), pointer    :: sys
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(system_dftb_constructor)

    SAFE_ALLOCATE(sys)

    call system_dftb_init(sys, namespace)

    POP_SUB(system_dftb_constructor)
  end function system_dftb_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine system_dftb_init(this, namespace)
    class(system_dftb_t), target, intent(inout) :: this
    type(namespace_t),            intent(in)    :: namespace

    integer :: ii, jj, ispec, il, ierr
    character(len=MAX_PATH_LEN) :: slako_dir
    character(len=1), allocatable  :: max_ang_mom(:)
    character(len=LABEL_LEN) :: this_max_ang_mom, this_label
    integer :: n_maxang_block
    type(block_t) :: blk
    character(len=200) :: envelope_expression, phase_expression
    FLOAT :: omega0, initial_temp
    CMPLX :: pol(MAX_DIM)
    type(tdf_t) :: ff, phi

#ifdef HAVE_DFTBPLUS
    type(TDftbPlusInput) :: input
    type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pAnalysis
    type(fnode), pointer :: pParserOpts
#endif
#ifdef HAVE_DFTBPLUS_DEVEL
    type(fnode), pointer :: pElecDyn, pPerturb, pLaser
#endif

    PUSH_SUB(system_dftb_init)

    this%namespace = namespace

    call messages_print_stress(stdout, "DFTB+ System", namespace=namespace)

    call space_init(this%space, namespace)
    if(this%space%periodic_dim > 0) then
      call messages_not_implemented('DFTB+ for periodic systems')
    end if
    if(this%space%dim /= 3) then
      call messages_not_implemented('DFTB+ for Dimensions /= 3')
    end if

    this%geo => geometry_t(namespace, this%space)
    this%n_atom = this%geo%natoms
    SAFE_ALLOCATE(this%coords(3, this%n_atom))
    SAFE_ALLOCATE(this%acc(3, this%n_atom))
    SAFE_ALLOCATE(this%vel(3, this%n_atom))
    SAFE_ALLOCATE(this%tot_force(3, this%n_atom))
    SAFE_ALLOCATE(this%prev_tot_force(3, this%n_atom))
    SAFE_ALLOCATE(this%gradients(3, this%n_atom))
    SAFE_ALLOCATE(this%species(this%n_atom))
    SAFE_ALLOCATE(this%mass(this%n_atom))
    SAFE_ALLOCATE(this%atom_charges(this%n_atom, 1))
    SAFE_ALLOCATE(this%labels(this%geo%nspecies))
    SAFE_ALLOCATE(max_ang_mom(this%geo%nspecies))

    ispec = 1
    this%species(1) = 1
    this%labels(1) = trim(this%geo%atom(1)%label)

    do ii = 1, this%n_atom
      this%coords(1:3,ii) = this%geo%atom(ii)%x(1:3)
      ! mass is read from the default pseudopotential files
      this%mass(ii) = species_mass(this%geo%atom(ii)%species)
      if ((ii > 1) .and. .not. (any(this%labels(1:ispec) == this%geo%atom(ii)%label))) then
        ispec = ispec + 1
        this%labels(ispec) = trim(this%geo%atom(ii)%label)
      end if
      do jj = 1, ispec
        if (trim(this%geo%atom(ii)%label) == trim(this%labels(jj))) then
          this%species(ii) = jj
        end if
      end do
    end do
    this%vel = M_ZERO
    this%tot_force = M_ZERO

    !%Variable MaxAngularMomentum
    !%Type block
    !%Section DFTBPlusInterface
    !%Description
    !% Specifies the highest angular momentum for each atom type. All orbitals up
    !% to that angular momentum will be included in the calculation.
    !% Possible values for the angular momenta are s, p, d, f.
    !% These are examples:
    !%
    !% <tt>%MaxAngularMomentum
    !% <br>&nbsp;&nbsp;'O'   | 'p'
    !% <br>&nbsp;&nbsp;'H'   | 's'
    !% <br>%</tt>
    !%End
    n_maxang_block = 0
    if(parse_block(namespace, 'MaxAngularMomentum', blk) == 0) then
      n_maxang_block = parse_block_n(blk)
      if (n_maxang_block /= this%geo%nspecies) then
        call messages_input_error(namespace, "MaxAngularMomentum", "Wrong number of species.")
      end if

      do ii = 1, n_maxang_block
        call parse_block_string(blk, ii-1, 0, this_label)
        call parse_block_string(blk, ii-1, 1, this_max_ang_mom)
        if (any(["s","p","d","f"] == trim(this_max_ang_mom))) then
          call messages_input_error(namespace, "MaxAngularMomentum", "Wrong maximum angular momentum for element"//trim(this_label))
        end if
        do jj = 1, this%geo%nspecies
          if (trim(adjustl(this_label)) == trim(adjustl(this%labels(jj)))) then
            max_ang_mom(jj) = trim(adjustl(this_max_ang_mom))
          end if
        end do
      end do
    end if
    call parse_block_end(blk)

    !%Variable SlakoDir
    !%Type string
    !%Default "./"
    !%Section Execution::IO
    !%Description
    !% Folder containing the Slako files
    !%End
    call parse_variable(namespace, 'SlakoDir', './', slako_dir)


    ! Dynamics variables

    call ion_dynamics_init(this%ions, namespace, this%geo)

    call parse_variable(namespace, 'TDDynamics', BO, this%dynamics)
    call messages_print_var_option(stdout, 'TDDynamics', this%dynamics)
    if (this%dynamics == BO) then
      call ion_dynamics_unfreeze(this%ions)
    end if

    !%Variable InitialIonicTemperature
    !%Type float
    !%Default 0.0
    !%Section DFTBPlusInterface
    !%Description
    !% If this variable is present, the ions will have initial velocities
    !% velocities to the atoms following a Boltzmann distribution with
    !% this temperature (in Kelvin). Used only if <tt>TDDynamics = Ehrenfest</tt>
    !% and  <tt>MoveIons = yes</tt>.
    !%End
    call parse_variable(namespace, 'InitialIonicTemperature', M_zero, initial_temp, unit = unit_kelvin)

    allocate(this%ext_lasers)
    this%ext_lasers%no_lasers = 0
    if(parse_block(namespace, 'TDExternalFields', blk) == 0) then
      this%laser_field = .true.
      ! No call to safe_deallocate macro here, as it gives an ICE with gfortran
      this%ext_lasers%no_lasers = parse_block_n(blk)
      SAFE_ALLOCATE(this%ext_lasers%lasers(1:this%ext_lasers%no_lasers))

      do il = 1, this%ext_lasers%no_lasers

        pol(1:MAX_DIM) = M_z0
        call parse_block_cmplx(blk, il-1, 0, pol(1))
        call parse_block_cmplx(blk, il-1, 1, pol(2))
        call parse_block_cmplx(blk, il-1, 2, pol(3))
        call parse_block_float(blk, il-1, 3, omega0)
        omega0 = units_to_atomic(units_inp%energy, omega0)

        call laser_set_frequency(this%ext_lasers%lasers(il), omega0)
        pol(1:3) = pol(1:3)/sqrt(sum(abs(pol(1:3))**2))
        call laser_set_polarization(this%ext_lasers%lasers(il), pol)

        call parse_block_string(blk, il-1, 4, envelope_expression)
        call tdf_read(ff, namespace, trim(envelope_expression), ierr)
        call laser_set_f(this%ext_lasers%lasers(il), ff)

        ! Check if there is a phase.
        if(parse_block_cols(blk, il-1) > 5) then
          call parse_block_string(blk, il-1, 5, phase_expression)
          call tdf_read(phi, namespace, trim(phase_expression), ierr)
          call laser_set_phi(this%ext_lasers%lasers(il), phi)
          if (ierr /= 0) then
            write(message(1),'(3A)') 'Error in the "', trim(envelope_expression), '" field defined in the TDExternalFields block:'
            write(message(2),'(3A)') 'Time-dependent phase function "', trim(phase_expression), '" not found.'
            call messages_warning(2, namespace=namespace)
          end if
        else
          call laser_set_empty_phi(this%ext_lasers%lasers(il))
        end if
      end do

      call parse_block_end(blk)
    else
      this%laser_field = .false.
    end if

#ifdef HAVE_DFTBPLUS
    call TDftbPlus_init(this%dftbp, mpicomm=mpi_world%comm)

    call this%dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call setChild(pRoot, "Geometry", pGeo)
    call setChildValue(pGeo, "Periodic", .false.)
    call setChildValue(pGeo, "TypeNames", this%labels(1:this%geo%nspecies))
    call setChildValue(pGeo, "TypesAndCoordinates", reshape(this%species, [1, size(this%species)]), this%coords)
    call setChild(pRoot, "Hamiltonian", pHam)
    call setChild(pHam, "Dftb", pDftb)
    call setChildValue(pDftb, "Scc", .true.)

    !%Variable SccTolerance
    !%Type float
    !%Section DFTBPlusInterface
    !%Description
    !% Self-consistent-charges convergence tolerance. Once this
    !% tolerance has been achieved the SCC cycle will stop.
    !%End
    call parse_variable(namespace, 'SccTolerance', CNST(1e-9), this%scc_tolerance)
    call messages_print_var_value(stdout, 'SccTolerance', this%scc_tolerance)
    call setChildValue(pDftb, "SccTolerance", this%scc_tolerance)

    ! sub-block inside hamiltonian for the maximum angular momenta
    call setChild(pDftb, "MaxAngularMomentum", pMaxAng)
    ! explicitly set the maximum angular momenta for the species
    do ii = 1, this%geo%nspecies
      call setChildValue(pMaxAng, this%labels(ii), max_ang_mom(ii))
    end do

    ! get the SK data
    ! You should provide the skfiles as found in the external/slakos/origin/mio-1-1/ folder. These can
    ! be downloaded with the utils/get_opt_externals script
    call setChild(pDftb, "SlaterKosterFiles", pSlakos)
    call setChild(pSlakos, "Type2FileNames", pType2Files)
    call setChildValue(pType2Files, "Prefix", slako_dir)
    call setChildValue(pType2Files, "Separator", "-")
    call setChildValue(pType2Files, "Suffix", ".skf")

    !  set up analysis options
    call setChild(pRoot, "Analysis", pAnalysis)
    call setChildValue(pAnalysis, "CalculateForces", .true.)

    call setChild(pRoot, "ParserOptions", pParserOpts)
    call setChildValue(pParserOpts, "ParserVersion", 5)
#endif

    if (this%dynamics == EHRENFEST) then
#ifdef HAVE_DFTBPLUS_DEVEL
      call setChild(pRoot, "ElectronDynamics", pElecDyn)
      call setChildValue(pElecDyn, "IonDynamics", ion_dynamics_ions_move(this%ions))
      if (ion_dynamics_ions_move(this%ions)) then
        call setChildValue(pElecDyn, "InitialTemperature", initial_temp)
      end if

      ! initialize with wrong arguments for the moment, will be overriden later
      call setChildValue(pElecDyn, "Steps", 1)
      call setChildValue(pElecDyn, "TimeStep", CNST(1.0))
      call setChild(pElecDyn, "Perturbation", pPerturb)
      if (this%laser_field) then
        call setChild(pPerturb, "Laser", pLaser)
        call setChildValue(pLaser, "PolarizationDirection", [ CNST(1.0) , CNST(0.0) , CNST(0.0) ])
        call setChildValue(pLaser, "LaserEnergy", CNST(1.0))
        call setChildValue(pElecDyn, "FieldStrength", CNST(1.0))
      else
        call setChild(pPerturb, "None", pLaser)
      end if
#else
      message(1) = "DFTB Ehrenfest dynamics enabled only in DFTB development library"
      call messages_fatal(1)
#endif
    end if

#ifdef HAVE_DFTBPLUS
    message(1) = 'Input tree in HSD format:'
    call messages_info(1)
    call dumpHsd(input%hsdTree, stdout)

    ! initialise the DFTB+ calculator
    call this%dftbp%setupCalculator(input)
    call this%dftbp%setGeometry(this%coords)
#endif

    POP_SUB(system_dftb_init)
  end subroutine system_dftb_init

  ! ---------------------------------------------------------
  subroutine system_dftb_init_interaction(this, interaction)
    class(system_dftb_t), target, intent(inout) :: this
    class(interaction_t),                intent(inout) :: interaction

    PUSH_SUB(system_dftb_init_interaction)

    select type (interaction)
    class default
      message(1) = "Trying to initialize an unsupported interaction by DFTB+."
      call messages_fatal(1)
    end select

    POP_SUB(system_dftb_init_interaction)
  end subroutine system_dftb_init_interaction

  ! ---------------------------------------------------------
  subroutine system_dftb_initial_conditions(this, from_scratch)
    class(system_dftb_t), intent(inout) :: this
    logical,                 intent(in)    :: from_scratch

    PUSH_SUB(system_dftb_initial_conditions)

    select case (this%dynamics)
    case (BO)
#ifdef HAVE_DFTBPLUS
      call this%dftbp%getGradients(this%gradients)
      this%tot_force = -this%gradients
#endif
    case (EHRENFEST)
#ifdef HAVE_DFTBPLUS_DEVEL
      call this%dftbp%getEnergy(this%energy)
      call this%dftbp%initializeTimeProp(this%prop%dt, this%laser_field)
#endif
    end select

    POP_SUB(system_dftb_initial_conditions)
  end subroutine system_dftb_initial_conditions

  ! ---------------------------------------------------------
  subroutine system_dftb_do_td(this, operation)
    class(system_dftb_t),    intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    integer :: ii, jj, il
    type(tdf_t) :: ff, phi
    CMPLX :: amp, pol(MAX_DIM)
    FLOAT :: time, omega

    PUSH_SUB(system_dftb_do_td)

    select case(this%dynamics)
    ! Born-Oppenheimer dynamics
    case (BO)

      select case (operation%id)
      case (SKIP)
        ! Do nothing
      case (STORE_CURRENT_STATUS)
        ! Do nothing

      case (VERLET_START)
        SAFE_ALLOCATE(this%prev_acc(1:this%space%dim, this%n_atom, 1))
        do jj = 1, this%n_atom
          this%acc(1:this%space%dim, jj) = this%tot_force(1:this%space%dim, jj) / this%mass(jj)
        end do

      case (VERLET_FINISH)
        SAFE_DEALLOCATE_A(this%prev_acc)

      case (VERLET_UPDATE_POS)
        do jj = 1, this%n_atom
          this%coords(1:this%space%dim, jj) = this%coords(1:this%space%dim, jj) + this%prop%dt * this%vel(1:this%space%dim, jj) &
               + M_HALF * this%prop%dt**2 * this%acc(1:this%space%dim, jj)
        end do
        this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name=QUANTITY_LABEL(POSITION), & 
                                            clock = this%quantities(POSITION)%clock, action="tick") )

      case (VERLET_COMPUTE_ACC)
        do ii = size(this%prev_acc, dim=3) - 1, 1, -1
          this%prev_acc(1:this%space%dim, 1:this%n_atom, ii + 1) = this%prev_acc(1:this%space%dim, 1:this%n_atom, ii)
        end do
        this%prev_acc(1:this%space%dim, 1:this%n_atom, 1) = this%acc(1:this%space%dim, 1:this%n_atom)
#ifdef HAVE_DFTBPLUS
        call this%dftbp%setGeometry(this%coords)
        call this%dftbp%getGradients(this%gradients)
        this%tot_force = -this%gradients
#endif
        do jj = 1, this%n_atom
          this%acc(1:this%space%dim, jj) = this%tot_force(1:this%space%dim, jj) / this%mass(jj)
        end do

      case (VERLET_COMPUTE_VEL)
        this%vel(1:this%space%dim, 1:this%n_atom) = this%vel(1:this%space%dim, 1:this%n_atom) &
            + M_HALF * this%prop%dt * (this%prev_acc(1:this%space%dim, 1:this%n_atom, 1) + &
              this%acc(1:this%space%dim, 1:this%n_atom))
        this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name=QUANTITY_LABEL(VELOCITY), & 
                                            clock = this%quantities(VELOCITY)%clock, action="tick") )


      case default
        message(1) = "Unsupported TD operation."
        call messages_fatal(1, namespace=this%namespace)
      end select

    ! Ehrenfest dynamics
    case (EHRENFEST)

      select case (operation%id)
      case (SKIP)
        ! Do nothing
      case (STORE_CURRENT_STATUS)
        ! Do nothing
      case (VERLET_START)
        !Do nothing
      case (VERLET_FINISH)
        !Do nothing
      case (VERLET_UPDATE_POS)
        this%field = M_zero
        time = this%clock%time()
        do il = 1, this%ext_lasers%no_lasers
          ! get properties of laser
          call laser_get_f(this%ext_lasers%lasers(il), ff)
          call laser_get_phi(this%ext_lasers%lasers(il), phi)
          omega = laser_carrier_frequency(this%ext_lasers%lasers(il))
          pol = laser_polarization(this%ext_lasers%lasers(il))
          ! calculate electric field from laser
          amp = tdf(ff, time) * exp(M_zI * (omega*time + tdf(phi, time)))
          this%field(1:3) = this%field(1:3) + TOFLOAT(amp*pol(1:3))
        end do
#ifdef HAVE_DFTBPLUS_DEVEL
        call this%dftbp%setTdElectricField(this%field)
        call this%dftbp%doOneTdStep(this%clock%get_tick(), atomNetCharges=this%atom_charges, coord=this%coords,&
             force=this%tot_force, energy=this%energy)
#endif
      case (VERLET_COMPUTE_ACC)
        !Do nothing
      case (VERLET_COMPUTE_VEL)
        !Do nothing
      case default
        message(1) = "Unsupported TD operation."
        call messages_fatal(1, namespace=this%namespace)
      end select

    end select
   POP_SUB(system_dftb_do_td)
  end subroutine system_dftb_do_td

  ! ---------------------------------------------------------
  logical function system_dftb_is_tolerance_reached(this, tol) result(converged)
    class(system_dftb_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(system_dftb_is_tolerance_reached)

    ! this routine is never called at present, no reason to be here
    ASSERT(.false.)
    converged = .false.

    POP_SUB(system_dftb_is_tolerance_reached)
  end function system_dftb_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine system_dftb_iteration_info(this)
    class(system_dftb_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(system_dftb_iteration_info)

    write(message(1),'(2X,A,1X,A)') "DFTB+ System:", trim(this%namespace%get())

    write(fmt,'("(4X,A,1X,",I2,"e14.6)")') this%space%dim
    write(message(2),fmt) "Coordinates: ", (this%coords(idir, 1)/P_Ang, idir = 1, this%space%dim)
    write(message(3),fmt) "Velocity:    ", (this%vel(idir, 1), idir = 1, this%space%dim)
    write(message(4),fmt) "Force:       ", (this%tot_force(idir, 1), idir = 1, this%space%dim)
    write(message(5),'(4x,A,I8.7)')  'Clock tick:      ', this%clock%get_tick()
    write(message(6),'(4x,A,e14.6)') 'Simulation time: ', this%clock%time()
    call messages_info(6)

    POP_SUB(system_dftb_iteration_info)
  end subroutine system_dftb_iteration_info

  ! ---------------------------------------------------------
  subroutine system_dftb_output_start(this)
    class(system_dftb_t), intent(inout) :: this

    PUSH_SUB(system_dftb_output_start)

    ! Create output handle
    call io_mkdir('td.general', this%namespace)
    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_init(this%output_handle(1), 0, this%prop%dt, trim(io_workpath("td.general/coordinates", this%namespace)))
      call write_iter_init(this%output_handle(2), 0, this%prop%dt, trim(io_workpath("td.general/forces", this%namespace)))
    end if

    ! Output info for first iteration
    call this%output_write()

    POP_SUB(system_dftb_output_start)
  end subroutine system_dftb_output_start

  ! ---------------------------------------------------------
  subroutine system_dftb_output_finish(this)
    class(system_dftb_t), intent(inout) :: this

    PUSH_SUB(system_dftb_output_finish)

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle(1))
      call write_iter_end(this%output_handle(2))
    end if

    POP_SUB(system_dftb_output_finish)
  end subroutine system_dftb_output_finish

  ! ---------------------------------------------------------
  subroutine system_dftb_output_write(this)
    class(system_dftb_t), intent(inout) :: this

    integer :: idir, iat, iout
    character(len=50) :: aux
    character(1) :: out_label(2)
    FLOAT :: tmp(MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(system_dftb_output_write)

    out_label(1) = "x"
    out_label(2) = "f"

    if (this%clock%get_tick() == 0) then
       ! header
      do iout = 1, 2
        call write_iter_clear(this%output_handle(iout))
        call write_iter_string(this%output_handle(iout),'#####################################################################')
        call write_iter_nl(this%output_handle(iout))
        call write_iter_string(this%output_handle(iout),'# HEADER')
        call write_iter_nl(this%output_handle(iout))

        ! first line: column names
        call write_iter_header_start(this%output_handle(iout))

        do iat = 1, this%n_atom
          do idir = 1, this%space%dim
            write(aux, '(a1,a1,i3,a1,i3,a1)') out_label(iout),'(', iat, ',', idir, ')'
            call write_iter_header(this%output_handle(iout), aux)
          end do
        end do
        call write_iter_nl(this%output_handle(iout))

        ! second line: units
        call write_iter_string(this%output_handle(iout), '#[Iter n.]')
        call write_iter_header(this%output_handle(iout), '[' // trim(units_abbrev(units_out%time)) // ']')
      end do

      call write_iter_string(this%output_handle(1), &
        'Position in '   // trim(units_abbrev(units_out%length)))
      call write_iter_string(this%output_handle(2), &
        ', Force in '    // trim(units_abbrev(units_out%force)))

      do iout = 1, 2
        call write_iter_nl(this%output_handle(iout))
        call write_iter_string(this%output_handle(iout),'#######################################################################')
        call write_iter_nl(this%output_handle(iout))
      end do
    end if

    call write_iter_start(this%output_handle(1))
    call write_iter_start(this%output_handle(2))

    do iat = 1, this%n_atom
    ! Position
      tmp(1:this%space%dim) = units_from_atomic(units_out%length, this%coords(1:this%space%dim, iat))
      call write_iter_double(this%output_handle(1), tmp, this%space%dim) 
    ! Force
      tmp(1:this%space%dim) = units_from_atomic(units_out%force, this%tot_force(1:this%space%dim, iat))
      call write_iter_double(this%output_handle(2), tmp, this%space%dim)
    end do

    call write_iter_nl(this%output_handle(1))
    call write_iter_nl(this%output_handle(2))

    POP_SUB(system_dftb_output_write)
  end subroutine system_dftb_output_write

  ! ---------------------------------------------------------
  subroutine system_dftb_update_quantity(this, iq)
    class(system_dftb_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(system_dftb_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_dftb_update_quantity)
  end subroutine system_dftb_update_quantity

  ! ---------------------------------------------------------
  subroutine system_dftb_update_exposed_quantity(partner, iq)
    class(system_dftb_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(system_dftb_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_dftb_update_exposed_quantity)
  end subroutine system_dftb_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine system_dftb_copy_quantities_to_interaction(partner, interaction)
    class(system_dftb_t),          intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(system_dftb_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(system_dftb_copy_quantities_to_interaction)
  end subroutine system_dftb_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine system_dftb_update_interactions_start(this)
    class(system_dftb_t), intent(inout) :: this

    PUSH_SUB(system_dftb_update_interactions_start)

    ! Store previous force, as it is used as SCF criterium
    this%prev_tot_force(1:this%space%dim, 1:this%n_atom) = this%tot_force(1:this%space%dim, 1:this%n_atom)

    POP_SUB(system_dftb_update_interactions_start)
  end subroutine system_dftb_update_interactions_start

  ! ---------------------------------------------------------
  subroutine system_dftb_update_interactions_finish(this)
    class(system_dftb_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(system_dftb_update_interactions_finish)

    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class default
        message(1) = "Interactions not implemented for DFTB+ systems."
        call messages_fatal(1)
      end select
    end do

    POP_SUB(system_dftb_update_interactions_finish)
  end subroutine system_dftb_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine system_dftb_finalize(this)
    type(system_dftb_t), intent(inout) :: this

    PUSH_SUB(system_dftb_finalize)

    SAFE_DEALLOCATE_A(this%coords)
    SAFE_DEALLOCATE_A(this%acc)
    SAFE_DEALLOCATE_A(this%vel)
    SAFE_DEALLOCATE_A(this%tot_force)
    SAFE_DEALLOCATE_A(this%prev_tot_force)
    SAFE_DEALLOCATE_A(this%gradients)
    SAFE_DEALLOCATE_A(this%species)
    SAFE_DEALLOCATE_A(this%mass)
    SAFE_DEALLOCATE_P(this%geo)
    call ion_dynamics_end(this%ions)

    ! No call to safe_deallocate macro here, as it gives an ICE with gfortran
    if (associated(this%ext_lasers)) then
      deallocate(this%ext_lasers)
    end if

    call system_end(this)

    POP_SUB(system_dftb_finalize)
  end subroutine system_dftb_finalize

end module system_dftb_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
