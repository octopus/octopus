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
  use force_interaction_oct_m
  use geometry_oct_m
  use global_oct_m
  use interaction_oct_m
  use lorentz_force_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_beeman_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::           &
    system_dftb_t,    &
    system_dftb_init

   type, extends(system_t) :: system_dftb_t
    FLOAT :: mass
    FLOAT :: pos(1:MAX_DIM)
    FLOAT :: vel(1:MAX_DIM)
    FLOAT :: acc(1:MAX_DIM)
    FLOAT, allocatable :: prev_acc(:,:) !< A storage of the prior times.
    FLOAT :: save_pos(1:MAX_DIM)   !< A storage for the SCF loops
    FLOAT :: save_vel(1:MAX_DIM)   !< A storage for the SCF loops
    FLOAT :: tot_force(1:MAX_DIM)
    FLOAT :: prev_tot_force(1:MAX_DIM) !< Used for the SCF convergence criterium
    FLOAT, allocatable :: prev_pos(:, :) !< Used for extrapolation
    FLOAT, allocatable :: prev_vel(:, :) !< Used for extrapolation
    FLOAT :: hamiltonian_elements(1:MAX_DIM)
    FLOAT :: scc_tolerance

    type(geometry_t) :: geo
    type(c_ptr) :: output_handle
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

    integer :: nAtom, ii
    integer, parameter :: nExtChrg = 2
    character(len=MAX_PATH_LEN) :: slako_dir

    ! H2O atom types
    integer, allocatable :: species(:)

    FLOAT, allocatable :: coords(:,:), gradients(:,:)

#ifdef HAVE_DFTBPLUS
    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input
    type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pAnalysis
    type(fnode), pointer :: pParserOpts
#endif

    PUSH_SUB(system_dftb_init)

    this%namespace = namespace

    call messages_print_stress(stdout, "DFTB+ System", namespace=namespace)

    call space_init(this%space, namespace)
    this%geo%space => this%space
    call geometry_init_xyz(this%geo, namespace)
    nAtom = this%geo%natoms
    SAFE_ALLOCATE(coords(3, nAtom))
    SAFE_ALLOCATE(gradients(3, nAtom))
    SAFE_ALLOCATE(species(nAtom))
    do ii = 1, nAtom
      coords(1:3,ii) = this%geo%atom(ii)%x(1:3)
    end do
    species = [1, 2, 2]

    !%Variable SlakoDir
    !%Type string
    !%Default "./"
    !%Section Execution::IO
    !%Description
    !% Folder containing the Slako files
    !%End
    call parse_variable(namespace, 'SlakoDir', './', slako_dir)


#ifdef HAVE_DFTBPLUS
    call TDftbPlus_init(dftbp, mpicomm=mpi_world%comm)

    call dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call setChild(pRoot, "Geometry", pGeo)
    call setChildValue(pGeo, "Periodic", .false.)
    call setChildValue(pGeo, "TypeNames", ["O", "H"])
    call setChildValue(pGeo, "TypesAndCoordinates", reshape(species, [1, size(species)]), coords)
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
    call setChildValue(pMaxAng, "O", "p")
    call setChildValue(pMaxAng, "H", "s")

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

    message(1) = 'Input tree in HSD format:'
    call messages_info(1)
    call dumpHsd(input%hsdTree, stdout)

    ! initialise the DFTB+ calculator
    call dftbp%setupCalculator(input)

    !call dftbp%setGeometry(coords)
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


    POP_SUB(system_dftb_initial_conditions)
  end subroutine system_dftb_initial_conditions

  ! ---------------------------------------------------------
  subroutine system_dftb_do_td(this, operation)
    class(system_dftb_t),    intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    !integer :: ii, sdim
    !LOAT, allocatable :: tmp_pos(:, :), tmp_vel(:, :)
    !LOAT :: factor

    PUSH_SUB(system_dftb_do_td)

    !sdim = this%space%dim

    select case (operation%id)
    case (SKIP)
      ! Do nothing
    case (STORE_CURRENT_STATUS)
      this%save_pos(1:this%space%dim) = this%pos(1:this%space%dim)
      this%save_vel(1:this%space%dim) = this%vel(1:this%space%dim)

    case (VERLET_START)
      SAFE_ALLOCATE(this%prev_acc(1:this%space%dim, 1))
      !this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim) / this%mass

    case (VERLET_FINISH)
      SAFE_DEALLOCATE_A(this%prev_acc)

    case (VERLET_UPDATE_POS)
      !this%pos(1:this%space%dim) = this%pos(1:this%space%dim) + this%prop%dt * this%vel(1:this%space%dim) &
      !                           + M_HALF * this%prop%dt**2 * this%acc(1:this%space%dim)

      !this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK

    case (VERLET_COMPUTE_ACC)
      !do ii = size(this%prev_acc, dim=2) - 1, 1, -1
        !this%prev_acc(1:this%space%dim, ii + 1) = this%prev_acc(1:this%space%dim, ii)
      !end do
      !this%prev_acc(1:this%space%dim, 1) = this%acc(1:this%space%dim)
      !this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim) / this%mass

    case (VERLET_COMPUTE_VEL)
      !this%vel(1:this%space%dim) = this%vel(1:this%space%dim) &
      !  + M_HALF * this%prop%dt * (this%prev_acc(1:this%space%dim, 1) + this%acc(1:this%space%dim))

      !this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

   POP_SUB(system_dftb_do_td)
  end subroutine system_dftb_do_td

  ! ---------------------------------------------------------
  logical function system_dftb_is_tolerance_reached(this, tol) result(converged)
    class(system_dftb_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(system_dftb_is_tolerance_reached)

    ! Here we put the criterion that acceleration change is below the tolerance
    converged = .false.
    if ( (sum((this%prev_tot_force(1:this%space%dim) - this%tot_force(1:this%space%dim))**2)/ this%mass) < tol**2) then
      converged = .true.
    end if

    if (debug%info) then
      write(message(1), '(a, e12.6, a, e12.6)') "Debug: -- Change in acceleration  ", &
        sqrt(sum((this%prev_tot_force(1:this%space%dim) - this%tot_force(1:this%space%dim))**2))/this%mass, &
        " and tolerance ", tol
      call messages_info(1)
    end if

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
    write(message(2),fmt) "Coordinates: ", (this%pos(idir), idir = 1, this%space%dim)
    write(message(3),fmt) "Velocity:    ", (this%vel(idir), idir = 1, this%space%dim)
    write(message(4),fmt) "Force:       ", (this%tot_force(idir), idir = 1, this%space%dim)
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
      call write_iter_init(this%output_handle, 0, this%prop%dt, trim(io_workpath("td.general/coordinates", this%namespace)))
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
      call write_iter_end(this%output_handle)
    end if

    POP_SUB(system_dftb_output_finish)
  end subroutine system_dftb_output_finish

  ! ---------------------------------------------------------
  subroutine system_dftb_output_write(this)
    class(system_dftb_t), intent(inout) :: this

    integer :: idir
    character(len=50) :: aux
    FLOAT :: tmp(MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(system_dftb_output_write)

    if (this%clock%get_tick() == 0) then
      ! header
      call write_iter_clear(this%output_handle)
      call write_iter_string(this%output_handle,'################################################################################')
      call write_iter_nl(this%output_handle)
      call write_iter_string(this%output_handle,'# HEADER')
      call write_iter_nl(this%output_handle)

      ! first line: column names
      call write_iter_header_start(this%output_handle)

      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'x(', idir, ')'
        call write_iter_header(this%output_handle, aux)
      end do
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'v(', idir, ')'
        call write_iter_header(this%output_handle, aux)
      end do
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'f(', idir, ')'
        call write_iter_header(this%output_handle, aux)
      end do
      call write_iter_nl(this%output_handle)

      ! second line: units
      call write_iter_string(this%output_handle, '#[Iter n.]')
      call write_iter_header(this%output_handle, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_string(this%output_handle, &
        'Position in '   // trim(units_abbrev(units_out%length))   //   &
        ', Velocity in '// trim(units_abbrev(units_out%velocity)) //   &
        ', Force in '    // trim(units_abbrev(units_out%force)))
      call write_iter_nl(this%output_handle)

      call write_iter_string(this%output_handle,'################################################################################')
      call write_iter_nl(this%output_handle)
    end if

    call write_iter_start(this%output_handle)

    ! Position
    !tmp(1:this%space%dim) = units_from_atomic(units_out%length, this%pos(1:this%space%dim))
    !call write_iter_double(this%output_handle, tmp, this%space%dim)
    ! Velocity
    !tmp(1:this%space%dim) = units_from_atomic(units_out%velocity, this%vel(1:this%space%dim))
    !call write_iter_double(this%output_handle, tmp, this%space%dim)
    ! Force
    !tmp(1:this%space%dim) = units_from_atomic(units_out%force, this%tot_force(1:this%space%dim))
    !call write_iter_double(this%output_handle, tmp, this%space%dim)

    call write_iter_nl(this%output_handle)

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
    case (MASS)
      ! The classical particle has a mass, but it is not necessary to update it, as it does not change with time.
      ! We still need to set its clock, so we set it to be in sync with the particle position.
      call this%quantities(iq)%clock%set_time(this%quantities(POSITION)%clock)
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
    case (MASS)
      ! The classical particle has a mass, but it does not require any update, as it does not change with time.
      ! We still need to set its clock, so we set it to be in sync with the particle position.
      call partner%quantities(iq)%clock%set_time(partner%quantities(POSITION)%clock)
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
    ! this%prev_tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim)

    POP_SUB(system_dftb_update_interactions_start)
  end subroutine system_dftb_update_interactions_start

  ! ---------------------------------------------------------
  subroutine system_dftb_update_interactions_finish(this)
    class(system_dftb_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(system_dftb_update_interactions_finish)

    ! Compute the total force acting on the classical particle
    ! this%tot_force(1:this%space%dim) = M_ZERO
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (force_interaction_t)
        ! this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + interaction%force(1:this%space%dim)
      end select
    end do

    POP_SUB(system_dftb_update_interactions_finish)
  end subroutine system_dftb_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine system_dftb_finalize(this)
    type(system_dftb_t), intent(inout) :: this

    PUSH_SUB(system_dftb_finalize)

    call system_end(this)

    POP_SUB(system_dftb_finalize)
  end subroutine system_dftb_finalize

end module system_dftb_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
