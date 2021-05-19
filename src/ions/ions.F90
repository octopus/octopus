!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 M. Oliveira
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

module ions_oct_m
  use atom_oct_m
  use classical_particles_oct_m
  use iso_c_binding
  use distributed_oct_m
  use global_oct_m
  use interaction_oct_m
  use io_oct_m
  use ion_interaction_oct_m
  use lalg_adv_oct_m
  use lattice_vectors_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use read_coords_oct_m
  use space_oct_m
  use species_oct_m
  use tdfunction_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m

  implicit none

  private
  public :: ions_t

  type, extends(classical_particles_t) :: ions_t
    ! Components are public by default

    type(lattice_vectors_t) :: latt

    integer                   :: natoms
    type(atom_t), allocatable :: atom(:)

    integer                             :: ncatoms  !< For QM+MM calculations
    type(atom_classical_t), allocatable :: catom(:)

    FLOAT   :: kinetic_energy      !< the ion kinetic energy
    type(distributed_t) :: atoms_dist

    !> Information about the species
    integer                      :: nspecies
    type(species_t), allocatable :: species(:)
    logical                      :: only_user_def          !< Do we want to treat only user-defined species?
    logical,         private     :: species_time_dependent !< For time-dependent user defined species

    logical                 :: force_total_enforce
    type(ion_interaction_t) :: ion_interaction
  
    !> variables for external forces over the ions
    logical,     private :: apply_global_force
    type(tdf_t), private :: global_force_function
  contains
    procedure :: copy => ions_copy
    generic   :: assignment(=) => copy
    procedure :: partition => ions_partition
    procedure :: init_interaction => ions_init_interaction
    procedure :: initial_conditions => ions_initial_conditions
    procedure :: iteration_info => ions_iteration_info
    procedure :: update_quantity => ions_update_quantity
    procedure :: update_exposed_quantity => ions_update_exposed_quantity
    procedure :: init_interaction_as_partner => ions_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => ions_copy_quantities_to_interaction
    procedure :: write_xyz => ions_write_xyz
    procedure :: read_xyz => ions_read_xyz
    procedure :: fold_atoms_into_cell => ions_fold_atoms_into_cell
    procedure :: min_distance => ions_min_distance
    procedure :: has_time_dependent_species => ions_has_time_dependent_species
    procedure :: val_charge => ions_val_charge
    procedure :: dipole => ions_dipole
    procedure :: center_of_mass => ions_center_of_mass
    procedure :: center_of_mass_vel => ions_center_of_mass_vel
    procedure :: center => ions_center
    procedure :: axis_large => ions_axis_large
    procedure :: axis_inertia => ions_axis_inertia
    procedure :: translate => ions_translate
    procedure :: rotate => ions_rotate
    procedure :: grid_defaults => ions_grid_defaults
    procedure :: grid_defaults_info => ions_grid_defaults_info
    procedure :: get_positions => ions_get_positions
    procedure :: set_positions => ions_set_positions
    procedure :: global_force => ions_global_force
    procedure :: write_crystal => ions_write_crystal
    final :: ions_finalize
  end type ions_t

  interface ions_t
    procedure ions_constructor
  end interface ions_t

contains

  ! ---------------------------------------------------------
  function ions_constructor(namespace, print_info) result(ions)
    type(namespace_t),          intent(in)    :: namespace
    logical,          optional, intent(in)    :: print_info
    class(ions_t), pointer :: ions

    type(read_coords_info) :: xyz
    integer :: ia, ierr
    character(len=100)  :: function_name
    FLOAT :: mindist
    FLOAT, parameter :: threshold = CNST(1e-5)

    PUSH_SUB(ions_constructor)

    SAFE_ALLOCATE(ions)

    ions%namespace = namespace

    call space_init(ions%space, namespace)

    call species_init_global(namespace)
    
    ! initialize geometry
    call read_coords_init(xyz)

    ! load positions of the atoms
    call read_coords_read('Coordinates', xyz, ions%space, namespace)

    if (xyz%n < 1) then
      message(1) = "Coordinates have not been defined."
      call messages_fatal(1, namespace=namespace)
    end if

    ! Initialize parent class
    call classical_particles_init(ions, xyz%n)

    ! copy information from xyz to ions
    ions%natoms = xyz%n
    SAFE_ALLOCATE(ions%atom(1:ions%natoms))
    do ia = 1, ions%natoms
      call atom_init(ions%atom(ia), xyz%atom(ia)%label)
      ions%pos(:,ia) = xyz%atom(ia)%x(1:ions%space%dim)
      if (bitand(xyz%flags, XYZ_FLAGS_MOVE) /= 0) then
        ions%fixed(ia) = .not. xyz%atom(ia)%move
      end if
    end do

    if (allocated(xyz%latvec)) then
      ! Build lattice vectors from the XSF input
      ions%latt = lattice_vectors_t(namespace, ions%space, xyz%latvec)
    else
      ! Build lattice vectors from input file
      ions%latt = lattice_vectors_t(namespace, ions%space)
    end if

    ! Convert coordinates to Cartesian in case we have reduced coordinates
    if (xyz%source == READ_COORDS_REDUCED) then
      do ia = 1, ions%natoms
        ions%pos(:, ia) = ions%latt%red_to_cart(ions%pos(:, ia))
      end do
    end if

    call read_coords_end(xyz)

    ! load positions of the classical atoms, if any
    call read_coords_init(xyz)
    ions%ncatoms = 0
    call read_coords_read('Classical', xyz, ions%space, namespace)
    if (xyz%source /= READ_COORDS_ERR) then ! found classical atoms
      if (.not. bitand(xyz%flags, XYZ_FLAGS_CHARGE) /= 0) then
        message(1) = "Need to know charge for the classical atoms."
        message(2) = "Please use a .pdb"
        call messages_fatal(2, namespace=namespace)
      end if
      ions%ncatoms = xyz%n
      write(message(1), '(a,i8)') 'Info: Number of classical atoms = ', ions%ncatoms
      call messages_info(1)
      if (ions%ncatoms>0) then
        SAFE_ALLOCATE(ions%catom(1:ions%ncatoms))
        do ia = 1, ions%ncatoms
          call atom_classical_init(ions%catom(ia), xyz%atom(ia)%label, xyz%atom(ia)%x, xyz%atom(ia)%charge)
        end do
      end if
      call read_coords_end(xyz)
    end if


    call ions_fold_atoms_into_cell(ions)
    call ions_init_species(ions, print_info=print_info)
    call distributed_nullify(ions%atoms_dist, ions%natoms)

    ! Set the masses. This needs to be done after initializing the species.
    do ia = 1, ions%natoms
      ions%mass(ia) = species_mass(ions%atom(ia)%species)
    end do

    ! Check that atoms are not too close
    if (ions%natoms > 1) then
      mindist = ions_min_distance(ions, real_atoms_only = .false.)
      if (mindist < threshold) then
        write(message(1), '(a)') "Some of the atoms seem to sit too close to each other."
        write(message(2), '(a)') "Please review your input files and the output geometry (in 'static/')."
        write(message(3), '(a, f12.6, 1x, a)') "Minimum distance = ", &
          units_from_atomic(units_out%length, mindist), trim(units_abbrev(units_out%length))
        call messages_warning(3, namespace=namespace)

        ! then write out the geometry, whether asked for or not in Output variable
        call io_mkdir(STATIC_DIR, namespace)
        call ions%write_xyz(trim(STATIC_DIR)//'/geometry')
      end if

      if (ions_min_distance(ions, real_atoms_only = .true.) < threshold) then
        message(1) = "It cannot be correct to run with physical atoms so close."
        call messages_fatal(1, namespace=namespace)
      end if
    end if


    call ion_interaction_init(ions%ion_interaction, namespace, ions%space, ions%natoms)

    !%Variable ForceTotalEnforce
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% (Experimental) If this variable is set to "yes", then the sum
    !% of the total forces will be enforced to be zero.
    !%End
    call parse_variable(namespace, 'ForceTotalEnforce', .false., ions%force_total_enforce)
    if(ions%force_total_enforce) call messages_experimental('ForceTotalEnforce')

    !%Variable TDGlobalForce
    !%Type string
    !%Section Time-Dependent
    !%Description
    !% If this variable is set, a global time-dependent force will be
    !% applied to the ions in the x direction during a time-dependent
    !% run. This variable defines the base name of the force, that
    !% should be defined in the <tt>TDFunctions</tt> block. This force
    !% does not affect the electrons directly.
    !%End
  
    if(parse_is_defined(namespace, 'TDGlobalForce')) then
  
      ions%apply_global_force = .true.
  
      call parse_variable(namespace, 'TDGlobalForce', 'none', function_name)
      call tdf_read(ions%global_force_function, namespace, trim(function_name), ierr)
  
      if(ierr /= 0) then
        call messages_write("You have enabled the GlobalForce option but Octopus could not find")
        call messages_write("the '"//trim(function_name)//"' function in the TDFunctions block.")
        call messages_fatal(namespace=namespace)
      end if
  
    else
  
      ions%apply_global_force = .false.
  
    end if

    POP_SUB(ions_constructor)
  end function ions_constructor

  ! ---------------------------------------------------------
  subroutine ions_init_species(ions, print_info)
    type(ions_t),      intent(inout) :: ions
    logical, optional, intent(in)    :: print_info

    logical :: print_info_, spec_user_defined
    integer :: i, j, k, ispin

    PUSH_SUB(ions_init_species)

    print_info_ = .true.
    if (present(print_info)) then
      print_info_ = print_info
    end if
    ! First, count the species
    ions%nspecies = 0
    atoms1:  do i = 1, ions%natoms
      do j = 1, i - 1
        if (atom_same_species(ions%atom(j), ions%atom(i))) cycle atoms1
      end do
      ions%nspecies = ions%nspecies + 1
    end do atoms1

    ! Allocate the species structure.
    SAFE_ALLOCATE(ions%species(1:ions%nspecies))

    ! Now, read the data.
    k = 0
    ions%only_user_def = .true.
    atoms2: do i = 1, ions%natoms
      do j = 1, i - 1
        if (atom_same_species(ions%atom(j), ions%atom(i))) cycle atoms2
      end do
      k = k + 1
      call species_init(ions%species(k), atom_get_label(ions%atom(j)), k)
      call species_read(ions%species(k), ions%namespace)
      ! these are the species which do not represent atoms
      ions%only_user_def = ions%only_user_def .and. .not. species_represents_real_atom(ions%species(k))
      
      if (species_is_ps(ions%species(k)) .and. ions%space%dim /= 3) then
        message(1) = "Pseudopotentials may only be used with Dimensions = 3."
        call messages_fatal(1, namespace=ions%namespace)
      end if

      if (species_type(ions%species(k)) == SPECIES_JELLIUM_SLAB) then
        if (ions%space%is_periodic() .and. ions%space%periodic_dim /= 2) then
          message(1) = "Periodic jelium slab can only be used if PeriodicDim = 2"
          call messages_fatal(1, namespace=ions%namespace)
        end if
      end if

    end do atoms2

    ! Reads the spin components. This is read here, as well as in states_init,
    ! to be able to pass it to the pseudopotential initializations subroutine.
    call parse_variable(ions%namespace, 'SpinComponents', 1, ispin)
    if (.not.varinfo_valid_option('SpinComponents', ispin)) call messages_input_error(ions%namespace, 'SpinComponents')
    ispin = min(2, ispin)

    if (print_info_) then
      call messages_print_stress(stdout, "Species", namespace=ions%namespace)
    end if
    do i = 1, ions%nspecies
      call species_build(ions%species(i), ions%namespace, ispin, ions%space%dim, print_info=print_info_)
    end do
    if (print_info_) then
      call messages_print_stress(stdout, namespace=ions%namespace)
    end if

    !%Variable SpeciesTimeDependent
    !%Type logical
    !%Default no
    !%Section System::Species
    !%Description
    !% When this variable is set, the potential defined in the block <tt>Species</tt> is calculated
    !% and applied to the Hamiltonian at each time step. You must have at least one <tt>species_user_defined</tt>
    !% type of species to use this.
    !%End
    call parse_variable(ions%namespace, 'SpeciesTimeDependent', .false., ions%species_time_dependent)
    ! we must have at least one user defined species in order to have time dependency
    do i = 1,ions%nspecies
      if (species_type(ions%species(i)) == SPECIES_USDEF) then
        spec_user_defined = .true.
      end if
    end do
    if (ions%species_time_dependent .and. .not. spec_user_defined) then
      call messages_input_error(ions%namespace, 'SpeciesTimeDependent')
    end if

    !  assign species
    do i = 1, ions%natoms
      do j = 1, ions%nspecies
        if (atom_same_species(ions%atom(i), ions%species(j))) then
          call atom_set_species(ions%atom(i), ions%species(j))
          exit
        end if
      end do
    end do

    POP_SUB(ions_init_species)
  end subroutine ions_init_species

  !--------------------------------------------------------------
  subroutine ions_copy(ions_out, ions_in)
    class(ions_t),     intent(out) :: ions_out
    class(ions_t),     intent(in)  :: ions_in

    PUSH_SUB(ions_copy)

    call classical_particles_copy(ions_out, ions_in)

    ions_out%latt = ions_in%latt

    ions_out%natoms = ions_in%natoms
    SAFE_ALLOCATE(ions_out%atom(1:ions_out%natoms))
    ions_out%atom = ions_in%atom

    ions_out%ncatoms = ions_in%ncatoms
    SAFE_ALLOCATE(ions_out%catom(1:ions_out%ncatoms))
    if (ions_in%ncatoms > 0) then
      ions_out%catom(1:ions_out%ncatoms) = ions_in%catom(1:ions_in%ncatoms)
    end if

    ions_out%nspecies = ions_in%nspecies
    SAFE_ALLOCATE(ions_out%species(1:ions_out%nspecies))
    ions_out%species = ions_in%species

    ions_out%only_user_def     = ions_in%only_user_def
    ions_out%kinetic_energy    = ions_in%kinetic_energy

    call distributed_copy(ions_in%atoms_dist, ions_out%atoms_dist)

    POP_SUB(ions_copy)
  end subroutine ions_copy

  ! ---------------------------------------------------------
  subroutine ions_partition(this, mc)
    class(ions_t),               intent(inout) :: this
    type(multicomm_t),           intent(in)    :: mc

    PUSH_SUB(ions_partition)

    call distributed_init(this%atoms_dist, this%natoms, mc%group_comm(P_STRATEGY_STATES), "atoms")

    call ion_interaction_init_parallelization(this%ion_interaction, this%natoms, mc)

    POP_SUB(ions_partition)
  end subroutine ions_partition

  ! ---------------------------------------------------------
  subroutine ions_init_interaction(this, interaction)
    class(ions_t),        target, intent(inout) :: this
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(ions_init_interaction)

    select type (interaction)
    class default
      call classical_particles_init_interaction(this, interaction)
    end select

    POP_SUB(ions_init_interaction)
  end subroutine ions_init_interaction

  ! ---------------------------------------------------------
  subroutine ions_initial_conditions(this, from_scratch)
    class(ions_t), intent(inout) :: this
    logical,       intent(in)    :: from_scratch

    PUSH_SUB(ions_initial_conditions)

    POP_SUB(ions_initial_conditions)
  end subroutine ions_initial_conditions

  ! ---------------------------------------------------------
  subroutine ions_iteration_info(this)
    class(ions_t), intent(in) :: this

    PUSH_SUB(ions_iteration_info)

    POP_SUB(ions_iteration_info)
  end subroutine ions_iteration_info

  ! ---------------------------------------------------------
  subroutine ions_update_quantity(this, iq)
    class(ions_t), intent(inout) :: this
    integer,       intent(in)    :: iq

    PUSH_SUB(ions_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      ! Other quantities should be handled by the parent class
      call classical_particles_update_quantity(this, iq)
    end select

    POP_SUB(ions_update_quantity)
  end subroutine ions_update_quantity

  ! ---------------------------------------------------------
  subroutine ions_update_exposed_quantity(partner, iq)
    class(ions_t), intent(inout) :: partner
    integer,       intent(in)    :: iq

    PUSH_SUB(ions_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      ! Other quantities should be handled by the parent class
      call classical_particles_update_exposed_quantity(partner, iq)
    end select

    POP_SUB(ions_update_exposed_quantity)
  end subroutine ions_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine ions_init_interaction_as_partner(partner, interaction)
    class(ions_t),                intent(in)    :: partner
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(ions_init_interaction_as_partner)

    select type (interaction)
    class default
      call classical_particles_init_interaction_as_partner(partner, interaction)
    end select

    POP_SUB(ions_init_interaction_as_partner)
  end subroutine ions_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine ions_copy_quantities_to_interaction(partner, interaction)
    class(ions_t),          intent(inout) :: partner
    class(interaction_t),   intent(inout) :: interaction

    PUSH_SUB(ions_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1, namespace=partner%namespace)
    end select

    POP_SUB(ions_copy_quantities_to_interaction)
  end subroutine ions_copy_quantities_to_interaction
  
  ! ---------------------------------------------------------
  subroutine ions_write_xyz(this, fname, append, comment)
    class(ions_t),              intent(in) :: this
    character(len=*),           intent(in) :: fname
    logical,          optional, intent(in) :: append
    character(len=*), optional, intent(in) :: comment

    integer :: iatom, idim, iunit
    character(len=6) position
    character(len=19) :: frmt

    if ( .not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(ions_write_xyz)

    position = 'asis'
    if (present(append)) then
      if (append) position = 'append'
    end if
    iunit = io_open(trim(fname)//'.xyz', this%namespace, action='write', position=position)

    write(iunit, '(i4)') this%natoms
    if (present(comment)) then
      write(iunit, '(1x,a)') comment
    else
      write(iunit, '(1x,a,a)') 'units: ', trim(units_abbrev(units_out%length_xyz_file))
    end if

    write(unit=frmt, fmt="(a5,i2.2,a4,i2.2,a6)") "(6x,a", LABEL_LEN, ",2x,", this%space%dim,"f12.6)"
    do iatom = 1, this%natoms
      write(unit=iunit, fmt=frmt) this%atom(iatom)%label, &
        (units_from_atomic(units_out%length_xyz_file, this%pos(idim, iatom)), idim=1, this%space%dim)
    end do
    call io_close(iunit)

    if (this%ncatoms > 0) then
      iunit = io_open(trim(fname)//'_classical.xyz', this%namespace, action='write', position=position)
      write(iunit, '(i4)') this%ncatoms
      write(iunit, '(1x)')
      do iatom = 1, this%ncatoms
        call atom_classical_write_xyz(this%catom(iatom), this%space%dim, iunit)
      end do
      call io_close(iunit)
    end if

    POP_SUB(ions_write_xyz)
  end subroutine ions_write_xyz

  ! ---------------------------------------------------------
  subroutine ions_read_xyz(this, fname, comment)
    class(ions_t),              intent(inout) :: this
    character(len=*),           intent(in)    :: fname
    character(len=*), optional, intent(in)    :: comment

    integer :: iatom, idir, iunit
    character(len=19) :: frmt, dum
    FLOAT :: tmp(this%space%dim)

    PUSH_SUB(ions_read_xyz)

    iunit = io_open(trim(fname)//'.xyz', this%namespace, action='read', position='rewind')

    read(iunit, '(i4)') this%natoms
    if (present(comment)) then
      read(iunit, *) 
    else
      read(iunit, *)  
    end if
    write(unit=frmt, fmt="(a5,i2.2,a4,i2.2,a6)") "(6x,a", LABEL_LEN, ",2x,", this%space%dim, "f12.6)"
    do iatom = 1, this%natoms
      read(unit=iunit, fmt=frmt) dum, (tmp(idir), idir=1, this%space%dim)

      this%pos(:, iatom) = units_to_atomic(units_out%length_xyz_file, tmp)
    end do
    call io_close(iunit)

    if (this%ncatoms > 0) then
      iunit = io_open(trim(fname)//'_classical.xyz', this%namespace, action='read', position='rewind')
      read(iunit, '(i4)') this%ncatoms
      read(iunit, *)
      do iatom = 1, this%ncatoms
        call atom_classical_read_xyz(this%catom(iatom), this%space%dim, iunit)
      end do
      call io_close(iunit)
    end if

    POP_SUB(ions_read_xyz)
  end subroutine ions_read_xyz

  ! ---------------------------------------------------------
  subroutine ions_fold_atoms_into_cell(this)
    class(ions_t),       intent(inout) :: this

    integer :: iatom

    PUSH_SUB(ions_fold_atoms_into_cell)

    do iatom = 1, this%natoms
      this%pos(:, iatom) = this%latt%fold_into_cell(this%pos(:, iatom))
    end do

    POP_SUB(ions_fold_atoms_into_cell)
  end subroutine ions_fold_atoms_into_cell

  ! ---------------------------------------------------------
  FLOAT function ions_min_distance(this, real_atoms_only) result(rmin)
    class(ions_t),      intent(in) :: this
    logical,  optional, intent(in) :: real_atoms_only

    integer :: iatom, jatom, idir
    FLOAT   :: xx(this%space%dim)
    logical :: real_atoms_only_
    type(species_t), pointer :: species

    PUSH_SUB(ions_min_distance)

    real_atoms_only_ = optional_default(real_atoms_only, .false.)

    rmin = huge(rmin)
    do iatom = 1, this%natoms
      call atom_get_species(this%atom(iatom), species)
      if(real_atoms_only_ .and. .not. species_represents_real_atom(species)) cycle
      do jatom = iatom + 1, this%natoms
        call atom_get_species(this%atom(iatom), species)
        if(real_atoms_only_ .and. .not. species_represents_real_atom(species)) cycle
        xx = abs(this%pos(:, iatom) - this%pos(:, jatom))
        if (this%space%is_periodic()) then
          xx = this%latt%cart_to_red(xx)
          do idir = 1, this%space%periodic_dim
            xx(idir) = xx(idir) - floor(xx(idir) + M_HALF)
          end do
          xx = this%latt%red_to_cart(xx)
        end if
        rmin = min(norm2(xx), rmin)
      end do
    end do

    if(.not. (this%only_user_def .and. real_atoms_only_)) then
      ! what if the nearest neighbors are periodic images?
      do idir = 1, this%space%periodic_dim
        rmin = min(rmin, norm2(this%latt%rlattice(:,idir)))
      end do
    end if

    POP_SUB(ions_min_distance)
  end function ions_min_distance

  ! ---------------------------------------------------------
  logical function ions_has_time_dependent_species(this) result(time_dependent)
    class(ions_t),     intent(in) :: this

    PUSH_SUB(ions_has_time_dependent_species)

    time_dependent = this%species_time_dependent

    POP_SUB(ions_has_time_dependent_species)
  end function ions_has_time_dependent_species

  ! ---------------------------------------------------------
  FLOAT function ions_val_charge(this, mask) result(val_charge)
    class(ions_t),              intent(in) :: this
    logical,          optional, intent(in) :: mask(:)

    integer :: iatom

    PUSH_SUB(ions_val_charge)

    val_charge = M_ZERO
    do iatom = 1, this%natoms
      if (present(mask)) then
        if (.not. mask(iatom)) cycle
      end if
      val_charge = val_charge - species_zval(this%atom(iatom)%species)
    end do

    POP_SUB(ions_val_charge)
  end function ions_val_charge

  ! ---------------------------------------------------------
  function ions_dipole(this, mask) result(dipole)
    class(ions_t),               intent(in) :: this
    logical,           optional, intent(in) :: mask(:)
    FLOAT :: dipole(this%space%dim)

    integer :: ia

    PUSH_SUB(ions_dipole)

    dipole = M_ZERO
    do ia = 1, this%natoms
      if (present(mask)) then
        if (.not. mask(ia)) cycle
      end if
      dipole = dipole + species_zval(this%atom(ia)%species)*this%pos(:, ia)
    end do
    dipole = P_PROTON_CHARGE*dipole

    POP_SUB(ions_dipole)
  end function ions_dipole

  ! ---------------------------------------------------------
  function ions_center_of_mass(this, mask, pseudo) result(pos)
    class(ions_t),               intent(in) :: this
    logical,           optional, intent(in) :: mask(:)
    logical,           optional, intent(in) :: pseudo !< calculate center considering all species to have equal mass.
    FLOAT :: pos(this%space%dim)

    FLOAT :: mass, total_mass
    integer :: ia

    PUSH_SUB(ions_center_of_mass)

    pos = M_ZERO
    total_mass = M_ZERO
    mass = M_ONE
    do ia = 1, this%natoms
      if (present(mask)) then
        if (.not. mask(ia)) cycle
      end if
      if (.not. optional_default(pseudo, .false.)) then
        mass = this%mass(ia)
      end if
      pos = pos + mass*this%pos(:, ia)
      total_mass = total_mass + mass
    end do
    pos = pos/total_mass

    POP_SUB(ions_center_of_mass)
  end function ions_center_of_mass

  ! ---------------------------------------------------------
  function ions_center_of_mass_vel(this) result(vel)
    class(ions_t),     intent(in) :: this
    FLOAT :: vel(this%space%dim)

    FLOAT :: mass, total_mass
    integer :: iatom

    PUSH_SUB(ions_center_of_mass_vel)

    vel = M_ZERO
    total_mass = M_ZERO
    do iatom = 1, this%natoms
      mass = this%mass(iatom)
      total_mass = total_mass + mass
      vel = vel + mass*this%vel(:, iatom)
    end do
    vel = vel/total_mass

    POP_SUB(ions_center_of_mass_vel)
  end function ions_center_of_mass_vel

  ! ---------------------------------------------------------
  function ions_center(this) result(pos)
    class(ions_t),     intent(in) :: this
    FLOAT :: pos(this%space%dim)

    FLOAT :: xmin(this%space%dim), xmax(this%space%dim)
    integer  :: iatom, idir

    PUSH_SUB(ions_center)

    xmin =  M_HUGE
    xmax = -M_HUGE
    do iatom = 1, this%natoms
      do idir = 1, this%space%dim
        if (this%pos(idir, iatom) > xmax(idir)) xmax(idir) = this%pos(idir, iatom)
        if (this%pos(idir, iatom) < xmin(idir)) xmin(idir) = this%pos(idir, iatom)
      end do
    end do

    pos = (xmax + xmin)/M_TWO

    POP_SUB(ions_center)
  end function ions_center

  ! ---------------------------------------------------------
  subroutine ions_axis_large(this, x, x2)
    class(ions_t),     intent(in)  :: this
    FLOAT,             intent(out) :: x(this%space%dim), x2(this%space%dim)

    integer  :: iatom, jatom
    FLOAT :: rmax, r, r2

    PUSH_SUB(ions_axis_large)

    ! first get the further apart atoms
    rmax = -M_HUGE
    do iatom = 1, this%natoms
      do jatom = 1, this%natoms/2 + 1
        r = norm2(this%pos(:, iatom) - this%pos(:, jatom))
        if (r > rmax) then
          rmax = r
          x = this%pos(:, iatom) - this%pos(:, jatom)
        end if
      end do
    end do
    x  = x /norm2(x)

    ! now let us find out what is the second most important axis
    rmax = -M_HUGE
    do iatom = 1, this%natoms
      r2 = sum(x * this%pos(:, iatom))
      r = norm2(this%pos(:, iatom) - r2*x)
      if (r > rmax) then
        rmax = r
        x2 = this%pos(:, iatom) - r2*x
      end if
    end do

    POP_SUB(ions_axis_large)
  end subroutine ions_axis_large

  ! ---------------------------------------------------------
  !> This subroutine assumes that the origin of the coordinates is the
  !! center of mass of the system
  subroutine ions_axis_inertia(this, x, x2, pseudo)
    class(ions_t),     intent(in)  :: this
    FLOAT,             intent(out) :: x(this%space%dim), x2(this%space%dim)
    logical,           intent(in)  :: pseudo !< calculate axis considering all species to have equal mass.

    FLOAT :: mass, tinertia(this%space%dim, this%space%dim), eigenvalues(this%space%dim)
    integer :: ii, jj, iatom
    type(unit_t) :: unit

    PUSH_SUB(ions_axis_inertia)

    ! first calculate the inertia tensor
    tinertia = M_ZERO
    mass = M_ONE
    do iatom = 1, this%natoms
      if (.not.pseudo) mass = this%mass(iatom)
      do ii = 1, this%space%dim
        tinertia(ii, :) = tinertia(ii, :) - mass*this%pos(ii, iatom)*this%pos(:, iatom)
        tinertia(ii, ii) = tinertia(ii, ii) + mass*sum(this%pos(:, iatom)**2)
      end do
    end do

    unit = units_out%length**2
    ! note: we always use amu for atomic masses, so no unit conversion to/from atomic is needed.
    if (pseudo) then
      write(message(1),'(a)') 'Moment of pseudo-inertia tensor [' // trim(units_abbrev(unit)) // ']'
    else
      write(message(1),'(a)') 'Moment of inertia tensor [amu*' // trim(units_abbrev(unit)) // ']'
    end if
    call messages_info(1)
    call output_tensor(stdout, tinertia, this%space%dim, unit, write_average = .true.)

    call lalg_eigensolve(this%space%dim, tinertia, eigenvalues)

    write(message(1),'(a,6f25.6)') 'Eigenvalues: ', units_from_atomic(unit, eigenvalues)
    call messages_info(1)

    ! make a choice to fix the sign of the axis.
    do ii = 1, 2
      jj = maxloc(abs(tinertia(:,ii)), dim = 1)
      if (tinertia(jj,ii) < M_ZERO) tinertia(:,ii) = -tinertia(:,ii)
    end do
    x  = tinertia(:,1)
    x2 = tinertia(:,2)

    POP_SUB(ions_axis_inertia)
  end subroutine ions_axis_inertia

  ! ---------------------------------------------------------
  subroutine ions_translate(this, xx)
    class(ions_t),     intent(inout) :: this
    FLOAT,             intent(in)    :: xx(this%space%dim)

    integer  :: iatom

    PUSH_SUB(ions_translate)

    do iatom = 1, this%natoms
      this%pos(:, iatom) = this%pos(:, iatom) - xx
    end do
    do iatom = 1, this%ncatoms
      this%catom(iatom)%x(1:this%space%dim) = this%catom(iatom)%x(1:this%space%dim) - xx
    end do

    POP_SUB(ions_translate)
  end subroutine ions_translate

  ! ---------------------------------------------------------
  subroutine ions_rotate(this, from, from2, to)
    class(ions_t),     intent(inout) :: this
    FLOAT,             intent(in)    :: from(this%space%dim)   !< assumed to be normalized
    FLOAT,             intent(in)    :: from2(this%space%dim)  !< assumed to be normalized
    FLOAT,             intent(in)    :: to(this%space%dim)     !< assumed to be normalized

    integer :: iatom, idim
    FLOAT :: m1(3, 3), m2(3, 3)
    FLOAT :: m3(3, 3), f2(3), per(3)
    FLOAT :: alpha, r

    PUSH_SUB(ions_rotate)

    if (this%space%dim /= 3) then
      call messages_not_implemented("ions_rotate in other than 3 dimensions", namespace=this%namespace)
    end if

    ! initialize matrices
    m1 = M_ZERO
    do idim = 1, this%space%dim
      m1(idim, idim) = M_ONE
    end do

    ! rotate the to-axis to the z-axis
    if (to(2) /= M_ZERO) then
      alpha = atan2(to(2), to(1))
      call rotate(m1, alpha, 3)
    end if
    alpha = atan2(sqrt(to(1)**2 + to(2)**2), to(3))
    call rotate(m1, -alpha, 2)

    ! get perpendicular to z and from
    f2 = matmul(m1, from)
    per(1) = -f2(2)
    per(2) =  f2(1)
    per(3) = M_ZERO
    r = norm2(per)
    if (r > M_ZERO) then
      per = per/r
    else
      per(2) = M_ONE
    end if

    ! rotate perpendicular axis to the y-axis
    m2 = M_ZERO; m2(1,1) = M_ONE; m2(2,2) = M_ONE; m2(3,3) = M_ONE
    alpha = atan2(per(1), per(2))
    call rotate(m2, -alpha, 3)

    ! rotate from => to (around the y-axis)
    m3 = M_ZERO; m3(1,1) = M_ONE; m3(2,2) = M_ONE; m3(3,3) = M_ONE
    alpha = acos(sum(from*to))
    call rotate(m3, -alpha, 2)

    ! join matrices
    m2 = matmul(transpose(m2), matmul(m3, m2))

    ! rotate around the z-axis to get the second axis
    per = matmul(m2, matmul(m1, from2))
    alpha = atan2(per(1), per(2))
    call rotate(m2, -alpha, 3) ! second axis is now y

    ! get combined transformation
    m1 = matmul(transpose(m1), matmul(m2, m1))

    ! now transform the coordinates
    ! it is written in this way to avoid what I consider a bug in the Intel compiler
    do iatom = 1, this%natoms
      f2 = this%pos(:, iatom)
      this%pos(:, iatom) = matmul(m1, f2)
    end do

    do iatom = 1, this%ncatoms
      f2 = this%catom(iatom)%x(1:this%space%dim)
      this%catom(iatom)%x(1:this%space%dim) = matmul(m1, f2)
    end do

    POP_SUB(ions_rotate)
  contains

    ! ---------------------------------------------------------
    subroutine rotate(m, angle, dir)
      FLOAT,   intent(inout) :: m(3, 3)
      FLOAT,   intent(in)    :: angle
      integer, intent(in)    :: dir

      FLOAT :: aux(3, 3), ca, sa

      PUSH_SUB(ions_rotate.rotate)

      ca = cos(angle)
      sa = sin(angle)

      aux = M_ZERO
      select case (dir)
      case (1)
        aux(1, 1) = M_ONE
        aux(2, 2) = ca
        aux(3, 3) = ca
        aux(2, 3) = sa
        aux(3, 2) = -sa
      case (2)
        aux(2, 2) = M_ONE
        aux(1, 1) = ca
        aux(3, 3) = ca
        aux(1, 3) = sa
        aux(3, 1) = -sa
      case (3)
        aux(3, 3) = M_ONE
        aux(1, 1) = ca
        aux(2, 2) = ca
        aux(1, 2) = sa
        aux(2, 1) = -sa
      end select

      m = matmul(aux, m)

      POP_SUB(ions_rotate.rotate)
    end subroutine rotate

  end subroutine ions_rotate

  ! ---------------------------------------------------------
  subroutine ions_grid_defaults(this, def_h, def_rsize)
    class(ions_t),     intent(in)  :: this
    FLOAT,             intent(out) :: def_h, def_rsize

    integer :: ispec

    PUSH_SUB(ions_grid_defaults)

    def_h     =  huge(def_h)
    def_rsize = -huge(def_rsize)
    do ispec = 1, this%nspecies
      def_h     = min(def_h,     species_def_h(this%species(ispec)))
      def_rsize = max(def_rsize, species_def_rsize(this%species(ispec)))
    end do

    POP_SUB(ions_grid_defaults)
  end subroutine ions_grid_defaults

  ! ---------------------------------------------------------
  subroutine ions_grid_defaults_info(this)
    class(ions_t),     intent(in) :: this

    integer :: ispec

    PUSH_SUB(ions_grid_defaults_info)

    do ispec = 1, this%nspecies
      call messages_write("Species '"//trim(species_label(this%species(ispec)))//"': spacing = ")
      if (species_def_h(this%species(ispec)) > M_ZERO) then
        call messages_write(species_def_h(this%species(ispec)), fmt = '(f7.3)')
        call messages_write(" b")
      else
        call messages_write(" unknown")
      end if
      call messages_write(", radius = ")
      if (species_def_rsize(this%species(ispec)) > M_ZERO) then
        call messages_write(species_def_rsize(this%species(ispec)), fmt = '(f5.1)')
        call messages_write(" b.")
      else
        call messages_write(" unknown.")
      end if
      call messages_info()
    end do

    POP_SUB(ions_grid_defaults_info)
  end subroutine ions_grid_defaults_info

  ! ---------------------------------------------------------
  subroutine ions_get_positions(this, q)
    class(ions_t),     intent(in)    :: this
    FLOAT,             intent(inout) :: q(:, :)

    PUSH_SUB(ions_get_positions)

    q = this%pos

    POP_SUB(ions_get_positions)
  end subroutine ions_get_positions

  ! ---------------------------------------------------------
  subroutine ions_set_positions(this, q)
    class(ions_t),     intent(inout) :: this
    FLOAT,             intent(in)    :: q(:, :)

    PUSH_SUB(ions_get_positions)

    this%pos = q

    POP_SUB(ions_get_positions)
  end subroutine ions_set_positions

  ! ---------------------------------------------------------
  function ions_global_force(this, time) result(force)
    class(ions_t),        intent(in)    :: this
    FLOAT,                intent(in)    :: time
    FLOAT :: force(this%space%dim)

    PUSH_SUB(ions_global_force)

    force = M_ZERO

    if (this%apply_global_force) then
      force(1) = units_to_atomic(units_inp%force, tdf(this%global_force_function, time))
    end if

    POP_SUB(ions_global_force)
  end function ions_global_force

  ! ----------------------------------------------------------------
  !> This subroutine creates a crystal by replicating the geometry and
  !! writes the result to dir//'crystal.xyz'
  subroutine ions_write_crystal(this, dir)
    class(ions_t),           intent(in) :: this 
    character(len=*),        intent(in) :: dir
    
    type(lattice_iterator_t) :: latt_iter
    FLOAT :: radius, pos(this%space%dim)
    integer :: iatom, icopy, iunit

    PUSH_SUB(ions_write_crystal)
    
    radius = maxval(M_HALF*norm2(this%latt%rlattice, dim=1))*(M_ONE + M_EPSILON)
    latt_iter = lattice_iterator_t(this%latt, radius)

    if(mpi_grp_is_root(mpi_world)) then
      
      iunit = io_open(trim(dir)//'/crystal.xyz', this%namespace, action='write')
      
      write(iunit, '(i9)') this%natoms*latt_iter%n_cells
      write(iunit, '(a)') '#generated by Octopus'
      
      do iatom = 1, this%natoms
        do icopy = 1, latt_iter%n_cells
          pos = units_from_atomic(units_out%length, this%pos(:, iatom) + latt_iter%get(icopy))
          write(iunit, '(a, 99f12.6)') this%atom(iatom)%label, pos
        end do
      end do

      call io_close(iunit)
    end if

    POP_SUB(ions_write_crystal)
  end subroutine ions_write_crystal

  ! ---------------------------------------------------------
  subroutine ions_finalize(ions)
    type(ions_t),     intent(inout) :: ions

    PUSH_SUB(ions_finalize)

    call classical_particles_end(ions)

    call distributed_end(ions%atoms_dist)

    call ion_interaction_end(ions%ion_interaction)

    SAFE_DEALLOCATE_A(ions%atom)
    ions%natoms=0
    SAFE_DEALLOCATE_A(ions%catom)
    ions%ncatoms=0

    call species_end(ions%nspecies, ions%species)
    SAFE_DEALLOCATE_A(ions%species)
    ions%nspecies=0

    call species_end_global()

    POP_SUB(ions_finalize)
  end subroutine ions_finalize

end module ions_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
