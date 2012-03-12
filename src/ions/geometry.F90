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

module geometry_m
  use atom_m
  use c_pointer_m
  use datasets_m
  use distributed_m
  use global_m
  use io_m
  use json_m
  use loct_m
  use loct_math_m
  use messages_m
  use multicomm_m
  use mpi_m
  use parser_m
  use profiling_m
  use space_m
  use species_m
  use string_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use xyz_file_m

  implicit none

  private
  public ::                          &
    geometry_t,                      &
    geometry_init,                   &
    geometry_init_xyz,               &
    geometry_init_species,           &
    geometry_init_from_data_object,  &
    geometry_partition,              &
    geometry_create_data_object,     &
    geometry_copy,                   &
    geometry_end,                    &
    geometry_dipole,                 &
    geometry_min_distance,           &
    geometry_mass,                   &
    cm_pos,                          &
    cm_vel,                          &
    geometry_write_xyz,              &
    loadPDB,                         &
    geometry_val_charge,             &
    geometry_grid_defaults,          &
    geometry_species_time_dependent

  integer, parameter, public :: &
    INTERACTION_COULOMB = 1,    &
    INTERACTION_LJ      = 2

  integer, parameter, public :: &
    LJ_EPSILON = 1,             &
    LJ_SIGMA   = 2

  type geometry_t
    type(space_t), pointer :: space
    integer                :: natoms
    type(atom_t), pointer  :: atom(:)

    integer :: ncatoms              !< For QM+MM calculations
    type(atom_classical_t), pointer :: catom(:)

    integer :: nspecies
    type(species_t), pointer :: species(:)

    logical :: only_user_def        !< Do we want to treat only user-defined species?
    logical :: species_time_dependent !< For time-dependent user defined species

    FLOAT :: kinetic_energy         !< the ion kinetic energy

    logical :: nlpp                 !< does any species have non-local pp?
    logical :: nlcc                 !< does any species have non-local core corrections?

    type(distributed_t) :: atoms_dist

    integer, pointer :: ionic_interaction_type(:, :)
    FLOAT,   pointer :: ionic_interaction_parameter(:, :, :)

    logical          :: reduced_coordinates !< If true the coordinates are stored in reduced coordinates and need to be converted.
  end type geometry_t

contains

  ! ---------------------------------------------------------
  subroutine geometry_init(geo, space, print_info)
    type(geometry_t),           intent(inout) :: geo
    type(space_t),    target,   intent(in)    :: space
    logical,          optional, intent(in)    :: print_info

    PUSH_SUB(geometry_init)

    geo%space => space

    ! initialize geometry
    call geometry_init_xyz(geo)
    call geometry_init_species(geo, print_info=print_info)
    call distributed_nullify(geo%atoms_dist, geo%natoms)

    POP_SUB(geometry_init)
  end subroutine geometry_init


  ! ---------------------------------------------------------------
  !> initializes the xyz positions of the atoms in the structure geo
  subroutine geometry_init_xyz(geo)
    type(geometry_t), intent(inout) :: geo
    !
    type(xyz_file_info) :: xyz
    integer             :: ia
    logical             :: move
    !
    PUSH_SUB(geometry_init_xyz)
    ! load positions of the atoms
    call xyz_file_init(xyz)

    !%Variable PDBCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If this variable is present, the program tries to read the atomic coordinates
    !% from the file specified by its value. The PDB (Protein Data Bank,
    !% <tt>http://www.rcsb.org/pdb/</tt>) format is quite complicated, and it goes 
    !% well beyond the scope of this manual. You can find a comprehensive
    !% description <a href='http://www.wwpdb.org/docs.html'>here</a>.
    !% From the plethora of instructions defined in the PDB standard, <tt>Octopus</tt>
    !% only reads two, <tt>ATOM</tt> and <tt>HETATOM</tt>. From these fields, it reads:
    !% <ul>
    !% <li> columns 13-16: The species; in fact <tt>Octopus</tt> only cares about the
    !% first letter - "CA" and "CB" will both refer to carbon - so elements whose
    !% chemical symbol has more than one letter cannot be represented in this way.
    !% So, if you want to run mercury (Hg), please use one of the other two methods
    !% to input the coordinates: <tt>XYZCoordinates</tt> or <tt>Coordinates</tt>.</li>
    !% <li> columns 18-21: The residue. If residue is <tt>QM</tt>, the atom is treated by quantum
    !% mechanics; otherwise it is simply treated as an external classical point charge.
    !% Its charge will be given by columns 61-65.</li>
    !% <li> columns 31-54: The Cartesian coordinates. The Fortran format is <tt>(3f8.3)</tt>.</li>
    !% <li> columns 61-65: Classical charge of the atom. The Fortran format is <tt>(f6.2)</tt>.</li>
    !% </ul>
    !%End

    !%Variable XYZCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If <tt>PDBCoordinates</tt> is not present, the program reads the atomic coordinates from
    !% the XYZ file specified by the variable <tt>XYZCoordinates</tt> -- in case this variable
    !% is present. The XYZ format is very simple: The first line of the file has an integer
    !% indicating the number of atoms. The second can contain comments that are simply ignored by
    !% <tt>Octopus</tt>. Then there follows one line per atom, containing the chemical species and
    !% the Cartesian coordinates of the atom.
    !%End

    !%Variable Coordinates
    !%Type block
    !%Section System::Coordinates
    !%Description
    !% If neither <tt>XYZCoordinates</tt> nor <tt>PDBCoordinates</tt> was found, <tt>Octopus</tt>
    !% tries to read the coordinates for the atoms from the block <tt>Coordinates</tt>. The
    !% format is quite straightforward:
    !%
    !% <tt>%Coordinates
    !% <br>&nbsp;&nbsp;'C' |      -0.56415 | 0.0 | 0.0 | no
    !% <br>&nbsp;&nbsp;'O' | &nbsp;0.56415 | 0.0 | 0.0 | no
    !% <br>%</tt>
    !%
    !% The first line defines a carbon atom at coordinates (-0.56415, 0.0, 0.0),
    !% that is <b>not</b> allowed to move during dynamical simulations. The second line has
    !% a similar meaning. This block obviously defines a carbon monoxide molecule, if the
    !% input units are <tt>eV_Angstrom</tt>. The number of coordinates for each species
    !% must be equal to the dimension of your space (generally 3).
    !% Note that in this way it is possible to fix some of the atoms (this
    !% is not possible when specifying the coordinates through a <tt>PDBCoordinates</tt> or
    !% <tt>XYZCoordinates</tt> file). It is always possible to fix <b>all</b> atoms using the <tt>MoveIons</tt> directive.
    !%End

    !%Variable ReducedCoordinates
    !%Type block
    !%Section System::Coordinates
    !%Description
    !% This block gives the atomic coordinates relative to the real
    !% space unit cell. The format is the same as the
    !% <tt>Coordinates</tt> block.
    !%
    !% Note that in Octopus the origin of coordinates is in the center
    !% of the cell, so the coordinates inside the cell are in the
    !% range [-0.5, 0.5).
    !%
    !% This block cannot be used with the <tt>minimum</tt> box shapes.
    !% 
    !%End

    call xyz_file_read('Coordinates', xyz, geo%space)

    ! copy information from xyz to geo
    geo%natoms = xyz%n
    nullify(geo%atom)
    if(geo%natoms>0)then
      SAFE_ALLOCATE(geo%atom(1:geo%natoms))
      do ia = 1, geo%natoms
        move=.true.
        if(iand(xyz%flags, XYZ_FLAGS_MOVE).ne.0)move=xyz%atom(ia)%move
        call atom_init(geo%atom(ia), xyz%atom(ia)%label, xyz%atom(ia)%x, move=move)
      end do
    end if

    geo%reduced_coordinates = xyz%file_type == XYZ_FILE_REDUCED

    call xyz_file_end(xyz)

    ! load positions of the classical atoms, if any
    call xyz_file_init(xyz)
    nullify(geo%catom)
    geo%ncatoms = 0
    call xyz_file_read('Classical', xyz, geo%space)
    if(xyz%file_type .ne. XYZ_FILE_ERR) then ! found classical atoms
      if(.not. iand(xyz%flags, XYZ_FLAGS_CHARGE) .ne. 0) then
        message(1) = "Need to know charge for the classical atoms."
        message(2) = "Please use a .pdb"
        call messages_fatal(2)
      end if
      geo%ncatoms = xyz%n
      write(message(1), '(a,i8)') 'Info: Number of classical atoms = ', geo%ncatoms
      call messages_info(1)
      if(geo%ncatoms>0)then
        SAFE_ALLOCATE(geo%catom(1:geo%ncatoms))
        do ia = 1, geo%ncatoms
          call atom_classical_init(geo%catom(ia), xyz%atom(ia)%label, xyz%atom(ia)%x, xyz%atom(ia)%charge)
        end do
      end if
      call xyz_file_end(xyz)
    end if

    if(geometry_atoms_are_too_close(geo)) then
      ! complain
      write(message(1), '(a)') "Some of the atoms seem to sit too close to each other."
      write(message(2), '(a)') "Please review your input files and the output geometry."
      ! then write out the geometry, whether asked for or not in Output variable
      call geometry_write_xyz(STATIC_DIR, "geometry", geo, geo%space%dim)
      call messages_fatal(2)
    end if

    POP_SUB(geometry_init_xyz)
  end subroutine geometry_init_xyz


  ! ---------------------------------------------------------
  subroutine geometry_init_species(geo, print_info)
    type(geometry_t),  intent(inout) :: geo
    logical, optional, intent(in)    :: print_info

    logical :: print_info_, spec_user_defined
    integer :: i, j, k, ispin

    PUSH_SUB(geometry_init_species)

    print_info_ = .true.
    if(present(print_info)) then
      print_info_ = print_info
    end if
    ! First, count the species
    geo%nspecies = 0
    atoms1:  do i = 1, geo%natoms
      do j = 1, i - 1
        if(atom_same_species(geo%atom(j), geo%atom(i))) cycle atoms1
      end do
      geo%nspecies = geo%nspecies + 1
    end do atoms1

    ! Allocate the species structure.
    SAFE_ALLOCATE(geo%species(1:geo%nspecies))

    ! Now, read the data.
    k = 0
    geo%only_user_def = .true.
    atoms2: do i = 1, geo%natoms
      do j = 1, i - 1
        if(atom_same_species(geo%atom(j), geo%atom(i))) cycle atoms2
      end do
      k = k + 1
      call species_set_label(geo%species(k), atom_get_label(geo%atom(j)))
      call species_set_index(geo%species(k), k)
      call species_read(geo%species(k))
      geo%only_user_def = (geo%only_user_def .and. (species_type(geo%species(k)) == SPEC_USDEF))
    end do atoms2

    ! Reads the spin components. This is read here, as well as in states_init,
    ! to be able to pass it to the pseudopotential initializations subroutine.
    call parse_integer(datasets_check('SpinComponents'), 1, ispin)
    if(.not.varinfo_valid_option('SpinComponents', ispin)) call input_error('SpinComponents')
    ispin = min(2, ispin)

    if(print_info_) then
      call messages_print_stress(stdout, "Species")
    end if
    do i = 1, geo%nspecies
      call species_init(geo%species(i), ispin, geo%space, print_info=print_info_)
    end do
    if(print_info_) then
      call messages_print_stress(stdout)
    end if

    !%Variable SpeciesTimeDependent
    !%Type logical
    !%Default no
    !%Section System::Species
    !%Description
    !% When this variable is set, the potential defined in the block <tt>Species</tt> is calculated
    !% and applied to the Hamiltonian at each time step. You must have at least one <tt>spec_user_defined</tt>
    !% type of species to use this.
    !%End
    call parse_logical(datasets_check('SpeciesTimeDependent'), .false., geo%species_time_dependent)
    ! we must have at least one user defined species in order to have time dependency
    do i = 1,geo%nspecies
      if(species_type(geo%species(i)) == SPEC_USDEF) then
        spec_user_defined = .true.
      end if
    end do
    if (geo%species_time_dependent .and. .not. spec_user_defined) then
      call input_error('SpeciesTimeDependent')
    end if

    !  assign species
    do i = 1, geo%natoms
      do j = 1, geo%nspecies
        if(atom_same_species(geo%atom(i), geo%species(j)))then
          call atom_set_species(geo%atom(i), geo%species(j))
          exit
        end if
      end do
    end do

    ! find out if we need non-local core corrections
    geo%nlcc = .false.
    geo%nlpp = .false.
    do i = 1, geo%nspecies
      geo%nlcc = (geo%nlcc.or.species_has_nlcc(geo%species(i)))
      geo%nlpp = (geo%nlpp .or. species_is_ps(geo%species(i)))
    end do

    call geometry_init_interaction(geo, print_info=print_info)

    POP_SUB(geometry_init_species)
  end subroutine geometry_init_species

  ! ---------------------------------------------------------

  subroutine geometry_init_interaction(geo, print_info)
    type(geometry_t),  intent(inout) :: geo
    logical, optional, intent(in)    :: print_info

    logical :: print_info_
    integer :: nrow, irow, idx1, idx2, ispecies
    type(block_t) :: blk
    character(len=LABEL_LEN)  :: label1, label2

    PUSH_SUB(geometry_init_interaction)

    print_info_ = .true.
    if(present(print_info)) then
      print_info_ = print_info
    end if

    SAFE_ALLOCATE(geo%ionic_interaction_type(1:geo%nspecies, 1:geo%nspecies))
    nullify(geo%ionic_interaction_parameter)

    ! coulomb interaction by default
    geo%ionic_interaction_type = INTERACTION_COULOMB

    !%Variable IonicInteraction
    !%Type block
    !%Section System::Species
    !%Description
    !% This block defines the type of classical interaction between
    !% ions. Each line represents the interaction between two types of
    !% species. The first two columns contain the element symbols, the
    !% next column is the type of interaction as defined below. The
    !% next columns are the parameters for the interaction (if
    !% any). Pairs not specified interact through Coulomb`s law.
    !%
    !% Note: In most cases there is no need to specify this block,
    !% since Coulomb interaction will be used by default.
    !%
    !%Option coulomb 1
    !% Particles interact according to Coulomb`s law. The interaction
    !% strength is given by the charge of the species. There are no
    !% parameters.
    !%Option lennard_jones 2
    !% (Experimental) The Lennard-Jones 12-6 model potential. It has
    !% the form <math>V(r) = 4\epsilon((\sigma/r)^12 -
    !% (\sigma/r)^6)</math>.  The next 2 columns contain the
    !% <math>\epsilon</math> and <math>\sigma</math> (given in the
    !% corresponding input file units).
    !%End

    if(parse_block(datasets_check('IonicInteraction'), blk) == 0) then
      call messages_experimental('non-Coulombian ionic interaction')
      nrow = parse_block_n(blk)

      !for the moment we consider two parameters for lj 12 6
      SAFE_ALLOCATE(geo%ionic_interaction_parameter(1:2, 1:geo%nspecies, 1:geo%nspecies))

      do irow = 0, nrow - 1
        ! get the labels
        call parse_block_string(blk, irow, 0, label1)
        call parse_block_string(blk, irow, 1, label2)

        ! and the index that corresponds to each species
        do ispecies = 1, geo%nspecies
          if(species_label(geo%species(ispecies)) == label1) idx1 = ispecies
          if(species_label(geo%species(ispecies)) == label2) idx2 = ispecies
        end do

        ! get the type of interaction
        call parse_block_integer(blk, irow, 2, geo%ionic_interaction_type(idx1, idx2))

        ! the interaction is symmetrical
        geo%ionic_interaction_type(idx2, idx1) = geo%ionic_interaction_type(idx1, idx2)

        select case(geo%ionic_interaction_type(idx1, idx2))
        case(INTERACTION_COULOMB)
          ! nothing to do
        case(INTERACTION_LJ)
          call parse_block_float(blk, irow, 3, geo%ionic_interaction_parameter(LJ_EPSILON, idx1, idx2), unit = units_inp%energy)
          call parse_block_float(blk, irow, 4, geo%ionic_interaction_parameter(LJ_SIGMA, idx1, idx2), unit = units_inp%length)

          ! interaction is symmetric
          geo%ionic_interaction_parameter(1:2, idx2, idx1) = geo%ionic_interaction_parameter(1:2, idx1, idx2)

          if(print_info_) then
            message(1) = 'Info: Interaction between '//trim(label1)//' and '//trim(label2)// &
              ' is given by the Lennard-Jones potential.'
            call messages_info(1)
          end if

        end select

      end do
    end if

    POP_SUB(geometry_init_interaction)
  end subroutine geometry_init_interaction

  ! ---------------------------------------------------------
  subroutine geometry_init_from_data_object(this, space, json)
    type(geometry_t),      intent(out) :: this
    type(space_t), target, intent(in)  :: space
    type(json_object_t),   intent(in)  :: json
    !
    type(json_object_t), pointer :: spec, atom
    type(json_array_t),  pointer :: species, atoms
    character(len=LABEL_LEN)     :: label
    integer                      :: i, j, ierr
    !
    PUSH_SUB(geometry_init_from_data_object)
    this%space=>space
    call json_get(json, "nspecies", this%nspecies, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "nspecies" from geometry data object.'
      call messages_fatal(1)
      return
    end if
    SAFE_ALLOCATE(this%species(this%nspecies))
    call json_get(json, "species", species, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "species" array from geometry data object.'
      call messages_fatal(1)
      return
    end if
    do i=1, this%nspecies
      call json_get(species, i, spec, ierr)
      if(ierr/=JSON_OK)then
        write(unit=message(1), fmt="(a,i3,a)") &
          'Could not read the ', i, 'th "species" element from geometry data object.'
        call messages_fatal(1)
        return
      end if
      call species_init_from_data_object(this%species(i), i, spec)
    end do
    call json_get(json, "natoms", this%natoms, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "natoms" from geometry data object.'
      call messages_fatal(1)
      return
    end if
    SAFE_ALLOCATE(this%atom(this%natoms))
    call json_get(json, "atom", atoms, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "atom" array from geometry data object.'
      call messages_fatal(1)
      return
    end if
    do i=1, this%natoms
      call json_get(atoms, i, atom, ierr)
      if(ierr/=JSON_OK)then
        write(unit=message(1), fmt="(a,i3,a)") &
          'Could not read the ', i, 'th "atom" element from geometry data object.'
        call messages_fatal(1)
        return
      end if
      do j=1, this%nspecies
        call json_get(atom, "label", label, ierr)
        if(ierr/=JSON_OK)then
          write(unit=message(1), fmt="(a,i3,a)") &
            'Could not read the ', i, 'th "atom" element "label" from geometry data object.'
          call messages_fatal(1)
          return
        end if
        if(trim(label)==trim(species_label(this%species(j))))then
          call atom_init_from_data_object(this%atom(i), this%species(j), atom)
          exit
        end if
      end do
    end do
    call json_get(json, "ncatoms", this%ncatoms, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "ncatoms" from geometry data object.'
      call messages_fatal(1)
      return
    end if
    SAFE_ALLOCATE(this%catom(this%ncatoms))
    call json_get(json, "catom", atoms, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "catom" array from geometry data object.'
      call messages_fatal(1)
      return
    end if
    do i=1, this%ncatoms
      call json_get(atoms, i, atom, ierr)
      if(ierr/=JSON_OK)then
        write(unit=message(1), fmt="(a,i3,a)") &
          'Could not read the ', i, 'th "catom" from geometry data object.'
        call messages_fatal(1)
        return
      end if
      call atom_classical_init_from_data_object(this%catom(i), atom)
    end do
    this%only_user_def=.false.
    this%species_time_dependent=.false.
    this%kinetic_energy=M_ZERO
    this%nlpp=.false.
    this%nlcc=.false.
    call distributed_nullify(this%atoms_dist, this%natoms)
    nullify(this%ionic_interaction_type)
    nullify(this%ionic_interaction_parameter)
    this%reduced_coordinates=.false.
    POP_SUB(geometry_init_from_data_object)
    return
  end subroutine geometry_init_from_data_object

  ! ---------------------------------------------------------
  subroutine geometry_create_data_object(this, json)
    type(geometry_t),    intent(in)  :: this
    type(json_object_t), intent(out) :: json
    !
    type(json_object_t), pointer :: spec, atom
    type(json_array_t),  pointer :: species, atoms
    integer                      :: i
    !
    PUSH_SUB(geometry_create_data_object)
    call json_init(json)
    call json_set(json, "nspecies", this%nspecies)
    SAFE_ALLOCATE(species)
    call json_init(species)
    do i=1, this%nspecies
      SAFE_ALLOCATE(spec)
      call species_create_data_object(this%species(i), spec)
      call json_append(species, spec)
      nullify(spec)
    end do
    call json_set(json, "species", species)
    nullify(species)
    call json_set(json, "natoms", this%natoms)
    SAFE_ALLOCATE(atoms)
    call json_init(atoms)
    do i=1, this%natoms
      SAFE_ALLOCATE(atom)
      call atom_create_data_object(this%atom(i), atom)
      call json_append(atoms, atom)
      nullify(atom)
    end do
    call json_set(json, "atom", atoms)
    nullify(atoms)
    call json_set(json, "ncatoms", this%ncatoms)
    SAFE_ALLOCATE(atoms)
    call json_init(atoms)
    do i=1, this%ncatoms
      SAFE_ALLOCATE(atom)
      call atom_classical_create_data_object(this%catom(i), atom)
      call json_append(atoms, atom)
      nullify(atom)
    end do
    call json_set(json, "catom", atoms)
    nullify(atoms)
    POP_SUB(geometry_create_data_object)
    return
  end subroutine geometry_create_data_object

  ! ---------------------------------------------------------
  subroutine geometry_partition(geo, mc)
    type(geometry_t),            intent(inout) :: geo
    type(multicomm_t),           intent(in)    :: mc

    PUSH_SUB(geometry_partition)

    call distributed_init(geo%atoms_dist, geo%natoms, mc%group_comm(P_STRATEGY_STATES), "atoms")

    POP_SUB(geometry_partition)
  end subroutine geometry_partition


  ! ---------------------------------------------------------
  subroutine loadPDB(iunit, geo)
    integer,          intent(in)    :: iunit
    type(geometry_t), intent(inout) :: geo

    FLOAT, dimension(3) :: x
    FLOAT               :: charge
    character(len=80)   :: record
    character(len=6)    :: record_name
    character(len=4)    :: atm
    character(len=3)    :: res
    integer             :: na, nca, ierr

    PUSH_SUB(loadPDB)

    ! First count number of atoms
    rewind(iunit)
    geo%natoms = 0
    geo%ncatoms = 0
    do
      read(iunit, '(a80)', iostat=ierr) record
      if(ierr/=0)exit
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(17x,a3)') res
        if(trim(res) == 'QM') then
          geo%natoms = geo%natoms + 1
        else
          geo%ncatoms = geo%ncatoms + 1
        end if
      end if
    end do

    SAFE_ALLOCATE(geo%atom(1:geo%natoms))
    SAFE_ALLOCATE(geo%catom(1:geo%ncatoms))

    ! read in the data
    rewind(iunit)
    na = 1; nca = 1
    do
      read(iunit, '(a80)', iostat=ierr) record
      if(ierr/=0)exit
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(12x,a4,1x,a3)') atm, res
        call str_trim(atm)
        if(trim(res) == 'QM') then
          read(record, '(30x,3f8.3)') x
          call atom_init(geo%atom(na), atm(1:1), x)
          na = na + 1
        else
          read(record, '(30x,3f8.3,6x,f6.2)') x, charge
          call atom_classical_init(geo%catom(nca), atm, x, charge)
          nca = nca + 1
        end if
      end if
    end do

    POP_SUB(loadPDB)
  end subroutine loadPDB

  ! ---------------------------------------------------------
  subroutine geometry_end(geo)
    type(geometry_t), intent(inout) :: geo

    PUSH_SUB(geometry_end)

    call distributed_end(geo%atoms_dist)

    SAFE_DEALLOCATE_P(geo%ionic_interaction_type)
    SAFE_DEALLOCATE_P(geo%ionic_interaction_parameter)
    SAFE_DEALLOCATE_P(geo%atom)
    geo%natoms=0
    SAFE_DEALLOCATE_P(geo%catom)
    geo%ncatoms=0

    call species_end(geo%nspecies, geo%species)
    SAFE_DEALLOCATE_P(geo%species)
    geo%nspecies=0

    POP_SUB(geometry_end)
  end subroutine geometry_end

  ! ---------------------------------------------------------
  ! This function returns .true. if two atoms are too close.
  logical function geometry_atoms_are_too_close(geo) result(too_close)
    type(geometry_t), intent(in) :: geo

    PUSH_SUB(geometry_atoms_are_too_close)

    too_close = (geometry_min_distance(geo) < CNST(1.0e-5) .and. geo%natoms > 1)

    POP_SUB(geometry_atoms_are_too_close)
  end function geometry_atoms_are_too_close

  ! ---------------------------------------------------------
  logical function geometry_species_time_dependent(geo) result(time_dependent)
    type(geometry_t), intent(in) :: geo

    PUSH_SUB(geometry_species_time_dependent)

    time_dependent = geo%species_time_dependent

    POP_SUB(geometry_species_time_dependent)
  end function geometry_species_time_dependent

  ! ---------------------------------------------------------
  subroutine geometry_dipole(geo, dipole)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: dipole(:)

    integer :: ia

    PUSH_SUB(geometry_dipole)

    dipole(1:geo%space%dim) = M_ZERO
    do ia = 1, geo%natoms
      dipole(1:geo%space%dim) = dipole(1:geo%space%dim) + &
        species_zval(geo%atom(ia)%spec)*geo%atom(ia)%x(1:geo%space%dim)
    end do
    dipole = P_PROTON_CHARGE*dipole

    POP_SUB(geometry_dipole)
  end subroutine geometry_dipole


  ! ---------------------------------------------------------
  FLOAT function geometry_min_distance(geo) result(rmin)
    type(geometry_t), intent(in)  :: geo

    integer :: i, j
    FLOAT :: r

    PUSH_SUB(geometry_min_distance)

    rmin = huge(rmin)
    do i = 1, geo%natoms
      do j = i + 1, geo%natoms
        r = atom_distance(geo%atom(i), geo%atom(j))
        if(r < rmin) then
          rmin = r
        end if
      end do
    end do

    POP_SUB(geometry_min_distance)
  end function geometry_min_distance

  ! ---------------------------------------------------------
  subroutine cm_pos(geo, pos)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: pos(:)

    FLOAT :: mass
    integer :: ia

    PUSH_SUB(cm_pos)

    pos = M_ZERO
    mass = M_ZERO
    do ia = 1, geo%natoms
      pos = pos + species_weight(geo%atom(ia)%spec) * geo%atom(ia)%x
      mass = mass + species_weight(geo%atom(ia)%spec)
    end do
    pos = pos/mass

    POP_SUB(cm_pos)
  end subroutine cm_pos


  ! ---------------------------------------------------------
  subroutine cm_vel(geo, vel)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: vel(:)

    FLOAT :: mass
    integer :: iatom

    PUSH_SUB(cm_vel)

    vel = M_ZERO
    mass = M_ZERO
    do iatom = 1, geo%natoms
      vel = vel + species_weight(geo%atom(iatom)%spec) * geo%atom(iatom)%v
      mass = mass + species_weight(geo%atom(iatom)%spec)
    end do
    vel = vel / mass

    POP_SUB(cm_vel)
  end subroutine cm_vel


  ! ---------------------------------------------------------
  subroutine geometry_write_xyz(dir, fname, geo, sbdim, append, comment)
    character(len=*),    intent(in) :: dir, fname
    type(geometry_t),    intent(in) :: geo
    integer,             intent(in) :: sbdim
    logical,             intent(in), optional :: append
    character(len=*),    intent(in), optional :: comment

    integer iatom, idir, iunit
    character(len=6) position

    if( .not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(atom_write_xyz)

    call io_mkdir(dir)
    position = 'asis'
    if(present(append)) then
      if(append) position = 'append'
    end if
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.xyz', action='write', position=position)

    write(iunit, '(i4)') geo%natoms
    if (present(comment)) then
      write(iunit, '(1x,a)') comment
    else
      write(iunit, '(1x,a,a)') 'units: ', trim(units_abbrev(units_out%length))
    endif
    do iatom = 1, geo%natoms
      call atom_write_xyz(geo%atom(iatom), sbdim, iunit)
    end do
    call io_close(iunit)

    if(geo%ncatoms > 0) then
      iunit = io_open(trim(dir)//'/'//trim(fname)//'_classical.xyz', action='write', position=position)
      write(iunit, '(i4)') geo%ncatoms
      write(iunit, '(1x)')
      do iatom = 1, geo%ncatoms
        call atom_classical_write_xyz(geo%catom(iatom), sbdim, iunit)
      end do
      call io_close(iunit)
    end if

    POP_SUB(atom_write_xyz)
  end subroutine geometry_write_xyz


  ! ---------------------------------------------------------
  subroutine geometry_val_charge(geo, val_charge)
    type(geometry_t), intent(in) :: geo
    FLOAT,           intent(out) :: val_charge

    integer :: iatom

    PUSH_SUB(geometry_val_charge)

    val_charge = M_ZERO
    do iatom = 1, geo%natoms
      val_charge = val_charge - species_zval(geo%atom(iatom)%spec)
    end do

    POP_SUB(geometry_val_charge)
  end subroutine geometry_val_charge


  ! ---------------------------------------------------------

  FLOAT pure function geometry_mass(geo) result(mass)
    type(geometry_t), intent(in) :: geo

    integer :: iatom

    mass = M_ZERO
    do iatom = 1, geo%natoms
      mass = mass + species_weight(geo%atom(iatom)%spec)
    end do

  end function geometry_mass

  ! ---------------------------------------------------------
  subroutine geometry_grid_defaults(geo, def_h, def_rsize)
    type(geometry_t), intent(in) :: geo
    FLOAT,           intent(out) :: def_h, def_rsize

    integer :: ispec

    PUSH_SUB(geometry_grid_defaults)

    def_h     =  huge(def_h)
    def_rsize = -huge(def_rsize)
    do ispec = 1, geo%nspecies
      def_h     = min(def_h,     species_def_h(geo%species(ispec)))
      def_rsize = max(def_rsize, species_def_rsize(geo%species(ispec)))
    end do

    POP_SUB(geometry_grid_defaults)
  end subroutine geometry_grid_defaults

  !--------------------------------------------------------------
  subroutine geometry_copy(geo_out, geo_in)
    type(geometry_t), intent(out) :: geo_out
    type(geometry_t), intent(in)  :: geo_in

    PUSH_SUB(geometry_copy)

    geo_out%natoms = geo_in%natoms
    SAFE_ALLOCATE(geo_out%atom(1:geo_out%natoms))
    geo_out%atom = geo_in%atom

    geo_out%ncatoms = geo_in%ncatoms
    SAFE_ALLOCATE(geo_out%catom(1:geo_out%ncatoms))
    if(geo_in%ncatoms.gt.0) then
      geo_out%catom(1:geo_out%ncatoms) = geo_in%catom(1:geo_in%ncatoms)
    end if

    geo_out%nspecies = geo_in%nspecies
    SAFE_ALLOCATE(geo_out%species(1:geo_out%nspecies))
    geo_out%species = geo_in%species

    geo_out%only_user_def     = geo_in%only_user_def
    geo_out%kinetic_energy    = geo_in%kinetic_energy
    geo_out%nlpp              = geo_in%nlpp
    geo_out%nlcc              = geo_in%nlcc

    call loct_pointer_copy(geo_out%ionic_interaction_type, geo_in%ionic_interaction_type)
    call loct_pointer_copy(geo_out%ionic_interaction_parameter, geo_in%ionic_interaction_parameter)

    call distributed_copy(geo_in%atoms_dist, geo_out%atoms_dist)

    POP_SUB(geometry_copy)
  end subroutine geometry_copy


end module geometry_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
