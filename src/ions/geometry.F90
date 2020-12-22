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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module geometry_oct_m
  use atom_oct_m
  use iso_c_binding
  use distributed_oct_m
  use global_oct_m
  use io_oct_m
  use lalg_adv_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use read_coords_oct_m
  use space_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                          &
    geometry_t,                      &
    geometry_init,                   &
    geometry_copy,                   &
    geometry_partition,              &
    geometry_write_xyz,              &
    geometry_read_xyz,               &
    geometry_min_distance,           &
    geometry_species_time_dependent, &
    geometry_val_charge,             &
    geometry_dipole,                 &
    geometry_center_of_mass,         &
    geometry_center_of_mass_vel,     &
    geometry_center,                 &
    geometry_axis_large,             &
    geometry_axis_inertia,           &
    geometry_translate,              &
    geometry_rotate,                 &
    geometry_grid_defaults,          &
    geometry_grid_defaults_info,     &
    geometry_get_positions,          &
    geometry_set_positions,          &
    geometry_end

      type geometry_t
    ! Components are public by default

    type(space_t), pointer :: space

    integer               :: natoms
    type(atom_t), pointer :: atom(:)

    integer                         :: ncatoms  !< For QM+MM calculations
    type(atom_classical_t), pointer :: catom(:)

    FLOAT   :: kinetic_energy      !< the ion kinetic energy
    logical :: reduced_coordinates !< If true the coordinates are stored in
                                   !! reduced coordinates and need to be converted.

    type(distributed_t) :: atoms_dist

    !> Information about the species
    integer                  :: nspecies
    type(species_t), pointer :: species(:)
    logical                  :: only_user_def          !< Do we want to treat only user-defined species?
    logical,         private :: species_time_dependent !< For time-dependent user defined species

    !> variables for passing info from XSF input to simul_box_init
    integer :: periodic_dim
    FLOAT :: lsize(MAX_DIM)
  end type geometry_t

contains

  ! ---------------------------------------------------------
  subroutine geometry_init(geo, namespace, space, print_info)
    type(geometry_t),           intent(inout) :: geo
    type(namespace_t),          intent(in)    :: namespace
    type(space_t),    target,   intent(in)    :: space
    logical,          optional, intent(in)    :: print_info

    PUSH_SUB(geometry_init)

    geo%space => space

    call species_init_global(namespace)
    
    ! initialize geometry
    call geometry_init_xyz(geo, namespace)
    call geometry_init_species(geo, namespace, print_info=print_info)
    call distributed_nullify(geo%atoms_dist, geo%natoms)

    POP_SUB(geometry_init)
  end subroutine geometry_init

  ! ---------------------------------------------------------------
  !> initializes the xyz positions of the atoms in the structure geo
  subroutine geometry_init_xyz(geo, namespace)
    type(geometry_t),  intent(inout) :: geo
    type(namespace_t), intent(in)    :: namespace

    type(read_coords_info) :: xyz
    integer :: ia
    logical :: move

    PUSH_SUB(geometry_init_xyz)

    call read_coords_init(xyz)

    ! load positions of the atoms
    call read_coords_read('Coordinates', xyz, geo%space, namespace)

    if (xyz%n < 1) then
      message(1) = "Coordinates have not been defined."
      call messages_fatal(1, namespace=namespace)
    end if

    ! copy information from xyz to geo
    geo%natoms = xyz%n
    nullify(geo%atom)
    if (geo%natoms>0) then
      SAFE_ALLOCATE(geo%atom(1:geo%natoms))
      do ia = 1, geo%natoms
        move=.true.
        if (bitand(xyz%flags, XYZ_FLAGS_MOVE) /= 0) move=xyz%atom(ia)%move
        call atom_init(geo%atom(ia), xyz%atom(ia)%label, xyz%atom(ia)%x, move=move)
      end do
    end if

    geo%reduced_coordinates = xyz%source == READ_COORDS_REDUCED
    geo%periodic_dim = xyz%periodic_dim
    geo%lsize(:) = xyz%lsize(:)

    call read_coords_end(xyz)

    ! load positions of the classical atoms, if any
    call read_coords_init(xyz)
    nullify(geo%catom)
    geo%ncatoms = 0
    call read_coords_read('Classical', xyz, geo%space, namespace)
    if (xyz%source /= READ_COORDS_ERR) then ! found classical atoms
      if (.not. bitand(xyz%flags, XYZ_FLAGS_CHARGE) /= 0) then
        message(1) = "Need to know charge for the classical atoms."
        message(2) = "Please use a .pdb"
        call messages_fatal(2, namespace=namespace)
      end if
      geo%ncatoms = xyz%n
      write(message(1), '(a,i8)') 'Info: Number of classical atoms = ', geo%ncatoms
      call messages_info(1)
      if (geo%ncatoms>0) then
        SAFE_ALLOCATE(geo%catom(1:geo%ncatoms))
        do ia = 1, geo%ncatoms
          call atom_classical_init(geo%catom(ia), xyz%atom(ia)%label, xyz%atom(ia)%x, xyz%atom(ia)%charge)
        end do
      end if
      call read_coords_end(xyz)
    end if

    POP_SUB(geometry_init_xyz)
  end subroutine geometry_init_xyz

  ! ---------------------------------------------------------
  subroutine geometry_init_species(geo, namespace, print_info)
    type(geometry_t),  intent(inout) :: geo
    type(namespace_t), intent(in)    :: namespace
    logical, optional, intent(in)    :: print_info

    logical :: print_info_, spec_user_defined
    integer :: i, j, k, ispin

    PUSH_SUB(geometry_init_species)

    print_info_ = .true.
    if (present(print_info)) then
      print_info_ = print_info
    end if
    ! First, count the species
    geo%nspecies = 0
    atoms1:  do i = 1, geo%natoms
      do j = 1, i - 1
        if (atom_same_species(geo%atom(j), geo%atom(i))) cycle atoms1
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
        if (atom_same_species(geo%atom(j), geo%atom(i))) cycle atoms2
      end do
      k = k + 1
      call species_init(geo%species(k), atom_get_label(geo%atom(j)), k)
      call species_read(geo%species(k), namespace)
      ! these are the species which do not represent atoms
      geo%only_user_def = geo%only_user_def .and. .not. species_represents_real_atom(geo%species(k))
      
      if (species_is_ps(geo%species(k)) .and. geo%space%dim /= 3) then
        message(1) = "Pseudopotentials may only be used with Dimensions = 3."
        call messages_fatal(1, namespace=namespace)
      end if
    end do atoms2

    ! Reads the spin components. This is read here, as well as in states_init,
    ! to be able to pass it to the pseudopotential initializations subroutine.
    call parse_variable(namespace, 'SpinComponents', 1, ispin)
    if (.not.varinfo_valid_option('SpinComponents', ispin)) call messages_input_error(namespace, 'SpinComponents')
    ispin = min(2, ispin)

    if (print_info_) then
      call messages_print_stress(stdout, "Species", namespace=namespace)
    end if
    do i = 1, geo%nspecies
      call species_build(geo%species(i), namespace, ispin, geo%space%dim, print_info=print_info_)
    end do
    if (print_info_) then
      call messages_print_stress(stdout, namespace=namespace)
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
    call parse_variable(namespace, 'SpeciesTimeDependent', .false., geo%species_time_dependent)
    ! we must have at least one user defined species in order to have time dependency
    do i = 1,geo%nspecies
      if (species_type(geo%species(i)) == SPECIES_USDEF) then
        spec_user_defined = .true.
      end if
    end do
    if (geo%species_time_dependent .and. .not. spec_user_defined) then
      call messages_input_error(namespace, 'SpeciesTimeDependent')
    end if

    !  assign species
    do i = 1, geo%natoms
      do j = 1, geo%nspecies
        if (atom_same_species(geo%atom(i), geo%species(j))) then
          call atom_set_species(geo%atom(i), geo%species(j))
          exit
        end if
      end do
    end do

    POP_SUB(geometry_init_species)
  end subroutine geometry_init_species

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
    if (geo_in%ncatoms > 0) then
      geo_out%catom(1:geo_out%ncatoms) = geo_in%catom(1:geo_in%ncatoms)
    end if

    geo_out%nspecies = geo_in%nspecies
    SAFE_ALLOCATE(geo_out%species(1:geo_out%nspecies))
    geo_out%species = geo_in%species

    geo_out%only_user_def     = geo_in%only_user_def
    geo_out%kinetic_energy    = geo_in%kinetic_energy

    call distributed_copy(geo_in%atoms_dist, geo_out%atoms_dist)

    POP_SUB(geometry_copy)
  end subroutine geometry_copy

  ! ---------------------------------------------------------
  subroutine geometry_partition(geo, mc)
    type(geometry_t),            intent(inout) :: geo
    type(multicomm_t),           intent(in)    :: mc

    PUSH_SUB(geometry_partition)

    call distributed_init(geo%atoms_dist, geo%natoms, mc%group_comm(P_STRATEGY_STATES), "atoms")

    POP_SUB(geometry_partition)
  end subroutine geometry_partition

  ! ---------------------------------------------------------
  subroutine geometry_write_xyz(geo, fname, namespace, append, comment)
    type(geometry_t),           intent(in) :: geo
    character(len=*),           intent(in) :: fname
    type(namespace_t),          intent(in) :: namespace
    logical,          optional, intent(in) :: append
    character(len=*), optional, intent(in) :: comment

    integer :: iatom, iunit
    character(len=6) position

    if ( .not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(atom_write_xyz)

    position = 'asis'
    if (present(append)) then
      if (append) position = 'append'
    end if
    iunit = io_open(trim(fname)//'.xyz', namespace, action='write', position=position)

    write(iunit, '(i4)') geo%natoms
    if (present(comment)) then
      write(iunit, '(1x,a)') comment
    else
      write(iunit, '(1x,a,a)') 'units: ', trim(units_abbrev(units_out%length_xyz_file))
    end if
    do iatom = 1, geo%natoms
      call atom_write_xyz(geo%atom(iatom), geo%space%dim, iunit)
    end do
    call io_close(iunit)

    if (geo%ncatoms > 0) then
      iunit = io_open(trim(fname)//'_classical.xyz', namespace, action='write', position=position)
      write(iunit, '(i4)') geo%ncatoms
      write(iunit, '(1x)')
      do iatom = 1, geo%ncatoms
        call atom_classical_write_xyz(geo%catom(iatom), geo%space%dim, iunit)
      end do
      call io_close(iunit)
    end if

    POP_SUB(atom_write_xyz)
  end subroutine geometry_write_xyz

  ! ---------------------------------------------------------
  subroutine geometry_read_xyz(geo, fname, namespace, comment)
    type(geometry_t),           intent(inout) :: geo
    character(len=*),           intent(in)    :: fname
    type(namespace_t),          intent(in)    :: namespace
    character(len=*), optional, intent(in)    :: comment

    integer :: iatom, iunit

    PUSH_SUB(geometry_read_xyz)

    iunit = io_open(trim(fname)//'.xyz', namespace, action='read', position='rewind')

    read(iunit, '(i4)') geo%natoms
    if (present(comment)) then
      read(iunit, *) 
    else
      read(iunit, *)  
    end if
    do iatom = 1, geo%natoms
      call atom_read_xyz(geo%atom(iatom), geo%space%dim, iunit)
    end do
    call io_close(iunit)

    if (geo%ncatoms > 0) then
      iunit = io_open(trim(fname)//'_classical.xyz', namespace, action='read', position='rewind')
      read(iunit, '(i4)') geo%ncatoms
      read(iunit, *)
      do iatom = 1, geo%ncatoms
        call atom_classical_read_xyz(geo%catom(iatom), geo%space%dim, iunit)
      end do
      call io_close(iunit)
    end if

    POP_SUB(geometry_read_xyz)
  end subroutine geometry_read_xyz

  !> Beware: this is wrong for periodic systems. Use simul_box_min_distance instead.
  ! ---------------------------------------------------------
  FLOAT function geometry_min_distance(geo, real_atoms_only) result(rmin)
    type(geometry_t),  intent(in) :: geo
    logical, optional, intent(in) :: real_atoms_only

    integer :: i, j
    FLOAT   :: r
    logical :: real_atoms_only_
    type(species_t), pointer :: species

    PUSH_SUB(geometry_min_distance)

    real_atoms_only_ = optional_default(real_atoms_only, .false.)

    rmin = huge(rmin)
    do i = 1, geo%natoms
      call atom_get_species(geo%atom(i), species)
      if (real_atoms_only_ .and. .not. species_represents_real_atom(species)) cycle
      do j = i + 1, geo%natoms
        call atom_get_species(geo%atom(i), species)
        if (real_atoms_only_ .and. .not. species_represents_real_atom(species)) cycle
        r = atom_distance(geo%atom(i), geo%atom(j))
        if (r < rmin) then
          rmin = r
        end if
      end do
    end do

    POP_SUB(geometry_min_distance)
  end function geometry_min_distance

  ! ---------------------------------------------------------
  logical function geometry_species_time_dependent(geo) result(time_dependent)
    type(geometry_t), intent(in) :: geo

    PUSH_SUB(geometry_species_time_dependent)

    time_dependent = geo%species_time_dependent

    POP_SUB(geometry_species_time_dependent)
  end function geometry_species_time_dependent

  ! ---------------------------------------------------------
  subroutine geometry_val_charge(geo, val_charge)
    type(geometry_t), intent(in) :: geo
    FLOAT,            intent(out) :: val_charge

    integer :: iatom

    PUSH_SUB(geometry_val_charge)

    val_charge = M_ZERO
    do iatom = 1, geo%natoms
      val_charge = val_charge - species_zval(geo%atom(iatom)%species)
    end do

    POP_SUB(geometry_val_charge)
  end subroutine geometry_val_charge

  ! ---------------------------------------------------------
  subroutine geometry_dipole(geo, dipole)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: dipole(:)

    integer :: ia

    PUSH_SUB(geometry_dipole)

    dipole(1:geo%space%dim) = M_ZERO
    do ia = 1, geo%natoms
      dipole(1:geo%space%dim) = dipole(1:geo%space%dim) + &
        species_zval(geo%atom(ia)%species)*geo%atom(ia)%x(1:geo%space%dim)
    end do
    dipole = P_PROTON_CHARGE*dipole

    POP_SUB(geometry_dipole)
  end subroutine geometry_dipole

  ! ---------------------------------------------------------
  function geometry_center_of_mass(geo, pseudo) result(pos)
    type(geometry_t),           intent(in) :: geo
    logical,          optional, intent(in) :: pseudo !< calculate center considering all species to have equal mass.
    FLOAT :: pos(MAX_DIM)

    FLOAT :: mass, total_mass
    integer :: ia

    PUSH_SUB(geometry_center_of_mass)

    pos = M_ZERO
    total_mass = M_ZERO
    mass = M_ONE
    do ia = 1, geo%natoms
      if (.not. optional_default(pseudo, .false.)) then
        mass = species_mass(geo%atom(ia)%species)
      end if
      pos(1:geo%space%dim) = pos(1:geo%space%dim) + mass*geo%atom(ia)%x(1:geo%space%dim)
      total_mass = total_mass + mass
    end do
    pos(1:geo%space%dim) = pos(1:geo%space%dim)/total_mass

    POP_SUB(geometry_center_of_mass)
  end function geometry_center_of_mass

  ! ---------------------------------------------------------
  function geometry_center_of_mass_vel(geo) result(vel)
    type(geometry_t), intent(in) :: geo
    FLOAT :: vel(MAX_DIM)

    FLOAT :: mass, total_mass
    integer :: iatom

    PUSH_SUB(geometry_center_of_mass_vel)

    vel = M_ZERO
    total_mass = M_ZERO
    do iatom = 1, geo%natoms
      mass = species_mass(geo%atom(iatom)%species)
      total_mass = total_mass + mass
      vel(1:geo%space%dim) = vel(1:geo%space%dim) + mass*geo%atom(iatom)%v(1:geo%space%dim)
    end do
    vel(1:geo%space%dim) = vel(1:geo%space%dim)/total_mass

    POP_SUB(geometry_center_of_mass_vel)
  end function geometry_center_of_mass_vel

  ! ---------------------------------------------------------
  function geometry_center(geo) result(pos)
    type(geometry_t), intent(in) :: geo
    FLOAT :: pos(MAX_DIM)

    FLOAT :: xmin(MAX_DIM), xmax(MAX_DIM)
    integer  :: i, j

    PUSH_SUB(geometry_center)

    xmin =  CNST(1e10)
    xmax = -CNST(1e10)
    do i = 1, geo%natoms
      do j = 1, geo%space%dim
        if (geo%atom(i)%x(j) > xmax(j)) xmax(j) = geo%atom(i)%x(j)
        if (geo%atom(i)%x(j) < xmin(j)) xmin(j) = geo%atom(i)%x(j)
      end do
    end do

    pos = M_ZERO
    pos(1:geo%space%dim) = (xmax(1:geo%space%dim) + xmin(1:geo%space%dim))/M_TWO

    POP_SUB(geometry_center)
  end function geometry_center

  ! ---------------------------------------------------------
  subroutine geometry_axis_large(geo, x, x2)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: x(MAX_DIM), x2(MAX_DIM)

    integer  :: i, j
    FLOAT :: rmax, r, r2

    PUSH_SUB(geometry_axis_large)

    ! first get the further apart atoms
    rmax = -CNST(1e10)
    do i = 1, geo%natoms
      do j = 1, geo%natoms/2 + 1
        r = sqrt(sum((geo%atom(i)%x(1:geo%space%dim)-geo%atom(j)%x(1:geo%space%dim))**2))
        if (r > rmax) then
          rmax = r
          x = geo%atom(i)%x - geo%atom(j)%x
        end if
      end do
    end do
    x  = x /sqrt(sum(x(1:geo%space%dim)**2))

    ! now let us find out what is the second most important axis
    rmax = -CNST(1e10)
    do i = 1, geo%natoms
      r2 = sum(x(1:geo%space%dim) * geo%atom(i)%x(1:geo%space%dim))
      r = sqrt(sum((geo%atom(i)%x(1:geo%space%dim) - r2*x(1:geo%space%dim))**2))
      if (r > rmax) then
        rmax = r
        x2 = geo%atom(i)%x - r2*x
      end if
    end do

    POP_SUB(geometry_axis_large)
  end subroutine geometry_axis_large

  ! ---------------------------------------------------------
  !> This subroutine assumes that the origin of the coordinates is the
  !! center of mass of the system
  subroutine geometry_axis_inertia(geo, x, x2, pseudo)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: x(MAX_DIM), x2(MAX_DIM)
    logical,          intent(in)  :: pseudo !< calculate axis considering all species to have equal mass.

    FLOAT :: mass, tinertia(MAX_DIM, MAX_DIM), eigenvalues(MAX_DIM)
    integer :: ii, jj, iatom
    type(unit_t) :: unit

    PUSH_SUB(geometry_axis_inertia)

    ! first calculate the inertia tensor
    tinertia = M_ZERO
    mass = M_ONE
    do iatom = 1, geo%natoms
      if (.not.pseudo) mass = species_mass(geo%atom(iatom)%species)
      do ii = 1, geo%space%dim
        do jj = 1, geo%space%dim
          tinertia(ii, jj) = tinertia(ii, jj) - mass*geo%atom(iatom)%x(ii)*geo%atom(iatom)%x(jj)
        end do
        tinertia(ii, ii) = tinertia(ii, ii) + mass*sum(geo%atom(iatom)%x(:)**2)
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
    call output_tensor(stdout, tinertia, geo%space%dim, unit, write_average = .true.)

    call lalg_eigensolve(geo%space%dim, tinertia, eigenvalues)

    write(message(1),'(a,6f25.6)') 'Eigenvalues: ', &
      (units_from_atomic(unit, eigenvalues(jj)), jj = 1, geo%space%dim)
    call messages_info(1)

    ! make a choice to fix the sign of the axis.
    do ii = 1, 2
      jj = maxloc(abs(tinertia(:,ii)), dim = 1)
      if (tinertia(jj,ii) < M_ZERO) tinertia(:,ii) = -tinertia(:,ii)
    end do
    x  = tinertia(:,1)
    x2 = tinertia(:,2)

    POP_SUB(geometry_axis_inertia)
  end subroutine geometry_axis_inertia

  ! ---------------------------------------------------------
  subroutine geometry_translate(geo, x)
    type(geometry_t), intent(inout) :: geo
    FLOAT,            intent(in)    :: x(MAX_DIM)

    integer  :: iatom

    PUSH_SUB(geometry_translate)

    do iatom = 1, geo%natoms
      geo%atom(iatom)%x(1:geo%space%dim) = geo%atom(iatom)%x(1:geo%space%dim) - x(1:geo%space%dim)
    end do
    do iatom = 1, geo%ncatoms
      geo%catom(iatom)%x(1:geo%space%dim) = geo%catom(iatom)%x(1:geo%space%dim) - x(1:geo%space%dim)
    end do

    POP_SUB(geometry_translate)
  end subroutine geometry_translate

  ! ---------------------------------------------------------
  subroutine geometry_rotate(geo, namespace, from, from2, to)
    type(geometry_t),  intent(inout) :: geo
    type(namespace_t), intent(in)    :: namespace
    FLOAT,             intent(in)    :: from(MAX_DIM)   !< assumed to be normalized
    FLOAT,             intent(in)    :: from2(MAX_DIM)  !< assumed to be normalized
    FLOAT,             intent(in)    :: to(MAX_DIM)     !< assumed to be normalized

    integer :: iatom, idim
    FLOAT :: m1(MAX_DIM, MAX_DIM), m2(MAX_DIM, MAX_DIM)
    FLOAT :: m3(MAX_DIM, MAX_DIM), f2(MAX_DIM), per(MAX_DIM)
    FLOAT :: alpha, r

    PUSH_SUB(geometry_rotate)

    if (geo%space%dim /= 3) then
      call messages_not_implemented("geometry_rotate in other than 3 dimensions", namespace=namespace)
    end if

    ! initialize matrices
    m1 = M_ZERO
    do idim = 1, MAX_DIM
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
    r = sqrt(sum(per(1:3)**2))
    if (r > M_ZERO) then
      per(1:3) = per(1:3)/r
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
    do iatom = 1, geo%natoms
      f2 = geo%atom(iatom)%x
      geo%atom(iatom)%x = matmul(m1, f2)
    end do

    do iatom = 1, geo%ncatoms
      f2 = geo%catom(iatom)%x
      geo%catom(iatom)%x = matmul(m1, f2)
    end do

    POP_SUB(geometry_rotate)
  contains

    ! ---------------------------------------------------------
    subroutine rotate(m, angle, dir)
      FLOAT,   intent(inout) :: m(MAX_DIM, MAX_DIM)
      FLOAT,   intent(in)    :: angle
      integer, intent(in)    :: dir

      FLOAT :: aux(MAX_DIM, MAX_DIM), ca, sa

      PUSH_SUB(geometry_rotate.rotate)

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

      POP_SUB(rotate)
    end subroutine rotate

  end subroutine geometry_rotate

  ! ---------------------------------------------------------
  subroutine geometry_grid_defaults(geo, def_h, def_rsize)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: def_h, def_rsize

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

  ! ---------------------------------------------------------
  subroutine geometry_grid_defaults_info(geo)
    type(geometry_t), intent(in) :: geo

    integer :: ispec

    PUSH_SUB(geometry_grid_defaults_info)

    do ispec = 1, geo%nspecies
      call messages_write("Species '"//trim(species_label(geo%species(ispec)))//"': spacing = ")
      if (species_def_h(geo%species(ispec)) > CNST(0.0)) then
        call messages_write(species_def_h(geo%species(ispec)), fmt = '(f7.3)')
        call messages_write(" b")
      else
        call messages_write(" unknown")
      end if
      call messages_write(", radius = ")
      if (species_def_rsize(geo%species(ispec)) > CNST(0.0)) then
        call messages_write(species_def_rsize(geo%species(ispec)), fmt = '(f5.1)')
        call messages_write(" b.")
      else
        call messages_write(" unknown.")
      end if
      call messages_info()
    end do

    POP_SUB(geometry_grid_defaults_info)
  end subroutine geometry_grid_defaults_info

  ! ---------------------------------------------------------
  subroutine geometry_get_positions(geo, q)
    type(geometry_t), intent(in)    :: geo
    FLOAT,            intent(inout) :: q(:, :)

    integer :: iatom

    PUSH_SUB(geometry_get_positions)

    do iatom = 1, geo%natoms
      q(iatom, 1:geo%space%dim) = geo%atom(iatom)%x(1:geo%space%dim)
    end do

    POP_SUB(geometry_get_positions)
  end subroutine geometry_get_positions

  ! ---------------------------------------------------------
  subroutine geometry_set_positions(geo, q)
    type(geometry_t), intent(inout) :: geo
    FLOAT,            intent(in)    :: q(:, :)

    integer :: iatom

    PUSH_SUB(geometry_get_positions)

    do iatom = 1, geo%natoms
      geo%atom(iatom)%x(1:geo%space%dim) = q(iatom, 1:geo%space%dim)
    end do

    POP_SUB(geometry_get_positions)
  end subroutine geometry_set_positions

  ! ---------------------------------------------------------
  subroutine geometry_end(geo)
    type(geometry_t), intent(inout) :: geo

    PUSH_SUB(geometry_end)

    call distributed_end(geo%atoms_dist)

    SAFE_DEALLOCATE_P(geo%atom)
    geo%natoms=0
    SAFE_DEALLOCATE_P(geo%catom)
    geo%ncatoms=0

    call species_end(geo%nspecies, geo%species)
    SAFE_DEALLOCATE_P(geo%species)
    geo%nspecies=0

    call species_end_global()

    POP_SUB(geometry_end)
  end subroutine geometry_end

end module geometry_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
