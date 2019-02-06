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
  use loct_pointer_oct_m
  use loct_math_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use read_coords_oct_m
  use space_oct_m
  use species_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                          &
    geometry_t,                      &
    geometry_nullify,                &
    geometry_init,                   &
    geometry_init_xyz,               &
    geometry_init_species,           &
    geometry_partition,              &
    geometry_copy,                   &
    geometry_end,                    &
    geometry_dipole,                 &
    geometry_min_distance,           &
    geometry_mass,                   &
    cm_pos,                          &
    cm_vel,                          &
    geometry_write_xyz,              &
    geometry_read_xyz,               &
    geometry_val_charge,             &
    geometry_grid_defaults,          &
    geometry_grid_defaults_info,     &
    geometry_species_time_dependent, &
    geometry_get_positions,          &
    geometry_set_positions

  type geometry_t
    type(space_t), pointer :: space
    integer                :: natoms
    type(atom_t), pointer  :: atom(:)
    integer :: nspecies
    type(species_t), pointer :: species(:)
    logical :: species_time_dependent !< For time-dependent user defined species
    FLOAT :: kinetic_energy         !< the ion kinetic energy
    logical :: nlpp                 !< does any species have non-local pp?
    logical :: nlcc                 !< does any species have non-local core corrections?
    type(distributed_t) :: atoms_dist
    logical          :: reduced_coordinates !< If true the coordinates are stored in
                                            !! reduced coordinates and need to be converted.
    !> variables for passing info from XSF input to simul_box_init
    integer :: periodic_dim
    FLOAT :: lsize(MAX_DIM)
  end type geometry_t

contains

  ! ---------------------------------------------------------
  subroutine geometry_nullify(this)
    type(geometry_t), intent(inout) :: this

    PUSH_SUB(geometry_nullify)

    nullify(this%space, this%atom, this%species)
    this%natoms=0
    this%nspecies=0
    this%species_time_dependent=.false.
    this%kinetic_energy=M_ZERO
    this%nlpp=.false.
    this%nlcc=.false.
    call distributed_nullify(this%atoms_dist, 0)
    this%reduced_coordinates=.false.
    this%periodic_dim=0
    this%lsize=M_ZERO

    POP_SUB(geometry_nullify)
  end subroutine geometry_nullify

  ! ---------------------------------------------------------
  subroutine geometry_init(geo, space, print_info)
    type(geometry_t),           intent(inout) :: geo
    type(space_t),    target,   intent(in)    :: space
    logical,          optional, intent(in)    :: print_info

    PUSH_SUB(geometry_init)

    geo%space => space

    call species_init_global()
    
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

    type(read_coords_info) :: xyz
    integer             :: ia
    logical             :: move

    PUSH_SUB(geometry_init_xyz)

    call read_coords_init(xyz)

    ! load positions of the atoms
    call read_coords_read('Coordinates', xyz, geo%space)

    if(xyz%n < 1) then
      message(1) = "Coordinates have not been defined."
      call messages_fatal(1)
    end if

    ! copy information from xyz to geo
    geo%natoms = xyz%n
    nullify(geo%atom)
    if(geo%natoms>0)then
      SAFE_ALLOCATE(geo%atom(1:geo%natoms))
      do ia = 1, geo%natoms
        move=.true.
        if(bitand(xyz%flags, XYZ_FLAGS_MOVE) /= 0) move=xyz%atom(ia)%move
        call atom_init(geo%atom(ia), xyz%atom(ia)%label, xyz%atom(ia)%x, move=move)
      end do
    end if

    geo%reduced_coordinates = xyz%source == READ_COORDS_REDUCED
    geo%periodic_dim = xyz%periodic_dim
    geo%lsize(:) = xyz%lsize(:)

    call read_coords_end(xyz)

    POP_SUB(geometry_init_xyz)
  end subroutine geometry_init_xyz


  ! ---------------------------------------------------------
  subroutine geometry_init_species(geo, print_info)
    type(geometry_t),  intent(inout) :: geo
    logical, optional, intent(in)    :: print_info

    logical :: print_info_
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
    atoms2: do i = 1, geo%natoms
      do j = 1, i - 1
        if(atom_same_species(geo%atom(j), geo%atom(i))) cycle atoms2
      end do
      k = k + 1
      call species_init(geo%species(k), atom_get_label(geo%atom(j)), k)
      call species_read(geo%species(k))
    end do atoms2

    ! Reads the spin components. This is read here, as well as in states_init,
    ! to be able to pass it to the pseudopotential initializations subroutine.
    call parse_variable('SpinComponents', 1, ispin)
    if(.not.varinfo_valid_option('SpinComponents', ispin)) call messages_input_error('SpinComponents')
    ispin = min(2, ispin)

    if(print_info_) then
      call messages_print_stress(stdout, "Species")
    end if
    do i = 1, geo%nspecies
      call species_build(geo%species(i), ispin, geo%space%dim, print_info=print_info_)
    end do
    if(print_info_) then
      call messages_print_stress(stdout)
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

    POP_SUB(geometry_init_species)
  end subroutine geometry_init_species

  ! ---------------------------------------------------------

  subroutine geometry_partition(geo, mc)
    type(geometry_t),            intent(inout) :: geo
    type(multicomm_t),           intent(in)    :: mc

    PUSH_SUB(geometry_partition)

    call distributed_init(geo%atoms_dist, geo%natoms, mc%group_comm(P_STRATEGY_STATES), "atoms")

    POP_SUB(geometry_partition)
  end subroutine geometry_partition


  ! ---------------------------------------------------------
  subroutine geometry_end(geo)
    type(geometry_t), intent(inout) :: geo

    PUSH_SUB(geometry_end)

    call distributed_end(geo%atoms_dist)

    SAFE_DEALLOCATE_P(geo%atom)
    geo%natoms=0

    call species_end(geo%nspecies, geo%species)
    SAFE_DEALLOCATE_P(geo%species)
    geo%nspecies=0

    call species_end_global()
    
    POP_SUB(geometry_end)
  end subroutine geometry_end

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
        species_zval(geo%atom(ia)%species)*geo%atom(ia)%x(1:geo%space%dim)
    end do
    dipole = P_PROTON_CHARGE*dipole

    POP_SUB(geometry_dipole)
  end subroutine geometry_dipole

  !> Beware: this is wrong for periodic systems. Use simul_box_min_distance instead.
  ! ---------------------------------------------------------
  FLOAT function geometry_min_distance(geo) result(rmin)
    type(geometry_t),  intent(in) :: geo

    integer :: i, j
    FLOAT   :: r
    type(species_t), pointer :: species

    PUSH_SUB(geometry_min_distance)

    rmin = huge(rmin)
    do i = 1, geo%natoms
      call atom_get_species(geo%atom(i), species)
      do j = i + 1, geo%natoms
        call atom_get_species(geo%atom(i), species)
        r = atom_distance(geo%atom(i), geo%atom(j))
        if(r < rmin) then
          rmin = r
        end if
      end do
    end do

    nullify(species)

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
      pos = pos + species_mass(geo%atom(ia)%species) * geo%atom(ia)%x
      mass = mass + species_mass(geo%atom(ia)%species)
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
      vel = vel + species_mass(geo%atom(iatom)%species) * geo%atom(iatom)%v
      mass = mass + species_mass(geo%atom(iatom)%species)
    end do
    vel = vel / mass

    POP_SUB(cm_vel)
  end subroutine cm_vel


  ! ---------------------------------------------------------
  subroutine geometry_write_xyz(geo, fname, append, comment)
    type(geometry_t),    intent(in) :: geo
    character(len=*),    intent(in) :: fname
    logical,             intent(in), optional :: append
    character(len=*),    intent(in), optional :: comment

    integer :: iatom, iunit
    character(len=6) position

    if( .not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(atom_write_xyz)

    position = 'asis'
    if(present(append)) then
      if(append) position = 'append'
    end if
    iunit = io_open(trim(fname)//'.xyz', action='write', position=position)

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

    POP_SUB(atom_write_xyz)
  end subroutine geometry_write_xyz

  ! ---------------------------------------------------------
  subroutine geometry_read_xyz(geo, fname, comment)
    type(geometry_t),    intent(inout) :: geo
    character(len=*),    intent(in) :: fname
    character(len=*),    intent(in), optional :: comment

    integer :: iatom, iunit

    PUSH_SUB(geometry_read_xyz)

    iunit = io_open(trim(fname)//'.xyz', action='read', position='rewind')

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

    POP_SUB(geometry_read_xyz)
  end subroutine geometry_read_xyz

  ! ---------------------------------------------------------

  subroutine geometry_val_charge(geo, val_charge)
    type(geometry_t), intent(in) :: geo
    FLOAT,           intent(out) :: val_charge

    integer :: iatom

    PUSH_SUB(geometry_val_charge)

    val_charge = M_ZERO
    do iatom = 1, geo%natoms
      val_charge = val_charge - species_zval(geo%atom(iatom)%species)
    end do

    POP_SUB(geometry_val_charge)
  end subroutine geometry_val_charge


  ! ---------------------------------------------------------

  FLOAT pure function geometry_mass(geo) result(mass)
    type(geometry_t), intent(in) :: geo

    integer :: iatom

    mass = M_ZERO
    do iatom = 1, geo%natoms
      mass = mass + species_mass(geo%atom(iatom)%species)
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

  ! ---------------------------------------------------------

  subroutine geometry_grid_defaults_info(geo)
    type(geometry_t), intent(in) :: geo
  
    integer :: ispec

    PUSH_SUB(geometry_grid_defaults_info)

    do ispec = 1, geo%nspecies
      call messages_write("Species '"//trim(species_label(geo%species(ispec)))//"': spacing = ")
      if(species_def_h(geo%species(ispec)) > CNST(0.0)) then
        call messages_write(species_def_h(geo%species(ispec)), fmt = '(f7.3)')
        call messages_write(" b")
      else
        call messages_write(" unknown")
      end if
      call messages_write(", radius = ")
      if(species_def_rsize(geo%species(ispec)) > CNST(0.0)) then
        call messages_write(species_def_rsize(geo%species(ispec)), fmt = '(f5.1)')
        call messages_write(" b.")
      else
        call messages_write(" unknown.")
      end if
      call messages_info()
    end do

    POP_SUB(geometry_grid_defaults_info)
  end subroutine geometry_grid_defaults_info

  !--------------------------------------------------------------
  subroutine geometry_copy(geo_out, geo_in)
    type(geometry_t), intent(out) :: geo_out
    type(geometry_t), intent(in)  :: geo_in

    PUSH_SUB(geometry_copy)

    geo_out%natoms = geo_in%natoms
    SAFE_ALLOCATE(geo_out%atom(1:geo_out%natoms))
    geo_out%atom = geo_in%atom

    geo_out%nspecies = geo_in%nspecies
    SAFE_ALLOCATE(geo_out%species(1:geo_out%nspecies))
    geo_out%species = geo_in%species

    geo_out%kinetic_energy    = geo_in%kinetic_energy
    geo_out%nlpp              = geo_in%nlpp
    geo_out%nlcc              = geo_in%nlcc

    call distributed_copy(geo_in%atoms_dist, geo_out%atoms_dist)

    POP_SUB(geometry_copy)
  end subroutine geometry_copy

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


end module geometry_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
