#include "global.h"

module igeometry_m

  use global_m
  use messages_m
  use profiling_m

  use iatom_m,       only: atom_init, atom_create_data_object, atom_classical_create_data_object
  use distributed_m, only: distributed_nullify
  use json_m,        only: JSON_OK, json_object_t, json_array_t, json_init, json_get, json_set, json_append
  use space_m,       only: space_t
  use species_m,     only: LABEL_LEN, species_label, species_init_from_data_object, species_create_data_object

  use geometry_m, only: &
    geometry_t,         &
    geometry_copy,      &
    geometry_end

  implicit none

  private
  public ::                      &
    geometry_t,                  &
    geometry_init,               &
    geometry_create_data_object, &
    geometry_copy,               &
    geometry_end

contains

  ! ---------------------------------------------------------
  subroutine geometry_init(this, space, json)
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
    call json_get(json, "nspecies", this%nspecies, ierr=ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "nspecies" from geometry data object.'
      call messages_fatal(1)
      return
    end if
    SAFE_ALLOCATE(this%species(this%nspecies))
    call json_get(json, "species", species, ierr=ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "species" array from geometry data object.'
      call messages_fatal(1)
      return
    end if
    do i=1, this%nspecies
      call json_get(species, i, spec, ierr=ierr)
      if(ierr/=JSON_OK)then
        write(unit=message(1), fmt="(a,i3,a)") &
          'Could not read the ', i, 'th "species" element from geometry data object.'
        call messages_fatal(1)
        return
      end if
      call species_init_from_data_object(this%species(i), i, spec)
    end do
    call json_get(json, "natoms", this%natoms, ierr=ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "natoms" from geometry data object.'
      call messages_fatal(1)
      return
    end if
    SAFE_ALLOCATE(this%atom(this%natoms))
    call json_get(json, "atom", atoms, ierr=ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "atom" array from geometry data object.'
      call messages_fatal(1)
      return
    end if
    do i=1, this%natoms
      call json_get(atoms, i, atom, ierr=ierr)
      if(ierr/=JSON_OK)then
        write(unit=message(1), fmt="(a,i3,a)") &
          'Could not read the ', i, 'th "atom" element from geometry data object.'
        call messages_fatal(1)
        return
      end if
      do j=1, this%nspecies
        call json_get(atom, "label", label, ierr=ierr)
        if(ierr/=JSON_OK)then
          write(unit=message(1), fmt="(a,i3,a)") &
            'Could not read the ', i, 'th "atom" element "label" from geometry data object.'
          call messages_fatal(1)
          return
        end if
        if(trim(label)==trim(species_label(this%species(j))))then
          call atom_init(this%atom(i), this%species(j), atom)
          exit
        end if
      end do
    end do
    call json_get(json, "ncatoms", this%ncatoms, ierr=ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "ncatoms" from geometry data object.'
      call messages_fatal(1)
      return
    end if
    SAFE_ALLOCATE(this%catom(this%ncatoms))
    call json_get(json, "catom", atoms, ierr=ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "catom" array from geometry data object.'
      call messages_fatal(1)
      return
    end if
    do i=1, this%ncatoms
      call json_get(atoms, i, atom, ierr=ierr)
      if(ierr/=JSON_OK)then
        write(unit=message(1), fmt="(a,i3,a)") &
          'Could not read the ', i, 'th "catom" from geometry data object.'
        call messages_fatal(1)
        return
      end if
      call atom_init(this%catom(i), atom)
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
  end subroutine geometry_init

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

end module igeometry_m

!! Local Variables:
!! mode: f90
!! End:
