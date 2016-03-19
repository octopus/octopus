#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME
#undef LIST_INCLUDE_PREFIX
#undef LIST_INCLUDE_HEADER
#undef LIST_INCLUDE_BODY

#define LIST_TEMPLATE_NAME atom
#include "tlist_inc.F90"
#undef LIST_TEMPLATE_NAME

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#define DICT_TEMPLATE_NAME species
#include "tdict_inc.F90"
#undef DICT_TEMPLATE_NAME

module geo_build_oct_m

  use atom_oct_m
  use atom_list_oct_m
  use basis_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use space_oct_m
  use species_oct_m
  use species_dict_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::             &
    GEO_BUILD_NAME_LEN

  public ::          &
    GEO_BUILD_OK,    &
    GEO_BUILD_ERROR

  public ::      &
    geo_build_t

  public ::           &
    geo_build_len,    &
    geo_build_init,   &
    geo_build_extend, &
    geo_build_export, &
    geo_build_get,    &
    geo_build_copy,   &
    geo_build_end

  integer, parameter :: GEO_BUILD_NAME_LEN = LABEL_LEN

  integer, parameter :: GEO_BUILD_OK    = 0
  integer, parameter :: GEO_BUILD_ERROR =-1

  type :: geo_build_t
    private
    type(space_t), pointer :: space  =>null()
    type(atom_list_t)      :: list
    type(species_dict_t)   :: dict
  end type geo_build_t

  type :: geo_build_iterator_t
    private
    type(geo_build_t),    pointer :: self =>null()
    type(atom_list_iterator_t)    :: aitr
    type(species_dict_iterator_t) :: sitr
  end type geo_build_iterator_t

  interface geo_build__new__
    module procedure geo_build__new__species
    module procedure geo_build__new__atom
  end interface geo_build__new__

  interface geo_build__idel__
    module procedure geo_build__idel__species
    module procedure geo_build__idel__atom
  end interface geo_build__idel__

  interface geo_build__del__
    module procedure geo_build__del__species
    module procedure geo_build__del__atom
  end interface geo_build__del__

  interface geo_build__extend__
    module procedure geo_build__extend__geo
    module procedure geo_build__extend__geometry
    module procedure geo_build__extend__config
  end interface geo_build__extend__

  interface geo_build_init
    module procedure geo_build_init_type
    module procedure geo_build_init_copy
  end interface geo_build_init

  interface geo_build_extend
    module procedure geo_build_extend_geo
    module procedure geo_build_extend_geometry
    module procedure geo_build_extend_config
    module procedure geo_build_extend_geo_config
    module procedure geo_build_extend_geometry_config
    module procedure geo_build_extend_config_config
    module procedure geo_build_extend_geo_list
    module procedure geo_build_extend_geometry_list
    module procedure geo_build_extend_config_list
  end interface geo_build_extend

  interface geo_build_get
    module procedure geo_build_get_space
  end interface geo_build_get

  interface geo_build_iterator_next
    module procedure geo_build_iterator_next_atom
    module procedure geo_build_iterator_next_species
  end interface geo_build_iterator_next

contains

  ! ---------------------------------------------------------
  subroutine geo_build__new__species(this, label, spec)
    type(geo_build_t), intent(inout) :: this
    character(len=*),  intent(in)    :: label
    type(species_t),  pointer        :: spec

    PUSH_SUB(geo_build__new__species)

    nullify(spec)
    SAFE_ALLOCATE(spec)
    call species_dict_set(this%dict, trim(adjustl(label)), spec)


    POP_SUB(geo_build__new__species)
  end subroutine geo_build__new__species

  ! ---------------------------------------------------------
  subroutine geo_build__idel__species(spec)
    type(species_t), pointer :: spec

    PUSH_SUB(geo_build__idel__species)

    if(associated(spec))then
      call species_end(spec)
      SAFE_DEALLOCATE_P(spec)
    end if
    nullify(spec)

    POP_SUB(geo_build__idel__species)
  end subroutine geo_build__idel__species

  ! ---------------------------------------------------------
  subroutine geo_build__del__species(this, spec)
    type(geo_build_t), intent(inout) :: this
    type(species_t),  pointer        :: spec

    type(species_t), pointer :: ispc
    integer                  :: ierr

    PUSH_SUB(geo_build__del__species)

    nullify(ispc)
    if(associated(spec))then
      call species_dict_del(this%dict, species_label(spec), ispc, ierr)
      ASSERT(ierr==SPECIES_DICT_OK)
      ASSERT(associated(ispc))
      ASSERT(associated(ispc,spec))
      nullify(ispc)
      call geo_build__idel__(spec)
    end if
    nullify(spec)

    POP_SUB(geo_build__del__species)
  end subroutine geo_build__del__species

  ! ---------------------------------------------------------
  subroutine geo_build__new__atom(this, atom)
    type(geo_build_t), intent(inout) :: this
    type(atom_t),     pointer        :: atom

    PUSH_SUB(geo_build__new__atom)

    nullify(atom)
    SAFE_ALLOCATE(atom)
    call atom_list_append(this%list, atom)

    POP_SUB(geo_build__new__atom)
  end subroutine geo_build__new__atom

  ! ---------------------------------------------------------
  subroutine geo_build__idel__atom(atom)
    type(atom_t), pointer :: atom

    PUSH_SUB(geo_build__idel__atom)

    if(associated(atom))then
      call atom_end(atom)
      SAFE_DEALLOCATE_P(atom)
    end if
    nullify(atom)

    POP_SUB(geo_build__idel__atom)
  end subroutine geo_build__idel__atom

  ! ---------------------------------------------------------
  subroutine geo_build__del__atom(this, atom)
    type(geo_build_t), intent(inout) :: this
    type(atom_t),     pointer        :: atom

    type(atom_t), pointer :: iatm
    integer               :: ierr

    PUSH_SUB(geo_build__del__atom)

    nullify(iatm)
    if(associated(atom))then
      call atom_list_del(this%list, atom, ierr)
      ASSERT(ierr==ATOM_LIST_OK)
      call geo_build__idel__(atom)
    end if
    nullify(atom)

    POP_SUB(geo_build__del__atom)
  end subroutine geo_build__del__atom

  ! ---------------------------------------------------------
  subroutine geo_build_iadd_species_from_config(this, config)
    type(geo_build_t),   intent(inout) :: this
    type(json_object_t), intent(in)    :: config

    character(len=GEO_BUILD_NAME_LEN) :: label
    type(species_t),          pointer :: spec
    integer                           :: ierr

    PUSH_SUB(geo_build_iadd_species_from_config)

    nullify(spec)
    call json_get(config, "label", label, ierr)
    ASSERT(ierr==JSON_OK)
    if(.not.species_dict_has_key(this%dict, trim(adjustl(label))))then
      call geo_build__new__(this, label, spec)
      call species_init_from_data_object(spec, 0, config)
      nullify(spec)
    end if

    POP_SUB(geo_build_iadd_species_from_config)
  end subroutine geo_build_iadd_species_from_config

  ! ---------------------------------------------------------
  subroutine geo_build_iadd_species_from_species(this, that)
    type(geo_build_t), intent(inout) :: this
    type(species_t),   intent(in)    :: that

    type(species_t), pointer :: spec

    PUSH_SUB(geo_build_iadd_species_from_species)

    nullify(spec)
    if(.not.species_dict_has_key(this%dict, species_label(that)))then
      call geo_build__new__(this, species_label(that), spec)
      call species_init(spec, "", 0)
      call species_copy(spec, that)
      nullify(spec)
    end if

    POP_SUB(geo_build_iadd_species_from_species)
  end subroutine geo_build_iadd_species_from_species

  ! ---------------------------------------------------------
  subroutine geo_build_iadd_atom_from_config(this, config, basis)
    type(geo_build_t),       intent(inout) :: this
    type(json_object_t),     intent(in)    :: config
    type(basis_t), optional, intent(in)    :: basis

    character(len=GEO_BUILD_NAME_LEN) :: label
    type(atom_t),             pointer :: atom
    type(species_t),          pointer :: spec
    integer                           :: ierr

    PUSH_SUB(geo_build_iadd_atom_from_config)

    nullify(atom, spec)
    call json_get(config, "label", label, ierr)
    ASSERT(ierr==JSON_OK)
    call species_dict_get(this%dict, label, spec, ierr)
    ASSERT(ierr==SPECIES_DICT_OK)
    ASSERT(associated(spec))
    call geo_build__new__(this, atom)
    call atom_init_from_data_object(atom, spec, config)
    if(present(basis)) call basis_to_external(basis, atom%x)

    POP_SUB(geo_build_iadd_atom_from_config)
  end subroutine geo_build_iadd_atom_from_config

  ! ---------------------------------------------------------
  subroutine geo_build_iadd_atom_from_atom(this, that, basis)
    type(geo_build_t),       intent(inout) :: this
    type(atom_t),            intent(in)    :: that
    type(basis_t), optional, intent(in)    :: basis

    real(kind=wp), dimension(MAX_DIM) :: x
    type(atom_t),             pointer :: atom
    type(species_t),          pointer :: spec
    integer                           :: ierr

    PUSH_SUB(geo_build_iadd_atom_from_atom)

    nullify(atom, spec)
    call species_dict_get(this%dict, that%label, spec, ierr)
    ASSERT(ierr==SPECIES_DICT_OK)
    ASSERT(associated(spec))
    call geo_build__new__(this, atom)
    if(present(basis))then
      call basis_to_external(basis, that%x, x)
    else
      x = that%x
    end if
    call atom_init(atom, trim(adjustl(that%label)), x, species=spec)
 
    POP_SUB(geo_build_iadd_atom_from_atom)
  end subroutine geo_build_iadd_atom_from_atom

  ! ---------------------------------------------------------
  subroutine geo_build__extend__geo(this, that, basis)
    type(geo_build_t),       intent(inout) :: this
    type(geo_build_t),       intent(in)    :: that
    type(basis_t), optional, intent(in)    :: basis

    type(geo_build_iterator_t) :: iter
    type(atom_t),      pointer :: atom
    type(species_t),   pointer :: spec
    integer                    :: ierr

    PUSH_SUB(geo_build__extend__geo)

    ASSERT(this%space==that%space)
    call geo_build_iterator_init(iter, that)
    do
      nullify(spec)
      call geo_build_iterator_next(iter, spec, ierr)
      if(ierr/=GEO_BUILD_OK)exit
      ASSERT(associated(spec))
      call geo_build_iadd_species_from_species(this, spec)
    end do
    nullify(spec)
    do
      nullify(atom)
      call geo_build_iterator_next(iter, atom, ierr)
      if(ierr/=GEO_BUILD_OK)exit
      ASSERT(associated(atom))
      call geo_build_iadd_atom_from_atom(this, atom, basis)
    end do
    nullify(atom)
    call geo_build_iterator_end(iter)

    POP_SUB(geo_build__extend__geo)
  end subroutine geo_build__extend__geo

  ! ---------------------------------------------------------
  subroutine geo_build__extend__geometry(this, that, basis)
    type(geo_build_t),       intent(inout) :: this
    type(geometry_t),        intent(in)    :: that
    type(basis_t), optional, intent(in)    :: basis

    integer :: indx

    PUSH_SUB(geo_build__extend__geometry)

    ASSERT(this%space==that%space)
    do indx = 1, that%nspecies
      call geo_build_iadd_species_from_species(this, that%species(indx))
    end do
    do indx = 1, that%natoms
      call geo_build_iadd_atom_from_atom(this, that%atom(indx), basis)
    end do

    POP_SUB(geo_build__extend__geometry)
  end subroutine geo_build__extend__geometry

  ! ---------------------------------------------------------
  subroutine geo_build__extend__config(this, config, basis)
    type(geo_build_t),       intent(inout) :: this
    type(json_object_t),     intent(in)    :: config
    type(basis_t), optional, intent(in)    :: basis

    type(json_array_iterator_t)  :: iter
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    integer                      :: nitm, ierr

    PUSH_SUB(geo_build__extend__config)

    nullify(cnfg, list)
    call json_get(config, "nspecies", nitm, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "species", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(json_len(list)==nitm)
    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call geo_build_iadd_species_from_config(this, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list)
    call json_get(config, "natoms", nitm, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "atom", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(json_len(list)==nitm)
    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call geo_build_iadd_atom_from_config(this, cnfg, basis)
    end do
    call json_end(iter)
    nullify(cnfg, list)

    POP_SUB(geo_build__extend__config)
  end subroutine geo_build__extend__config

  ! ---------------------------------------------------------
  elemental function geo_build_len(this) result(len)
    type(geo_build_t), intent(in) :: this

    integer :: len

    len = atom_list_len(this%list)

  end function geo_build_len

  ! ---------------------------------------------------------
  subroutine geo_build_init_type(this, space)
    type(geo_build_t),     intent(out) :: this
    type(space_t), target, intent(in)  :: space

    PUSH_SUB(geo_build_init_type)

    this%space => space
    call atom_list_init(this%list)
    call species_dict_init(this%dict)

    POP_SUB(geo_build_init_type)
  end subroutine geo_build_init_type

  ! ---------------------------------------------------------
  subroutine geo_build_init_copy(this, that)
    type(geo_build_t), intent(out) :: this
    type(geo_build_t), intent(in)  :: that

    PUSH_SUB(geo_build_init_copy)

    ASSERT(associated(that%space))
    call geo_build_init(this, that%space)

    POP_SUB(geo_build_init_copy)
  end subroutine geo_build_init_copy

  ! ---------------------------------------------------------
  subroutine geo_build_extend_geo(this, that)
    type(geo_build_t), intent(inout) :: this
    type(geo_build_t), intent(in)    :: that

    PUSH_SUB(geo_build_extend_geo)

    call geo_build__extend__(this, that)

    POP_SUB(geo_build_extend_geo)
  end subroutine geo_build_extend_geo

  ! ---------------------------------------------------------
  subroutine geo_build_extend_geometry(this, that)
    type(geo_build_t), intent(inout) :: this
    type(geometry_t),  intent(in)    :: that

    PUSH_SUB(geo_build_extend_geometry)

    call geo_build__extend__(this, that)

    POP_SUB(geo_build_extend_geometry)
  end subroutine geo_build_extend_geometry

  ! ---------------------------------------------------------
  subroutine geo_build_extend_config(this, that)
    type(geo_build_t),   intent(inout) :: this
    type(json_object_t), intent(in)    :: that

    PUSH_SUB(geo_build_extend_config)

    call geo_build__extend__(this, that)

    POP_SUB(geo_build_extend_config)
  end subroutine geo_build_extend_config

  ! ---------------------------------------------------------
  subroutine geo_build_extend_geo_config(this, that, config)
    type(geo_build_t),   intent(inout) :: this
    type(geo_build_t),   intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    type(basis_t) :: base

    PUSH_SUB(geo_build_extend_geo_config)

    call basis_init(base, this%space, config)
    call geo_build__extend__(this, that, base)
    call basis_end(base)

    POP_SUB(geo_build_extend_geo_config)
  end subroutine geo_build_extend_geo_config

  ! ---------------------------------------------------------
  subroutine geo_build_extend_geometry_config(this, that, config)
    type(geo_build_t),   intent(inout) :: this
    type(geometry_t),    intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    type(basis_t) :: base

    PUSH_SUB(geo_build_extend_geometry_config)

    call basis_init(base, this%space, config)
    call geo_build__extend__(this, that, base)
    call basis_end(base)

    POP_SUB(geo_build_extend_geometry_config)
  end subroutine geo_build_extend_geometry_config

  ! ---------------------------------------------------------
  subroutine geo_build_extend_config_config(this, that, config)
    type(geo_build_t),   intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    type(basis_t) :: base

    PUSH_SUB(geo_build_extend_config_config)

    call basis_init(base, this%space, config)
    call geo_build__extend__(this, that, base)
    call basis_end(base)

    POP_SUB(geo_build_extend_config_config)
  end subroutine geo_build_extend_config_config

  ! ---------------------------------------------------------
  subroutine geo_build_extend_geo_list(this, that, list)
    type(geo_build_t),  intent(inout) :: this
    type(geo_build_t),  intent(in)    :: that
    type(json_array_t), intent(in)    :: list

    type(json_object_t), pointer :: cnfg
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr

    PUSH_SUB(geo_build_extend_geo_list)

    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call geo_build_extend(this, that, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg)

    POP_SUB(geo_build_extend_geo_list)
  end subroutine geo_build_extend_geo_list

  ! ---------------------------------------------------------
  subroutine geo_build_extend_geometry_list(this, that, list)
    type(geo_build_t),  intent(inout) :: this
    type(geometry_t),   intent(in)    :: that
    type(json_array_t), intent(in)    :: list

    type(json_object_t), pointer :: cnfg
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr

    PUSH_SUB(geo_build_extend_geometry_list)

    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call geo_build_extend(this, that, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg)

    POP_SUB(geo_build_extend_geometry_list)
  end subroutine geo_build_extend_geometry_list

  ! ---------------------------------------------------------
  subroutine geo_build_extend_config_list(this, that, list)
    type(geo_build_t),   intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    type(json_array_t),  intent(in)    :: list

    type(json_object_t), pointer :: cnfg
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr

    PUSH_SUB(geo_build_extend_config_list)

    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call geo_build_extend(this, that, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg)

    POP_SUB(geo_build_extend_config_list)
  end subroutine geo_build_extend_config_list

  ! ---------------------------------------------------------
  subroutine geo_build_export(this, that)
    type(geo_build_t), intent(in)  :: this
    type(geometry_t),  intent(out) :: that

    type(geo_build_iterator_t) :: iter
    type(species_dict_t)       :: sdct
    type(atom_t),      pointer :: atom
    type(species_t),   pointer :: spec
    integer                    :: indx, ierr

    PUSH_SUB(geo_build_export)

    call geometry_nullify(that)
    that%space => this%space
    that%nspecies = species_dict_len(this%dict)
    ASSERT(that%nspecies>0)
    SAFE_ALLOCATE(that%species(that%nspecies))
    call species_dict_init(sdct, that%nspecies)
    call geo_build_iterator_init(iter, this)
    do indx = 1, that%nspecies
      nullify(spec)
      call geo_build_iterator_next(iter, spec, ierr)
      ASSERT(ierr==SPECIES_DICT_OK)
      ASSERT(associated(spec))
      call species_init(that%species(indx), "", 0)
      call species_copy(that%species(indx), spec, indx)
      call species_dict_set(sdct, species_label(that%species(indx)), that%species(indx))
    end do
    call geo_build_iterator_end(iter)
    that%natoms = atom_list_len(this%list)
    ASSERT(that%natoms>0)
    SAFE_ALLOCATE(that%atom(that%natoms))
    call geo_build_iterator_init(iter, this)
    do indx = 1, that%natoms
      nullify(atom, spec)
      call geo_build_iterator_next(iter, atom, ierr)
      ASSERT(ierr==ATOM_LIST_OK)
      ASSERT(associated(atom))
      call species_dict_get(sdct, atom_get_label(atom), spec, ierr)
      ASSERT(ierr==SPECIES_DICT_OK)
      ASSERT(associated(spec))
      call atom_init(that%atom(indx), atom_get_label(atom), atom%x, species=spec)
    end do
    call geo_build_iterator_end(iter)
    call species_dict_end(sdct)
    nullify(atom, spec)

    POP_SUB(geo_build_export)
  end subroutine geo_build_export

  ! ---------------------------------------------------------
  subroutine geo_build_get_space(this, that)
    type(geo_build_t),   target, intent(in) :: this
    type(space_t), pointer             :: that

    PUSH_SUB(geo_build_get_space)

    nullify(that)
    if(associated(this%space)) that => this%space

    POP_SUB(geo_build_get_space)
  end subroutine geo_build_get_space

  ! ---------------------------------------------------------
  subroutine geo_build_copy(this, that)
    type(geo_build_t), intent(inout) :: this
    type(geo_build_t), intent(in)    :: that

    PUSH_SUB(geo_build_copy)

    call geo_build_end(this)
    if(associated(that%space))then
      call geo_build_init(this, that%space)
      call geo_build_extend(this, that)
    end if

    POP_SUB(geo_build_copy)
  end subroutine geo_build_copy

  ! ---------------------------------------------------------
  subroutine geo_build__end__list(this)
    type(atom_list_t), intent(inout) :: this

    type(atom_t), pointer :: atom

    PUSH_SUB(geo_build__end__list)

    do
      nullify(atom)
      call atom_list_pop(this, atom)
      if(.not.associated(atom))exit
      call geo_build__idel__(atom)
    end do
    nullify(atom)
    call atom_list_end(this)

    POP_SUB(geo_build__end__list)
  end subroutine geo_build__end__list

  ! ---------------------------------------------------------
  subroutine geo_build__end__dict(this)
    type(species_dict_t), intent(inout) :: this

    type(species_t), pointer :: spec

    PUSH_SUB(geo_build__end__dict)

    do
      nullify(spec)
      call species_dict_pop(this, spec)
      if(.not.associated(spec))exit
      call geo_build__idel__(spec)
    end do
    nullify(spec)
    call species_dict_end(this)

    POP_SUB(geo_build__end__dict)
  end subroutine geo_build__end__dict

  ! ---------------------------------------------------------
  subroutine geo_build_end(this)
    type(geo_build_t), target, intent(inout) :: this

    PUSH_SUB(geo_build_end)

    nullify(this%space)
    call geo_build__end__list(this%list)
    call geo_build__end__dict(this%dict)

    POP_SUB(geo_build_end)
  end subroutine geo_build_end

  ! ---------------------------------------------------------
  subroutine geo_build_iterator_init(this, that)
    type(geo_build_iterator_t), intent(out) :: this
    type(geo_build_t),  target, intent(in)  :: that

    PUSH_SUB(geo_build_iterator_init)

    this%self => that
    call atom_list_init(this%aitr, that%list)
    call species_dict_init(this%sitr, that%dict)

    POP_SUB(geo_build_iterator_init)
  end subroutine geo_build_iterator_init

  ! ---------------------------------------------------------
  subroutine geo_build_iterator_next_atom(this, atom, ierr)
    type(geo_build_iterator_t), intent(inout) :: this
    type(atom_t),              pointer        :: atom
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(geo_build_iterator_next_atom)

    call atom_list_next(this%aitr, atom, ierr)

    POP_SUB(geo_build_iterator_next_atom)
  end subroutine geo_build_iterator_next_atom

  ! ---------------------------------------------------------
  subroutine geo_build_iterator_next_species(this, spec, ierr)
    type(geo_build_iterator_t), intent(inout) :: this
    type(species_t),           pointer        :: spec
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(geo_build_iterator_next_species)

    call species_dict_next(this%sitr, spec, ierr)

    POP_SUB(geo_build_iterator_next_species)
  end subroutine geo_build_iterator_next_species

  ! ---------------------------------------------------------
  subroutine geo_build_iterator_copy(this, that)
    type(geo_build_iterator_t), intent(inout) :: this
    type(geo_build_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(geo_build_iterator_copy)

    this%self => that%self
    call atom_list_copy(this%aitr, that%aitr)
    call species_dict_copy(this%sitr, that%sitr)

    POP_SUB(geo_build_iterator_copy)
  end subroutine geo_build_iterator_copy

  ! ---------------------------------------------------------
  subroutine geo_build_iterator_end(this)
    type(geo_build_iterator_t), intent(inout) :: this

    PUSH_SUB(geo_build_iterator_end)

    nullify(this%self)
    call species_dict_end(this%sitr)
    call atom_list_end(this%aitr)
    call species_dict_end(this%sitr)

    POP_SUB(geo_build_iterator_end)
  end subroutine geo_build_iterator_end

end module geo_build_oct_m

!! Local Variables:
!! mode: f90
!! End:

