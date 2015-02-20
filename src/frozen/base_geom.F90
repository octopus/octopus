#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME
#undef LIST_INCLUDE_PREFIX
#undef LIST_INCLUDE_HEADER
#undef LIST_INCLUDE_BODY

#define LIST_TEMPLATE_NAME atom
#include "tlist.F90"
#undef LIST_TEMPLATE_NAME

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#define DICT_TEMPLATE_NAME species
#include "tdict.F90"
#undef DICT_TEMPLATE_NAME

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_KEY_TYPE_MODULE_NAME
#undef HASH_KEY_FUNCTION_NAME
#undef HASH_KEY_FUNCTION_MODULE_NAME
#undef HASH_VAL_TEMPLATE_NAME
#undef HASH_VAL_TYPE_NAME
#undef HASH_VAL_TYPE_MODULE_NAME
#undef HASH_INCLUDE_PREFIX
#undef HASH_INCLUDE_HEADER
#undef HASH_INCLUDE_BODY

#define HASH_TEMPLATE_NAME base_geom
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_geom

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_geom_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_object_t, json_hash
  use kinds_m,  only: wp

  use atom_m,        only: atom_t, atom_set_species, atom_get_label, atom_end
  use basis_m,       only: basis_t, basis_init, basis_to_external, basis_end
  use distributed_m, only: distributed_nullify
  use json_m,        only: JSON_OK
  use json_m,        only: json_object_t, json_init, json_get, json_end
  use json_m,        only: json_array_t, json_array_iterator_t, json_len, json_next
  use space_m,       only: operator(==), space_t
  use species_m

  use config_dict_m, only: &
    CONFIG_DICT_OK,        &
    CONFIG_DICT_NAME_LEN

  use config_dict_m, only: &
    config_dict_t,         &
    config_dict_init,      &
    config_dict_set,       &
    config_dict_get,       &
    config_dict_copy,      &
    config_dict_end

  use geometry_m, only:             &
    geometry_t,                     &
    geometry_init_from_data_object, &
    geometry_copy,                  &
    geometry_end

  use atom_list_m,  only: &
    ATOM_LIST_OK,         &
    atom_list_t,          &
    atom_list_len,        &
    atom_list_init,       &
    atom_list_next,       &
    atom_list_append,     &
    atom_list_push,       &
    atom_list_pop,        &
    atom_list_copy,       &
    atom_list_end

  use atom_list_m,  only: &
    atom_list_iterator_t

  use species_dict_m, only: &
    SPECIES_DICT_OK,        &
    species_dict_t,         &
    species_dict_len,       &
    species_dict_init,      &
    species_dict_next,      &
    species_dict_pop,       &
    species_dict_set,       &
    species_dict_get,       &
    species_dict_copy,      &
    species_dict_end

  use species_dict_m, only:  &
    species_dict_iterator_t

  implicit none

  private
  public ::            &
    base_geom__init__, &
    base_geom__add__

  public ::         &
    base_geom_init, &
    base_geom_next, &
    base_geom_get,  &
    base_geom_copy, &
    base_geom_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_geom_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(geometry_t),    pointer :: pgeo
    type(geometry_t)             :: igeo
    type(atom_list_t)            :: list
    type(species_dict_t)         :: sdct
    type(config_dict_t)          :: cdct
    type(base_geom_hash_t)       :: hash
  end type base_geom_t

  type, public :: base_geom_iterator_t
    private
    type(atom_list_iterator_t)      :: aitr
    type(species_dict_iterator_t)   :: sitr
    type(base_geom_hash_iterator_t) :: hitr
  end type base_geom_iterator_t

  interface base_geom__init__
    module procedure base_geom__init__begin
    module procedure base_geom__init__finish
  end interface base_geom__init__

  interface base_geom_init
    module procedure base_geom_init_geom
    module procedure base_geom_iterator_init
  end interface base_geom_init

  interface base_geom_next
    module procedure base_geom_iterator_next_atom
    module procedure base_geom_iterator_next_species
    module procedure base_geom_iterator_next_config_geom
    module procedure base_geom_iterator_next_config
    module procedure base_geom_iterator_next_geom
    module procedure base_geom_iterator_next_geometry
  end interface base_geom_next

  interface base_geom_get
    module procedure base_geom_get_config
    module procedure base_geom_get_space
    module procedure base_geom_get_geometry
  end interface base_geom_get

  interface base_geom_copy
    module procedure base_geom_copy_geom
    module procedure base_geom_iterator_copy
  end interface base_geom_copy

  interface base_geom_end
    module procedure base_geom_end_geom
    module procedure base_geom_iterator_end
  end interface base_geom_end

  integer, public, parameter :: BASE_GEOM_OK          = BASE_GEOM_HASH_OK
  integer, public, parameter :: BASE_GEOM_KEY_ERROR   = BASE_GEOM_HASH_KEY_ERROR
  integer, public, parameter :: BASE_GEOM_EMPTY_ERROR = BASE_GEOM_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_geom_build_geom(this, geo, basis)
    type(base_geom_t),       intent(inout) :: this
    type(geometry_t),        intent(in)    :: geo
    type(basis_t), optional, intent(in)    :: basis
    !
    type(atom_t),    pointer :: atom
    type(species_t), pointer :: spec
    integer                  :: indx
    !
    PUSH_SUB(base_geom_build_geom)
    ASSERT(this%space==geo%space)
    do indx=1,geo%natoms
      SAFE_ALLOCATE(atom)
      atom=geo%atom(indx)
      if(present(basis))&
        call basis_to_external(basis, geo%atom(indx)%x, atom%x)
      call atom_list_push(this%list, atom)
    end do
    do indx=1,geo%nspecies
      SAFE_ALLOCATE(spec)
      spec=geo%species(indx)
      call species_dict_set(this%sdct, species_label(spec), spec)
    end do
    POP_SUB(base_geom_build_geom)
    return
  end subroutine base_geom_build_geom

  ! ---------------------------------------------------------
  subroutine base_geom_build_geometry(this, geo, space)
    type(base_geom_t),     intent(in)  :: this
    type(geometry_t),      intent(out) :: geo
    type(space_t), target, intent(in)  :: space
    !
    type(base_geom_iterator_t) :: iter
    type(atom_t),      pointer :: atom
    type(species_t),   pointer :: spec
    integer                    :: indx, ierr
    !
    PUSH_SUB(base_geom_build_geometry)
    !call geometry_end(geo)
    geo%natoms=atom_list_len(this%list)
    ASSERT(geo%natoms>0)
    SAFE_ALLOCATE(geo%atom(geo%natoms))
    geo%space=>space
    geo%nspecies=species_dict_len(this%sdct)
    ASSERT(geo%nspecies>0)
    SAFE_ALLOCATE(geo%species(geo%nspecies))
    call base_geom_init(iter, this)
    do indx=1, geo%nspecies
      nullify(spec)
      call base_geom_next(iter, spec, ierr)
      ASSERT(ierr==0)
      call species_set_index(spec, indx)
      geo%species(indx)=spec
    end do
    do indx=1, geo%natoms
      nullify(atom, spec)
      call base_geom_next(iter, atom, ierr)
      ASSERT(ierr==0)
      call species_dict_get(this%sdct, atom_get_label(atom), spec, ierr)
      ASSERT(ierr==SPECIES_DICT_OK)
      geo%atom(indx)=atom
      call atom_set_species(geo%atom(indx), geo%species(species_index(spec)))
    end do
    call base_geom_end(iter)
    nullify(atom, spec)
    geo%ncatoms=0
    geo%catom=>null()
    geo%only_user_def=.false.
    geo%species_time_dependent=.false.
    geo%kinetic_energy=M_ZERO
    geo%nlpp=.false.
    geo%nlcc=.false.
    call distributed_nullify(geo%atoms_dist, geo%natoms)
    nullify(geo%ionic_interaction_type)
    nullify(geo%ionic_interaction_parameter)
    geo%reduced_coordinates=.false.
    POP_SUB(base_geom_build_geometry)
    return
  end subroutine base_geom_build_geometry

  ! ---------------------------------------------------------
  subroutine base_geom__init__begin(this, space, config)
    type(base_geom_t),           intent(out) :: this
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(base_geom__init__begin)
    this%config=>config
    this%space=>space
    nullify(this%pgeo)
    call geometry_init_from_data_object(this%igeo, this%space, config)
    call atom_list_init(this%list)
    call species_dict_init(this%sdct, this%igeo%nspecies)
    call config_dict_init(this%cdct)
    call base_geom_hash_init(this%hash)
    call base_geom_build_geom(this, this%igeo)
    POP_SUB(base_geom__init__begin)
    return
  end subroutine base_geom__init__begin

  ! ---------------------------------------------------------
  subroutine base_geom__init__finish(this)
    type(base_geom_t), target, intent(inout) :: this
    !
    PUSH_SUB(base_geom__init__finish)
    if(base_geom_hash_len(this%hash)>0)then
      if(associated(this%pgeo))then
        if(.not.associated(this%pgeo,this%igeo))then
          call geometry_end(this%pgeo)
          SAFE_DEALLOCATE_P(this%pgeo)
        end if
      end if
      nullify(this%pgeo)
      SAFE_ALLOCATE(this%pgeo)
      call base_geom_build_geometry(this, this%pgeo, this%space)
    else
      this%pgeo=>this%igeo
    end if
    POP_SUB(base_geom__init__finish)
    return
  end subroutine base_geom__init__finish

  ! ---------------------------------------------------------
  subroutine base_geom_init_geom(this, space, config)
    type(base_geom_t),   intent(out) :: this
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config
    !
    PUSH_SUB(base_geom_init_geom)
    call base_geom__init__begin(this, space, config)
    call base_geom__init__finish(this)
    POP_SUB(base_geom_init_geom)
    return
  end subroutine base_geom_init_geom

  ! ---------------------------------------------------------
  subroutine base_geom__add__(this, that, config)
    type(base_geom_t),   intent(inout) :: this
    type(base_geom_t),   intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    type(geometry_t),           pointer :: geo
    type(json_object_t),        pointer :: cnfg
    type(json_array_t),         pointer :: list
    character(len=CONFIG_DICT_NAME_LEN) :: name
    type(json_array_iterator_t)         :: iter
    type(basis_t)                       :: base
    integer                             :: ierr
    !
    PUSH_SUB(base_geom__add__)
    nullify(geo, cnfg, list)
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%cdct, trim(adjustl(name)), config)
    call base_geom_hash_set(this%hash, config, that)
    call base_geom_get(that, geo)
    ASSERT(associated(geo))
    call json_get(config, "positions", list, ierr)
    ASSERT(ierr==JSON_OK)
    if(json_len(list)>0)then
      call json_init(iter, list)
      do
        nullify(cnfg)
        call json_next(iter, cnfg, ierr)
        if(ierr/=JSON_OK)exit
        call basis_init(base, this%space, cnfg)
        call base_geom_build_geom(this, geo, base)
        call basis_end(base)
      end do
      call json_end(iter)
      nullify(cnfg)
    else
      call base_geom_build_geom(this, geo)
    end if
    POP_SUB(base_geom__add__)
    return
  end subroutine base_geom__add__

  ! ---------------------------------------------------------
  subroutine base_geom__get__(this, name, that)
    type(base_geom_t),  intent(inout) :: this
    character(len=*),   intent(in)    :: name
    type(base_geom_t), pointer        :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_geom__get__)
    nullify(that)
    call config_dict_get(this%cdct, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_geom_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_GEOM_OK)nullify(that)
    end if
    POP_SUB(base_geom__get__)
    return
  end subroutine base_geom__get__

  ! ---------------------------------------------------------
  subroutine base_geom_get_config(this, that)
    type(base_geom_t),    target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(base_geom_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_geom_get_config)
    return
  end subroutine base_geom_get_config

  ! ---------------------------------------------------------
  subroutine base_geom_get_space(this, that)
    type(base_geom_t), target, intent(in) :: this
    type(space_t),    pointer             :: that
    !
    PUSH_SUB(base_geom_get_space)
    nullify(that)
    if(associated(this%space))&
      that=>this%space
    POP_SUB(base_geom_get_space)
    return
  end subroutine base_geom_get_space

  ! ---------------------------------------------------------
  subroutine base_geom_get_geometry(this, that)
    type(base_geom_t), target, intent(in) :: this
    type(geometry_t), pointer             :: that
    !
    PUSH_SUB(base_geom_get_geometry)
    nullify(that)
    if(associated(this%pgeo))&
      that=>this%pgeo
    POP_SUB(base_geom_get_geometry)
    return
  end subroutine base_geom_get_geometry

  ! ---------------------------------------------------------
  subroutine base_geom__copy__list(this, that)
    type(atom_list_t), intent(inout) :: this
    type(atom_list_t), intent(in)    :: that
    !
    type(atom_list_iterator_t) :: iter
    type(atom_t),      pointer :: oatm, iatm
    integer                    :: ierr
    !
    PUSH_SUB(base_geom__copy__list)
    call atom_list_init(iter, that)
    do
      nullify(oatm, iatm)
      call atom_list_next(iter, iatm, ierr)
      if(ierr/=ATOM_LIST_OK)exit
      SAFE_ALLOCATE(oatm)
      oatm=iatm
      call atom_list_append(this, oatm)
    end do
    call atom_list_end(iter)
    POP_SUB(base_geom__copy__list)
    return
  end subroutine base_geom__copy__list

  ! ---------------------------------------------------------
  subroutine base_geom__copy__dict(this, that)
    type(species_dict_t), intent(inout) :: this
    type(species_dict_t), intent(in)    :: that
    !
    type(species_dict_iterator_t) :: iter
    type(species_t),      pointer :: ospc, ispc
    character(len=LABEL_LEN)      :: name
    integer                       :: ierr
    !
    PUSH_SUB(base_geom__copy__dict)
    call species_dict_init(iter, that)
    do
      nullify(ospc, ispc)
      call species_dict_next(iter, name, ispc, ierr)
      if(ierr/=SPECIES_DICT_OK)exit
      SAFE_ALLOCATE(ospc)
      ospc=ispc
      call species_dict_set(this, trim(adjustl(name)), ospc)
    end do
    call species_dict_end(iter)
    nullify(ospc, ispc)
    POP_SUB(base_geom__copy__dict)
    return
  end subroutine base_geom__copy__dict

  ! ---------------------------------------------------------
  subroutine base_geom_copy_geom(this, that)
    type(base_geom_t), target, intent(inout) :: this
    type(base_geom_t), target, intent(in)    :: that
    !
    PUSH_SUB(base_geom_copy_geom)
    call base_geom_end_geom(this)
    this%config=>that%config
    this%space=>that%space
    call geometry_copy(this%igeo, that%igeo)
    if(associated(that%pgeo))then
      if(associated(that%pgeo, that%igeo))then
        this%pgeo=>this%igeo
      else
        SAFE_ALLOCATE(this%pgeo)
        call geometry_copy(this%pgeo, that%pgeo)
      end if
    end if
    call base_geom__copy__list(this%list, that%list)
    call base_geom__copy__dict(this%sdct, that%sdct)
    call config_dict_copy(this%cdct, that%cdct)
    call base_geom_hash_copy(this%hash, that%hash)
    POP_SUB(base_geom_copy_geom)
    return
  end subroutine base_geom_copy_geom

  ! ---------------------------------------------------------
  subroutine base_geom__end__list(this)
    type(atom_list_t), intent(inout) :: this
    !
    type(atom_t), pointer :: atom
    !
    PUSH_SUB(base_geom__end__list)
    nullify(atom)
    do
      call atom_list_pop(this, atom)
      if(.not.associated(atom))exit
      call atom_end(atom)
      SAFE_DEALLOCATE_P(atom)
      nullify(atom)
    end do
    call atom_list_end(this)
    POP_SUB(base_geom__end__list)
    return
  end subroutine base_geom__end__list

  ! ---------------------------------------------------------
  subroutine base_geom__end__dict(this)
    type(species_dict_t), intent(inout) :: this
    !
    type(species_t), pointer :: spec
    !
    PUSH_SUB(base_geom__end__dict)
    nullify(spec)
    do
      call species_dict_pop(this, spec)
      if(.not.associated(spec))exit
      call species_end(spec)
      SAFE_DEALLOCATE_P(spec)
      nullify(spec)
    end do
    call species_dict_end(this)
    POP_SUB(base_geom__end__dict)
    return
  end subroutine base_geom__end__dict

  ! ---------------------------------------------------------
  subroutine base_geom_end_geom(this)
    type(base_geom_t), target, intent(inout) :: this
    !
    PUSH_SUB(base_geom_end_geom)
    nullify(this%config, this%space)
    if(associated(this%pgeo))then
      if(.not.associated(this%pgeo, this%igeo))then
        call geometry_end(this%pgeo)
        SAFE_DEALLOCATE_P(this%pgeo)
      end if
    end if
    call geometry_end(this%igeo)
    nullify(this%pgeo)
    call base_geom__end__list(this%list)
    call base_geom__end__dict(this%sdct)
    call config_dict_end(this%cdct)
    call base_geom_hash_end(this%hash)
    POP_SUB(base_geom_end_geom)
    return
  end subroutine base_geom_end_geom

  ! ! ---------------------------------------------------------
  ! recursive subroutine base_geom_end_geom(this)
  !   type(base_geom_t), intent(inout) :: this
  !   !
  !   type(base_geom_iterator_t) :: iter
  !   type(base_geom_t), pointer :: subs
  !   integer                    :: ierr
  !   !
  !   P!USH_SUB(base_geom_end_geom)
  !   nullify(subs)
  !   call base_geom_init(iter, this)
  !   do
  !     nullify(subs)
  !     call base_geom_next(iter, subs, ierr)
  !     if(ierr/=BASE_GEOM_OK)exit
  !     call base_geom_end(subs)
  !   end do
  !   call base_geom_end(iter)
  !   nullify(subs)
  !   call base_geom__end__(this)
  !   P!OP_SUB(base_geom_end_geom)
  !   return
  ! end subroutine base_geom_end_geom

  ! ---------------------------------------------------------
  subroutine base_geom_iterator_init(this, that)
    type(base_geom_iterator_t), intent(out) :: this
    type(base_geom_t),          intent(in)  :: that
    !
    PUSH_SUB(base_geom_iterator_init)
    call atom_list_init(this%aitr, that%list)
    call species_dict_init(this%sitr, that%sdct)
    call base_geom_hash_init(this%hitr, that%hash)
    POP_SUB(base_geom_iterator_init)
    return
  end subroutine base_geom_iterator_init

  ! ---------------------------------------------------------
  subroutine base_geom_iterator_next_atom(this, atom, ierr)
    type(base_geom_iterator_t), intent(inout) :: this
    type(atom_t),              pointer        :: atom
    integer,          optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_geom_iterator_next_atom)
    call atom_list_next(this%aitr, atom, ierr)
    POP_SUB(base_geom_iterator_next_atom)
    return
  end subroutine base_geom_iterator_next_atom

  ! ---------------------------------------------------------
  subroutine base_geom_iterator_next_species(this, spec, ierr)
    type(base_geom_iterator_t), intent(inout) :: this
    type(species_t),           pointer        :: spec
    integer,          optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_geom_iterator_next_species)
    call species_dict_next(this%sitr, spec, ierr)
    POP_SUB(base_geom_iterator_next_species)
    return
  end subroutine base_geom_iterator_next_species

  ! ---------------------------------------------------------
  subroutine base_geom_iterator_next_config_geom(this, config, geom, ierr)
    type(base_geom_iterator_t), intent(inout) :: this
    type(json_object_t),       pointer        :: config
    type(base_geom_t),         pointer        :: geom
    integer,          optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_geom_iterator_next_config_geom)
    call base_geom_hash_next(this%hitr, config, geom, ierr)
    POP_SUB(base_geom_iterator_next_config_geom)
    return
  end subroutine base_geom_iterator_next_config_geom

  ! ---------------------------------------------------------
  subroutine base_geom_iterator_next_config(this, config, ierr)
    type(base_geom_iterator_t), intent(inout) :: this
    type(json_object_t),       pointer        :: config
    integer,          optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_geom_iterator_next_config)
    call base_geom_hash_next(this%hitr, config, ierr)
    POP_SUB(base_geom_iterator_next_config)
    return
  end subroutine base_geom_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_geom_iterator_next_geom(this, geom, ierr)
    type(base_geom_iterator_t), intent(inout) :: this
    type(base_geom_t),         pointer        :: geom
    integer,          optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_geom_iterator_next_geom)
    call base_geom_hash_next(this%hitr, geom, ierr)
    POP_SUB(base_geom_iterator_next_geom)
    return
  end subroutine base_geom_iterator_next_geom

  ! ---------------------------------------------------------
  subroutine base_geom_iterator_next_geometry(this, geo, ierr)
    type(base_geom_iterator_t), intent(inout) :: this
    type(geometry_t),          pointer        :: geo
    integer,          optional, intent(out)   :: ierr
    !
    type(base_geom_t), pointer :: geom
    !
    PUSH_SUB(base_geom_iterator_next_geometry)
    nullify(geo, geom)
    call base_geom_hash_next(this%hitr, geom, ierr)
    call base_geom_get(geom, geo)
    nullify(geom)
    POP_SUB(base_geom_iterator_next_geometry)
    return
  end subroutine base_geom_iterator_next_geometry

  ! ---------------------------------------------------------
  subroutine base_geom_iterator_copy(this, that)
    type(base_geom_iterator_t), intent(out) :: this
    type(base_geom_iterator_t), intent(in)  :: that
    !
    PUSH_SUB(base_geom_iterator_copy)
    call atom_list_copy(this%aitr, that%aitr)
    call species_dict_copy(this%sitr, that%sitr)
    call base_geom_hash_copy(this%hitr, that%hitr)
    POP_SUB(base_geom_iterator_copy)
    return
  end subroutine base_geom_iterator_copy

  ! ---------------------------------------------------------
  elemental subroutine base_geom_iterator_end(this)
    type(base_geom_iterator_t), intent(inout) :: this
    !
    call atom_list_end(this%aitr)
    call species_dict_end(this%sitr)
    call base_geom_hash_end(this%hitr)
    return
  end subroutine base_geom_iterator_end

end module base_geom_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

