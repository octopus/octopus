#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

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

#define HASH_TEMPLATE_NAME base_geometry
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_geometry

module base_geometry_oct_m

  use atom_oct_m
  use config_dict_oct_m
  use geo_build_oct_m
  use geo_intrf_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use space_oct_m
  use species_oct_m

#define LIST_TEMPLATE_NAME base_geometry
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_geometry
#define EXTENDED_TYPE
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef EXTENDED_TYPE
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::          &
    base_geometry_t

  public ::                &
    base_geometry__init__, &
    base_geometry__copy__, &
    base_geometry__end__

  public ::             &
    base_geometry_new,  &
    base_geometry_del,  &
    base_geometry_init, &
    base_geometry_sets, &
    base_geometry_gets, &
    base_geometry_get,  &
    base_geometry_copy, &
    base_geometry_end

  public ::                    &
    BASE_GEOMETRY_OK,          &
    BASE_GEOMETRY_EMPTY_ERROR

#define LIST_TEMPLATE_NAME base_geometry
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: BASE_GEOMETRY_OK          = CONFIG_DICT_OK
  integer, parameter :: BASE_GEOMETRY_EMPTY_ERROR = CONFIG_DICT_EMPTY_ERROR

  type :: base_geometry_t
    private
    type(json_object_t),   pointer :: config =>null()
    type(space_t),         pointer :: space  =>null()
    type(base_geometry_t), pointer :: prnt   =>null()
    type(geo_intrf_t)              :: igeo
    type(config_dict_t)            :: dict
    type(base_geometry_hash_t)     :: hash
    type(base_geometry_list_t)     :: list
  end type base_geometry_t

  type :: base_geometry_iterator_t
    private
    type(base_geometry_t), pointer :: self =>null()
    type(geo_intrf_iterator_t)     :: gitr
    type(config_dict_iterator_t)   :: iter
  end type base_geometry_iterator_t

  interface base_geometry__init__
    module procedure base_geometry__init__begin
    module procedure base_geometry__init__copy
    module procedure base_geometry__init__finish
  end interface base_geometry__init__

  interface base_geometry__copy__
    module procedure base_geometry__copy__begin
    module procedure base_geometry__copy__finish
  end interface base_geometry__copy__

  interface base_geometry_init
    module procedure base_geometry_init_type
    module procedure base_geometry_init_copy
  end interface base_geometry_init

  interface base_geometry_next
    module procedure base_geometry_iterator_next_atom
    module procedure base_geometry_iterator_next_species
    module procedure base_geometry_iterator_next_geometry
  end interface base_geometry_next

  interface base_geometry_gets
    module procedure base_geometry_gets_config
    module procedure base_geometry_gets_name
  end interface base_geometry_gets

  interface base_geometry_get
    module procedure base_geometry_get_config
    module procedure base_geometry_get_space
    module procedure base_geometry_get_geometry
    module procedure base_geometry_get_geo_intrf
  end interface base_geometry_get

  interface base_geometry_copy
    module procedure base_geometry_copy_type
  end interface base_geometry_copy

  interface base_geometry_end
    module procedure base_geometry_end_type
  end interface base_geometry_end

#define TEMPLATE_PREFIX base_geometry
#define EXTENDED_TYPE
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef EXTENDED_TYPE
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_geometry
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_geometry__new__(this)
    type(base_geometry_t), pointer :: this

    PUSH_SUB(base_geometry__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_geometry__new__)
  end subroutine base_geometry__new__

  ! ---------------------------------------------------------
  subroutine base_geometry__del__(this)
    type(base_geometry_t), pointer :: this

    PUSH_SUB(base_geometry__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_geometry__del__)
  end subroutine base_geometry__del__

  ! ---------------------------------------------------------
  subroutine base_geometry_new(this, that)
    type(base_geometry_t),  target, intent(inout) :: this
    type(base_geometry_t), pointer                :: that

    PUSH_SUB(base_geometry_new)

    nullify(that)
    call base_geometry__new__(that)
    that%prnt => this
    call base_geometry_list_push(this%list, that)

    POP_SUB(base_geometry_new)
  end subroutine base_geometry_new

  ! ---------------------------------------------------------
  subroutine base_geometry_del(this)
    type(base_geometry_t), pointer :: this

    PUSH_SUB(base_geometry_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_geometry_list_del(this%prnt%list, this)
        call base_geometry_end(this)
        call base_geometry__del__(this)
      end if
    end if

    POP_SUB(base_geometry_del)
  end subroutine base_geometry_del

  ! ---------------------------------------------------------
  subroutine base_geometry__init__begin(this, space, config)
    type(base_geometry_t),       intent(out) :: this
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_geometry__init__begin)

    this%config => config
    this%space => space
    call json_get(config, "molecule", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call geo_intrf_init(this%igeo, space, cnfg)
    call config_dict_init(this%dict)
    call base_geometry_hash_init(this%hash)
    call base_geometry_list_init(this%list)

    POP_SUB(base_geometry__init__begin)
  end subroutine base_geometry__init__begin

  ! ---------------------------------------------------------
  subroutine base_geometry__build__(this, that)
    type(base_geometry_t), intent(in)  :: this
    type(geometry_t),      intent(out) :: that

    type(base_geometry_iterator_t) :: iter
    type(base_geometry_t), pointer :: subs
    type(json_object_t),   pointer :: cnfg
    type(json_array_t),    pointer :: list
    type(geometry_t),      pointer :: pgeo
    type(geo_build_t)              :: bgeo
    integer                        :: ierr

    PUSH_SUB(base_geometry__build__)

    call geo_build_init(bgeo, this%space)
    call base_geometry_init(iter, this)
    do
      nullify(cnfg, list, subs, pgeo)
      call base_geometry_next(iter, cnfg, subs, ierr)
      if(ierr/=BASE_GEOMETRY_OK)exit
      call base_geometry_get(subs, pgeo)
      ASSERT(associated(pgeo))
      call json_get(cnfg, "positions", list, ierr)
      ASSERT(ierr==JSON_OK)
      if(json_len(list)>0)then
        call geo_build_extend(bgeo, pgeo, list)
      else
        call geo_build_extend(bgeo, pgeo)
      end if
    end do
    call base_geometry_end(iter)
    call geo_build_export(bgeo, that)
    call geo_build_end(bgeo)
    nullify(cnfg, list, subs, pgeo)

    POP_SUB(base_geometry__build__)
  end subroutine base_geometry__build__

  ! ---------------------------------------------------------
  subroutine base_geometry__init__finish(this)
    type(base_geometry_t), intent(inout) :: this

    type(geometry_t), pointer :: pgeo

    PUSH_SUB(base_geometry__init__finish)

    nullify(pgeo)
    if(.not.geo_intrf_assoc(this%igeo))then
      call geo_intrf_new(this%igeo, pgeo)
      call base_geometry__build__(this, pgeo)
      nullify(pgeo)
    end if

    POP_SUB(base_geometry__init__finish)
  end subroutine base_geometry__init__finish

  ! ---------------------------------------------------------
  subroutine base_geometry__init__copy(this, that)
    type(base_geometry_t), intent(out) :: this
    type(base_geometry_t), intent(in)  :: that

    PUSH_SUB(base_geometry__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%space))
    call base_geometry__init__(this, that%space, that%config)

    POP_SUB(base_geometry__init__copy)
  end subroutine base_geometry__init__copy

  ! ---------------------------------------------------------
  subroutine base_geometry_init_type(this, space, config)
    type(base_geometry_t), intent(out) :: this
    type(space_t),         intent(in)  :: space
    type(json_object_t),   intent(in)  :: config

    PUSH_SUB(base_geometry_init_type)

    call base_geometry__init__(this, space, config)
    call base_geometry__init__(this)

    POP_SUB(base_geometry_init_type)
  end subroutine base_geometry_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_geometry_init_copy(this, that)
    type(base_geometry_t), intent(out) :: this
    type(base_geometry_t), intent(in)  :: that

    type(base_geometry_iterator_t) :: iter
    type(base_geometry_t), pointer :: osub, isub
    type(json_object_t),   pointer :: cnfg
    integer                        :: ierr

    PUSH_SUB(base_geometry_init_copy)

    nullify(cnfg, osub, isub)
    call base_geometry__init__(this, that)
    call base_geometry_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_geometry_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_GEOMETRY_OK)exit
      call base_geometry_new(this, osub)
      call base_geometry_init(osub, isub)
      call base_geometry_sets(this, osub, cnfg)
    end do
    call base_geometry_end(iter)
    call base_geometry__init__(this)
    nullify(cnfg, osub, isub)

    POP_SUB(base_geometry_init_copy)
  end subroutine base_geometry_init_copy

  ! ---------------------------------------------------------
  subroutine base_geometry_sets(this, that, config)
    type(base_geometry_t), intent(inout) :: this
    type(base_geometry_t), intent(in)    :: that
    type(json_object_t),   intent(in)    :: config

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_geometry_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_geometry_hash_set(this%hash, config, that)

    POP_SUB(base_geometry_sets)
  end subroutine base_geometry_sets

  ! ---------------------------------------------------------
  subroutine base_geometry_gets_config(this, config, that)
    type(base_geometry_t),  intent(in) :: this
    type(json_object_t),    intent(in) :: config
    type(base_geometry_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_geometry_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_geometry_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_GEOMETRY_OK) nullify(that)

    POP_SUB(base_geometry_gets_config)
  end subroutine base_geometry_gets_config

  ! ---------------------------------------------------------
  subroutine base_geometry_gets_name(this, name, that)
    type(base_geometry_t),  intent(in) :: this
    character(len=*),       intent(in) :: name
    type(base_geometry_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_geometry_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_geometry_gets(this, config, that)

    POP_SUB(base_geometry_gets_name)
  end subroutine base_geometry_gets_name

  ! ---------------------------------------------------------
  subroutine base_geometry_get_config(this, that)
    type(base_geometry_t), target, intent(in) :: this
    type(json_object_t),  pointer             :: that

    PUSH_SUB(base_geometry_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_geometry_get_config)
  end subroutine base_geometry_get_config

  ! ---------------------------------------------------------
  subroutine base_geometry_get_space(this, that)
    type(base_geometry_t), target, intent(in) :: this
    type(space_t),        pointer             :: that

    PUSH_SUB(base_geometry_get_space)

    nullify(that)
    if(associated(this%space)) that => this%space

    POP_SUB(base_geometry_get_space)
  end subroutine base_geometry_get_space

  ! ---------------------------------------------------------
  subroutine base_geometry_get_geometry(this, that)
    type(base_geometry_t), intent(in) :: this
    type(geometry_t),     pointer     :: that

    PUSH_SUB(base_geometry_get_geometry)

    call geo_intrf_get(this%igeo, that)

    POP_SUB(base_geometry_get_geometry)
  end subroutine base_geometry_get_geometry

  ! ---------------------------------------------------------
  subroutine base_geometry_get_geo_intrf(this, that)
    type(base_geometry_t), target, intent(in) :: this
    type(geo_intrf_t),    pointer             :: that

    PUSH_SUB(base_geometry_get_geo_intrf)

    nullify(that)
    that => this%igeo

    POP_SUB(base_geometry_get_geo_intrf)
  end subroutine base_geometry_get_geo_intrf

  ! ---------------------------------------------------------
  subroutine base_geometry__copy__begin(this, that)
    type(base_geometry_t), intent(inout) :: this
    type(base_geometry_t), intent(in)    :: that

    PUSH_SUB(base_geometry__copy__begin)

    call base_geometry__end__(this)
    if(associated(that%config).and.associated(that%space))&
      call base_geometry__init__(this, that)

    POP_SUB(base_geometry__copy__begin)
  end subroutine base_geometry__copy__begin

  ! ---------------------------------------------------------
  subroutine base_geometry__copy__finish(this)
    type(base_geometry_t), intent(inout) :: this

    PUSH_SUB(base_geometry__copy__finish)

    if(associated(this%config).and.associated(this%space))&
      call base_geometry__init__(this)

    POP_SUB(base_geometry__copy__finish)
  end subroutine base_geometry__copy__finish

  ! ---------------------------------------------------------
  recursive subroutine base_geometry_copy_type(this, that)
    type(base_geometry_t), intent(inout) :: this
    type(base_geometry_t), intent(in)    :: that

    type(base_geometry_iterator_t)   :: iter
    type(base_geometry_t),   pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_geometry_copy_type)

    nullify(cnfg, osub, isub)
    call base_geometry_end(this)
    call base_geometry__copy__(this, that)
    call base_geometry_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_geometry_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_GEOMETRY_OK)exit
      call base_geometry_new(this, osub)
      call base_geometry_copy(osub, isub)
      call base_geometry_sets(this, osub, cnfg)
    end do
    call base_geometry_end(iter)
    call base_geometry__copy__(this)
    nullify(cnfg, osub, isub)

    POP_SUB(base_geometry_copy_type)
  end subroutine base_geometry_copy_type

  ! ---------------------------------------------------------
  subroutine base_geometry__end__(this)
    type(base_geometry_t), target, intent(inout) :: this

    PUSH_SUB(base_geometry__end__)

    nullify(this%config, this%prnt)
    call geo_intrf_end(this%igeo)
    call config_dict_end(this%dict)
    call base_geometry_hash_end(this%hash)
    call base_geometry_list_end(this%list)

    POP_SUB(base_geometry__end__)
  end subroutine base_geometry__end__

  ! ---------------------------------------------------------
  recursive subroutine base_geometry_end_type(this)
    type(base_geometry_t), intent(inout) :: this

    type(base_geometry_t), pointer :: subs

    PUSH_SUB(base_geometry_end_type)

    do
      nullify(subs)
      call base_geometry_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_geometry_end(subs)
      call base_geometry__del__(subs)
    end do
    nullify(subs)
    call base_geometry__end__(this)

    POP_SUB(base_geometry_end_type)
  end subroutine base_geometry_end_type

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator__init__(this, that)
    type(base_geometry_iterator_t), intent(out) :: this
    type(base_geometry_t),          intent(in)  :: that

    PUSH_SUB(base_geometry_iterator__init__)

    call geo_intrf_init(this%gitr, that%igeo)

    POP_SUB(base_geometry_iterator__init__)
  end subroutine base_geometry_iterator__init__

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator_next_atom(this, atom, ierr)
    type(base_geometry_iterator_t), intent(inout) :: this
    type(atom_t),              pointer        :: atom
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(base_geometry_iterator_next_atom)

    call geo_intrf_next(this%gitr, atom, ierr)

    POP_SUB(base_geometry_iterator_next_atom)
  end subroutine base_geometry_iterator_next_atom

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator_next_species(this, spec, ierr)
    type(base_geometry_iterator_t), intent(inout) :: this
    type(species_t),           pointer        :: spec
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(base_geometry_iterator_next_species)

    call geo_intrf_next(this%gitr, spec, ierr)

    POP_SUB(base_geometry_iterator_next_species)
  end subroutine base_geometry_iterator_next_species

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator_next_geometry(this, geo, ierr)
    type(base_geometry_iterator_t), intent(inout) :: this
    type(geometry_t),          pointer        :: geo
    integer,          optional, intent(out)   :: ierr

    type(base_geometry_t), pointer :: geometry
    integer                    :: jerr

    PUSH_SUB(base_geometry_iterator_next_geometry)

    nullify(geo, geometry)
    call base_geometry_next(this, geometry, jerr)
    if(jerr==BASE_GEOMETRY_OK) call base_geometry_get(geometry, geo)
    if(present(ierr)) ierr = jerr
    nullify(geometry)

    POP_SUB(base_geometry_iterator_next_geometry)
  end subroutine base_geometry_iterator_next_geometry

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator__copy__(this, that)
    type(base_geometry_iterator_t), intent(out) :: this
    type(base_geometry_iterator_t), intent(in)  :: that

    PUSH_SUB(base_geometry_iterator__copy__)

    call geo_intrf_copy(this%gitr, that%gitr)

    POP_SUB(base_geometry_iterator__copy__)
  end subroutine base_geometry_iterator__copy__

  ! ---------------------------------------------------------
  subroutine base_geometry_iterator__end__(this)
    type(base_geometry_iterator_t), intent(inout) :: this

    PUSH_SUB(base_geometry_iterator__end__)

    call geo_intrf_end(this%gitr)

    POP_SUB(base_geometry_iterator__end__)
  end subroutine base_geometry_iterator__end__

#define TEMPLATE_PREFIX base_geometry
#define EXTENDED_TYPE
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef EXTENDED_TYPE
#undef TEMPLATE_PREFIX

end module base_geometry_oct_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

