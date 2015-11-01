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

#define HASH_TEMPLATE_NAME base_system
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_system

module base_system_m

  use base_geometry_m
  use base_states_m
  use config_dict_m
  use geometry_m
  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m
  use simulation_m
  use space_m

#define LIST_TEMPLATE_NAME base_system
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_system
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                  &
    BASE_SYSTEM_OK,          &
    BASE_SYSTEM_KEY_ERROR,   &
    BASE_SYSTEM_EMPTY_ERROR

  public ::        &
    base_system_t

  public ::                &
    base_system__init__,   &
    base_system__start__,  &
    base_system__update__, &
    base_system__stop__,   &
    base_system__reset__,  &
    base_system__acc__,    &
    base_system__copy__,   &
    base_system__end__

  public ::             &
    base_system_new,    &
    base_system_del,    &
    base_system_init,   &
    base_system_start,  &
    base_system_update, &
    base_system_stop,   &
    base_system_sets,   &
    base_system_gets,   &
    base_system_set,    &
    base_system_get,    &
    base_system_copy,   &
    base_system_end

#define LIST_TEMPLATE_NAME base_system
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: BASE_SYSTEM_OK          = BASE_SYSTEM_HASH_OK
  integer, parameter :: BASE_SYSTEM_KEY_ERROR   = BASE_SYSTEM_HASH_KEY_ERROR
  integer, parameter :: BASE_SYSTEM_EMPTY_ERROR = BASE_SYSTEM_HASH_EMPTY_ERROR

  type :: base_system_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(base_system_t), pointer :: prnt   =>null()
    type(space_t)                :: space
    type(base_geometry_t)        :: geom
    type(base_states_t)          :: st
    type(config_dict_t)          :: dict
    type(base_system_hash_t)     :: hash
    type(base_system_list_t)     :: list
  end type base_system_t

  interface base_system__init__
    module procedure base_system__init__begin
    module procedure base_system__init__finish
    module procedure base_system__init__copy
  end interface base_system__init__

  interface base_system__copy__
    module procedure base_system__copy__begin
    module procedure base_system__copy__finish
  end interface base_system__copy__

  interface base_system_init
    module procedure base_system_init_system
    module procedure base_system_init_copy
  end interface base_system_init

  interface base_system_set
    module procedure base_system_set_simulation
  end interface base_system_set

  interface base_system_gets
    module procedure base_system_gets_config
    module procedure base_system_gets_name
  end interface base_system_gets

  interface base_system_get
    module procedure base_system_get_info
    module procedure base_system_get_config
    module procedure base_system_get_simulation
    module procedure base_system_get_space
    module procedure base_system_get_geom
    module procedure base_system_get_geometry
    module procedure base_system_get_states
  end interface base_system_get

  interface base_system_copy
    module procedure base_system_copy_system
  end interface base_system_copy

  interface base_system_end
    module procedure base_system_end_system
  end interface base_system_end

#define TEMPLATE_PREFIX base_system
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_system
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_system_new(this, that)
    type(base_system_t),  target, intent(inout) :: this
    type(base_system_t), pointer                :: that

    PUSH_SUB(base_system_new)

    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt => this
    call base_system_list_push(this%list, that)

    POP_SUB(base_system_new)
  end subroutine base_system_new

  ! ---------------------------------------------------------
  subroutine base_system__idel__(this)
    type(base_system_t), pointer :: this

    PUSH_SUB(base_system__idel__)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(base_system__idel__)
  end subroutine base_system__idel__

  ! ---------------------------------------------------------
  subroutine base_system_del(this)
    type(base_system_t), pointer :: this

    PUSH_SUB(base_system_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_system_list_del(this%prnt%list, this)
        call base_system_end(this)
        call base_system__idel__(this)
      end if
    end if

    POP_SUB(base_system_del)
  end subroutine base_system_del

  ! ---------------------------------------------------------
  subroutine base_system__inull__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__inull__)

    nullify(this%config, this%sim, this%prnt)

    POP_SUB(base_system__inull__)
  end subroutine base_system__inull__

  ! ---------------------------------------------------------
  subroutine base_system__iinit__(this, config)
    type(base_system_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_system__iinit__)

    nullify(cnfg)
    call base_system__inull__(this)
    this%config => config
    call json_get(this%config, "space", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call space_init(this%space, cnfg)
    nullify(cnfg)
    call config_dict_init(this%dict)
    call base_system_hash_init(this%hash)
    call base_system_list_init(this%list)

    POP_SUB(base_system__iinit__)
  end subroutine base_system__iinit__

  ! ---------------------------------------------------------
  subroutine base_system__init__begin(this, config)
    type(base_system_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_system__init__begin)

    nullify(cnfg)
    call base_system__iinit__(this, config)
    call json_get(this%config, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_geometry__init__(this%geom, this%space, cnfg)
    nullify(cnfg)
    call json_get(this%config, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_states__init__(this%st, cnfg)
    nullify(cnfg)

    POP_SUB(base_system__init__begin)
  end subroutine base_system__init__begin

  ! ---------------------------------------------------------
  subroutine base_system__init__finish(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__init__finish)

    call base_geometry__init__(this%geom)

    POP_SUB(base_system__init__finish)
  end subroutine base_system__init__finish

  ! ---------------------------------------------------------
  subroutine base_system__init__copy(this, that)
    type(base_system_t), intent(out) :: this
    type(base_system_t), intent(in)  :: that

    PUSH_SUB(base_system__init__copy)

    ASSERT(associated(that%config))
    call base_system__iinit__(this, that%config)
    call base_geometry__init__(this%geom, that%geom)
    call base_states__init__(this%st, that%st)

    POP_SUB(base_system__init__copy)
  end subroutine base_system__init__copy

  ! ---------------------------------------------------------
  subroutine base_system_init_system(this, config)
    type(base_system_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_system_init_system)

    call base_system__init__begin(this, config)
    call base_system__init__finish(this)

    POP_SUB(base_system_init_system)
  end subroutine base_system_init_system

  ! ---------------------------------------------------------
  recursive subroutine base_system_init_copy(this, that)
    type(base_system_t), intent(out) :: this
    type(base_system_t), intent(in)  :: that

    type(base_system_iterator_t) :: iter
    type(base_system_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_system_init_copy)

    nullify(cnfg, osub, isub)
    call base_system__init__(this, that)
    call base_system_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_system_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_SYSTEM_OK)exit
      call base_system_new(this, osub)
      call base_system_init(osub, isub)
      call base_system_sets(this, osub, cnfg)
    end do
    call base_system_end(iter)
    call base_system__init__(this)
    nullify(cnfg, osub, isub)

    POP_SUB(base_system_init_copy)
  end subroutine base_system_init_copy

  ! ---------------------------------------------------------
  subroutine base_system__istart__(this, sim)
    type(base_system_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_system__istart__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim

    POP_SUB(base_system__istart__)
  end subroutine base_system__istart__

  ! ---------------------------------------------------------
  subroutine base_system__start__(this, sim)
    type(base_system_t),          intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    PUSH_SUB(base_system__start__)

    if(present(sim))then
      call base_system__istart__(this, sim)
    else
      if(.not.associated(this%sim))then
        ASSERT(associated(this%prnt))
        ASSERT(associated(this%prnt%sim))
        call base_system__istart__(this, this%prnt%sim)
      end if
    end if
    call base_states__start__(this%st, sim)

    POP_SUB(base_system__start__)
  end subroutine base_system__start__

  ! ---------------------------------------------------------
  recursive subroutine base_system_start(this, sim)
    type(base_system_t),          intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    type(base_system_iterator_t) :: iter
    type(base_system_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_system_start)

    nullify(subs)
    call base_system_init(iter, this)
    do
      nullify(subs)
      call base_system_next(iter, subs, ierr)
      if(ierr/=BASE_SYSTEM_OK)exit
      call base_system_start(subs, sim)
    end do
    call base_system_end(iter)
    nullify(subs)
    call base_system__start__(this, sim)

    POP_SUB(base_system_start)
  end subroutine base_system_start

  ! ---------------------------------------------------------
  subroutine base_system__update__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_states__update__(this%st)

    POP_SUB(base_system__update__)
  end subroutine base_system__update__

  ! ---------------------------------------------------------
  recursive subroutine base_system_update(this)
    type(base_system_t), intent(inout) :: this

    type(base_system_iterator_t) :: iter
    type(base_system_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_system_update)

    nullify(subs)
    call base_system_init(iter, this)
    do
      nullify(subs)
      call base_system_next(iter, subs, ierr)
      if(ierr/=BASE_SYSTEM_OK)exit
      call base_system_update(subs)
    end do
    call base_system_end(iter)
    nullify(subs)
    call base_system__update__(this)

    POP_SUB(base_system_update)
  end subroutine base_system_update

  ! ---------------------------------------------------------
  subroutine base_system__stop__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_states__stop__(this%st)

    POP_SUB(base_system__stop__)
  end subroutine base_system__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_system_stop(this)
    type(base_system_t), intent(inout) :: this

    type(base_system_iterator_t) :: iter
    type(base_system_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_system_stop)

    nullify(subs)
    call base_system_init(iter, this)
    do
      nullify(subs)
      call base_system_next(iter, subs, ierr)
      if(ierr/=BASE_SYSTEM_OK)exit
      call base_system_stop(subs)
    end do
    call base_system_end(iter)
    nullify(subs)
    call base_system__stop__(this)

    POP_SUB(base_system_stop)
  end subroutine base_system_stop

  ! ---------------------------------------------------------
  subroutine base_system__reset__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_states__reset__(this%st)

    POP_SUB(base_system__reset__)
  end subroutine base_system__reset__

  ! ---------------------------------------------------------
  subroutine base_system__acc__(this, that)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_states__acc__(this%st, that%st)

    POP_SUB(base_system__acc__)
  end subroutine base_system__acc__

  ! ---------------------------------------------------------
  subroutine base_system__sets__(this, that, config)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    PUSH_SUB(base_system__sets__)

    ASSERT(this%space==that%space)
    call base_geometry_sets(this%geom, that%geom, config)
    call base_states_sets(this%st, that%st, config)

    POP_SUB(base_system__sets__)
  end subroutine base_system__sets__

  ! ---------------------------------------------------------
  subroutine base_system_sets(this, that, config)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_system_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_system_hash_set(this%hash, config, that)
    call base_system__sets__(this, that, config)

    POP_SUB(base_system_sets)
  end subroutine base_system_sets

  ! ---------------------------------------------------------
  subroutine base_system_gets_config(this, config, that)
    type(base_system_t),  intent(in) :: this
    type(json_object_t),  intent(in) :: config
    type(base_system_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_system_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_system_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_SYSTEM_OK) nullify(that)

    POP_SUB(base_system_gets_config)
  end subroutine base_system_gets_config

  ! ---------------------------------------------------------
  subroutine base_system_gets_name(this, name, that)
    type(base_system_t),  intent(in) :: this
    character(len=*),     intent(in) :: name
    type(base_system_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_system_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_system_gets(this, config, that)

    POP_SUB(base_system_gets_name)
  end subroutine base_system_gets_name

  ! ---------------------------------------------------------
  subroutine base_system_set_simulation(this, that)
    type(base_system_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that

    PUSH_SUB(base_system_set_simulation)

    ASSERT(.not.associated(this%sim))
    this%sim => that

    POP_SUB(base_system_set_simulation)
  end subroutine base_system_set_simulation

  ! ---------------------------------------------------------
  subroutine base_system_get_info(this, charge, nspin)
    type(base_system_t),     intent(in)  :: this
    real(kind=wp), optional, intent(out) :: charge
    integer,       optional, intent(out) :: nspin

    PUSH_SUB(base_system_get_info)

    call base_states_get(this%st, charge=charge, nspin=nspin)

    POP_SUB(base_system_get_info)
  end subroutine base_system_get_info
 
  ! ---------------------------------------------------------
  subroutine base_system_get_config(this, that)
    type(base_system_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_system_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_system_get_config)
  end subroutine base_system_get_config

  ! ---------------------------------------------------------
  subroutine base_system_get_simulation(this, that)
    type(base_system_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_system_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_system_get_simulation)
  end subroutine base_system_get_simulation

  ! ---------------------------------------------------------
  subroutine base_system_get_space(this, that)
    type(base_system_t), target, intent(in) :: this
    type(space_t),      pointer             :: that

    PUSH_SUB(base_system_get_space)

    that => this%space

    POP_SUB(base_system_get_space)
  end subroutine base_system_get_space

  ! ---------------------------------------------------------
  subroutine base_system_get_geom(this, that)
    type(base_system_t),     target, intent(in) :: this
    type(base_geometry_t),  pointer             :: that

    PUSH_SUB(base_system_get_geom)

    that => this%geom

    POP_SUB(base_system_get_geom)
  end subroutine base_system_get_geom

  ! ---------------------------------------------------------
  subroutine base_system_get_geometry(this, that)
    type(base_system_t), intent(in) :: this
    type(geometry_t),   pointer     :: that

    PUSH_SUB(base_system_get_geometry)

    call base_geometry_get(this%geom, that)

    POP_SUB(base_system_get_geometry)
  end subroutine base_system_get_geometry

  ! ---------------------------------------------------------
  subroutine base_system_get_states(this, that)
    type(base_system_t),  target, intent(in) :: this
    type(base_states_t), pointer             :: that

    PUSH_SUB(base_system_get_states)

    that => this%st

    POP_SUB(base_system_get_states)
  end subroutine base_system_get_states

  ! ---------------------------------------------------------
  subroutine base_system__icopy__(this, that)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system__icopy__)

    call base_system__iend__(this)
    if(associated(that%config))then
      call base_system__iinit__(this, that%config)
      if(associated(that%sim)) call base_system__istart__(this, that%sim)
    end if

    POP_SUB(base_system__icopy__)
  end subroutine base_system__icopy__

  ! ---------------------------------------------------------
  subroutine base_system__copy__begin(this, that)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system__copy__begin)

    call base_system__icopy__(this, that)
    call base_geometry__copy__(this%geom, that%geom)
    call base_states__copy__(this%st, that%st)

    POP_SUB(base_system__copy__begin)
  end subroutine base_system__copy__begin

  ! ---------------------------------------------------------
  subroutine base_system__copy__finish(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__copy__finish)

    call base_geometry__copy__(this%geom)

    POP_SUB(base_system__copy__finish)
  end subroutine base_system__copy__finish

  ! ---------------------------------------------------------
  recursive subroutine base_system_copy_system(this, that)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that

    type(base_system_iterator_t) :: iter
    type(base_system_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_system_copy_system)

    nullify(cnfg, osub, isub)
    call base_system_end(this)
    call base_system__copy__(this, that)
    call base_system_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_system_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_SYSTEM_OK)exit
      call base_system_new(this, osub)
      call base_system_copy(osub, isub)
      call base_system_sets(this, osub, cnfg)
    end do
    call base_system_end(iter)
    call base_system__copy__(this)
    nullify(cnfg, osub, isub)

    POP_SUB(base_system_copy_system)
  end subroutine base_system_copy_system

  ! ---------------------------------------------------------
  subroutine base_system__iend__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__iend__)

    call base_system__inull__(this)
    call space_end(this%space)
    call config_dict_end(this%dict)
    call base_system_hash_end(this%hash)
    call base_system_list_end(this%list)

    POP_SUB(base_system__iend__)
  end subroutine base_system__iend__

  ! ---------------------------------------------------------
  subroutine base_system__end__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__end__)

    call base_system__iend__(this)
    call base_geometry__end__(this%geom)
    call base_states__end__(this%st)

    POP_SUB(base_system__end__)
  end subroutine base_system__end__

  ! ---------------------------------------------------------
  recursive subroutine base_system_end_system(this)
    type(base_system_t), intent(inout) :: this

    type(base_system_t), pointer :: subs

    PUSH_SUB(base_system_end_system)

    do
      nullify(subs)
      call base_system_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_system_end(subs)
      call base_system__idel__(subs)
    end do
    nullify(subs)
    call base_system__end__(this)

    POP_SUB(base_system_end_system)
  end subroutine base_system_end_system

#define TEMPLATE_PREFIX base_system
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_system_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
