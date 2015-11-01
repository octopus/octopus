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

#define HASH_TEMPLATE_NAME base_model
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_model

module base_model_m

  use base_geometry_m
  use base_hamiltonian_m
  use base_system_m
  use config_dict_m
  use geometry_m
  use global_m
  use grid_m
  use json_m
  use messages_m
  use profiling_m
  use simulation_m
  use space_m

#define LIST_TEMPLATE_NAME base_model
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_model
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                 &
    BASE_MODEL_OK,          &
    BASE_MODEL_KEY_ERROR,   &
    BASE_MODEL_EMPTY_ERROR

  public ::       &
    base_model_t

  public ::               &
    base_model__init__,   &
    base_model__start__,  &
    base_model__update__, &
    base_model__stop__,   &
    base_model__reset__,  &
    base_model__acc__,    &
    base_model__copy__,   &
    base_model__end__

  public ::            &
    base_model_new,    &
    base_model_del,    &
    base_model_init,   &
    base_model_start,  &
    base_model_update, &
    base_model_stop,   &
    base_model_sets,   &
    base_model_gets,   &
    base_model_get,    &
    base_model_copy,   &
    base_model_end

#define LIST_TEMPLATE_NAME base_model
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: BASE_MODEL_OK          = BASE_MODEL_HASH_OK
  integer, parameter :: BASE_MODEL_KEY_ERROR   = BASE_MODEL_HASH_KEY_ERROR
  integer, parameter :: BASE_MODEL_EMPTY_ERROR = BASE_MODEL_HASH_EMPTY_ERROR

  type :: base_model_t
    private
    type(json_object_t), pointer :: config =>null()
    type(base_model_t),  pointer :: prnt   =>null()
    type(simulation_t)           :: sim
    type(base_system_t)          :: sys
    type(base_hamiltonian_t)     :: hm
    type(config_dict_t)          :: dict
    type(base_model_hash_t)      :: hash
    type(base_model_list_t)      :: list
  end type base_model_t

  interface base_model__init__
    module procedure base_model__init__begin
    module procedure base_model__init__finish
    module procedure base_model__init__copy
  end interface base_model__init__

  interface base_model__copy__
    module procedure base_model__copy__begin
    module procedure base_model__copy__finish
  end interface base_model__copy__

  interface base_model_init
    module procedure base_model_init_type
    module procedure base_model_init_copy
  end interface base_model_init

  interface base_model_gets
    module procedure base_model_gets_config
    module procedure base_model_gets_name
  end interface base_model_gets

  interface base_model_get
    module procedure base_model_get_config
    module procedure base_model_get_simulation
    module procedure base_model_get_system
    module procedure base_model_get_hamiltonian
  end interface base_model_get

  interface base_model_copy
    module procedure base_model_copy_type
  end interface base_model_copy

  interface base_model_end
    module procedure base_model_end_type
  end interface base_model_end

#define TEMPLATE_PREFIX base_model
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_model
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  subroutine base_model_new(this, that)
    type(base_model_t),  target, intent(inout) :: this
    type(base_model_t), pointer                :: that

    PUSH_SUB(base_model_new)

    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt => this
    call base_model_list_push(this%list, that)

    POP_SUB(base_model_new)
  end subroutine base_model_new

  ! ---------------------------------------------------------
  subroutine base_model__idel__(this)
    type(base_model_t), pointer :: this

    PUSH_SUB(base_model__idel__)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(base_model__idel__)
  end subroutine base_model__idel__

  ! ---------------------------------------------------------
  subroutine base_model_del(this)
    type(base_model_t), pointer :: this

    PUSH_SUB(base_model_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_model_list_del(this%prnt%list, this)
        call base_model_end(this)
        call base_model__idel__(this)
      end if
    end if

    POP_SUB(base_model_del)
  end subroutine base_model_del

  ! ---------------------------------------------------------
  subroutine base_model__inull__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__inull__)

    nullify(this%config, this%prnt)

    POP_SUB(base_model__inull__)
  end subroutine base_model__inull__

  ! ---------------------------------------------------------
  subroutine base_model__iinit__(this, config)
    type(base_model_t),          intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(base_model__iinit__)

    call base_model__inull__(this)
    this%config => config
    call config_dict_init(this%dict)
    call base_model_hash_init(this%hash)
    call base_model_list_init(this%list)

    POP_SUB(base_model__iinit__)
  end subroutine base_model__iinit__

  ! ---------------------------------------------------------
  subroutine base_model__init__begin(this, config)
    type(base_model_t),  intent(out) :: this
    type(json_object_t), intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_model__init__begin)

    nullify(cnfg)
    call base_model__iinit__(this, config)
    call json_get(this%config, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_system__init__(this%sys, cnfg)
    nullify(cnfg)

    POP_SUB(base_model__init__begin)
  end subroutine base_model__init__begin

  ! ---------------------------------------------------------
  subroutine base_model__init__finish(this)
    type(base_model_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(space_t),       pointer :: space
    type(geometry_t),    pointer :: geo
    integer                      :: ierr

    PUSH_SUB(base_model__init__finish)

    nullify(cnfg, space, geo)
    call base_system__init__(this%sys)
    call base_system_get(this%sys, space)
    ASSERT(associated(space))
    call base_system_get(this%sys, geo)
    ASSERT(associated(geo))
    call json_get(this%config, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call simulation_init(this%sim, geo, space, cnfg)
    nullify(cnfg, space, geo)
    call json_get(this%config, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_hamiltonian__init__(this%hm, this%sys, cnfg)
    call base_model__build__(this)
    nullify(cnfg)

    POP_SUB(base_model__init__finish)
  end subroutine base_model__init__finish

  ! ---------------------------------------------------------
  subroutine base_model__init__copy(this, that)
    type(base_model_t), intent(out) :: this
    type(base_model_t), intent(in)  :: that

    PUSH_SUB(base_model__init__copy)

    ASSERT(associated(that%config))
    call base_model__iinit__(this, that%config)
    call base_system__init__(this%sys, that%sys)

    POP_SUB(base_model__init__copy)
  end subroutine base_model__init__copy

  ! ---------------------------------------------------------
  subroutine base_model_init_type(this, config)
    type(base_model_t),  intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_model_init_type)

    call base_model__init__(this, config)
    call base_model__init__(this)

    POP_SUB(base_model_init_type)
  end subroutine base_model_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_model_init_copy(this, that)
    type(base_model_t), intent(out) :: this
    type(base_model_t), intent(in)  :: that

    type(base_model_iterator_t)  :: iter
    type(base_model_t),  pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_model_init_copy)

    call base_model__init__(this, that)
    call base_model_init(iter, that)
    do
      nullify(osub, isub, cnfg)
      call base_model_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_new(this, osub)
      call base_model_init(osub, isub)
      call base_model_sets(this, osub, cnfg)
    end do
    call base_model_end(iter)
    nullify(osub, isub, cnfg)
    call base_model__init__(this)

    POP_SUB(base_model_init_copy)
  end subroutine base_model_init_copy

  ! ---------------------------------------------------------
  subroutine base_model__build__(this)
    type(base_model_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(base_model_iterator_t)  :: iter
    type(base_model_t),  pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_model__build__)

    call base_model_init(iter, this)
    do
      nullify(cnfg, subs)
      call base_model_next(iter, cnfg, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call simulation_extend(this%sim, subs%sim, cnfg)
      call base_hamiltonian_sets(this%hm, subs%hm, cnfg)
    end do
    call base_model_end(iter)
    nullify(cnfg, subs)

    POP_SUB(base_model__build__)
  end subroutine base_model__build__

  ! ---------------------------------------------------------
  subroutine base_model__start__(this, grid)
    type(base_model_t),     intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid

    PUSH_SUB(base_model__start__)

    ASSERT(associated(this%config))
    call simulation_start(this%sim, grid)
    call base_system__start__(this%sys, this%sim)
    call base_hamiltonian__start__(this%hm, this%sim)

    POP_SUB(base_model__start__)
  end subroutine base_model__start__

  ! ---------------------------------------------------------
  recursive subroutine base_model_start(this, grid)
    type(base_model_t),     intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid

    type(base_model_iterator_t) :: iter
    type(base_model_t), pointer :: subs
    integer                     :: ierr

    PUSH_SUB(base_model_start)

    call base_model_init(iter, this)
    do
      nullify(subs)
      call base_model_next(iter, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_start(subs, grid)
    end do
    call base_model_end(iter)
    nullify(subs)
    call base_model__start__(this, grid)

    POP_SUB(base_model_start)
  end subroutine base_model_start

  ! ---------------------------------------------------------
  subroutine base_model__update__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__update__)

    call base_system__update__(this%sys)
    call base_hamiltonian__update__(this%hm)

    POP_SUB(base_model__update__)
  end subroutine base_model__update__

  ! ---------------------------------------------------------
  recursive subroutine base_model_update(this)
    type(base_model_t), intent(inout) :: this

    type(base_model_iterator_t) :: iter
    type(base_model_t), pointer :: subs
    integer                     :: ierr

    PUSH_SUB(base_model_update)

    call base_model_init(iter, this)
    do
      nullify(subs)
      call base_model_next(iter, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_update(subs)
    end do
    call base_model_end(iter)
    nullify(subs)
    call base_model__update__(this)

    POP_SUB(base_model_update)
  end subroutine base_model_update

  ! ---------------------------------------------------------
  subroutine base_model__stop__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__stop__)

    call base_system__stop__(this%sys)
    call base_hamiltonian__stop__(this%hm)

    POP_SUB(base_model__stop__)
  end subroutine base_model__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_model_stop(this)
    type(base_model_t), intent(inout) :: this

    type(base_model_iterator_t) :: iter
    type(base_model_t), pointer :: subs
    integer                     :: ierr

    PUSH_SUB(base_model_stop)

    call base_model_init(iter, this)
    do
      nullify(subs)
      call base_model_next(iter, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_stop(subs)
    end do
    call base_model_end(iter)
    nullify(subs)
    call base_model__stop__(this)

    POP_SUB(base_model_stop)
  end subroutine base_model_stop

  ! ---------------------------------------------------------
  subroutine base_model__reset__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__reset__)

    call base_system__reset__(this%sys)
    call base_hamiltonian__reset__(this%hm)

    POP_SUB(base_model__reset__)
  end subroutine base_model__reset__

  ! ---------------------------------------------------------
  subroutine base_model__acc__(this, that)
    type(base_model_t), intent(inout) :: this
    type(base_model_t), intent(in)    :: that

    PUSH_SUB(base_model__acc__)

    call base_system__acc__(this%sys, that%sys)
    call base_hamiltonian__acc__(this%hm, that%hm)

    POP_SUB(base_model__acc__)
  end subroutine base_model__acc__

  ! ---------------------------------------------------------
  subroutine base_model__sets__(this, that, config)
    type(base_model_t),  intent(inout) :: this
    type(base_model_t),  intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    PUSH_SUB(base_model__sets__)

    call base_system_sets(this%sys, that%sys, config)

    POP_SUB(base_model__sets__)
  end subroutine base_model__sets__

  ! ---------------------------------------------------------
  subroutine base_model_sets(this, that, config)
    type(base_model_t),  intent(inout) :: this
    type(base_model_t),  intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_model_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_model_hash_set(this%hash, config, that)
    call base_model__sets__(this, that, config)

    POP_SUB(base_model_sets)
  end subroutine base_model_sets

  ! ---------------------------------------------------------
  subroutine base_model_gets_config(this, config, that)
    type(base_model_t),  intent(in) :: this
    type(json_object_t), intent(in) :: config
    type(base_model_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_model_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_model_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_MODEL_OK) nullify(that)

    POP_SUB(base_model_gets_config)
  end subroutine base_model_gets_config

  ! ---------------------------------------------------------
  subroutine base_model_gets_name(this, name, that)
    type(base_model_t),  intent(in) :: this
    character(len=*),    intent(in) :: name
    type(base_model_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_model_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_model_gets(this, config, that)

    POP_SUB(base_model_gets_name)
  end subroutine base_model_gets_name

  ! ---------------------------------------------------------
  subroutine base_model_get_config(this, that)
    type(base_model_t),   target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_model_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_model_get_config)
  end subroutine base_model_get_config

  ! ---------------------------------------------------------
  subroutine base_model_get_simulation(this, that)
    type(base_model_t),  target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_model_get_simulation)

    that => this%sim

    POP_SUB(base_model_get_simulation)
  end subroutine base_model_get_simulation

  ! ---------------------------------------------------------
  subroutine base_model_get_system(this, that)
    type(base_model_t),   target, intent(in) :: this
    type(base_system_t), pointer             :: that

    PUSH_SUB(base_model_get_system)

    that => this%sys

    POP_SUB(base_model_get_system)
  end subroutine base_model_get_system

  ! ---------------------------------------------------------
  subroutine base_model_get_hamiltonian(this, that)
    type(base_model_t),        target, intent(in) :: this
    type(base_hamiltonian_t), pointer             :: that

    PUSH_SUB(base_model_get_hamiltonian)

    that => this%hm

    POP_SUB(base_model_get_hamiltonian)
  end subroutine base_model_get_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_model__icopy__(this, that)
    type(base_model_t), intent(inout) :: this
    type(base_model_t), intent(in)    :: that

    PUSH_SUB(base_model__icopy__)

    call base_model__iend__(this)
    if(associated(that%config)) call base_model__iinit__(this, that%config)

    POP_SUB(base_model__icopy__)
  end subroutine base_model__icopy__

  ! ---------------------------------------------------------
  subroutine base_model__copy__begin(this, that)
    type(base_model_t), intent(inout) :: this
    type(base_model_t), intent(in)    :: that

    PUSH_SUB(base_model__copy__begin)

    call base_model__icopy__(this, that)
    call simulation_copy(this%sim, that%sim)
    call base_system__copy__(this%sys, that%sys)
    call base_hamiltonian__copy__(this%hm, that%hm)

    POP_SUB(base_model__copy__begin)
  end subroutine base_model__copy__begin

  ! ---------------------------------------------------------
  subroutine base_model__copy__finish(this)
    type(base_model_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(base_model_iterator_t)  :: iter
    type(base_model_t),  pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_model__copy__finish)

    nullify(cnfg, subs)
    call base_system__copy__(this%sys)
    call base_model_init(iter, this)
    do
      nullify(cnfg, subs)
      call base_model_next(iter, cnfg, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_hamiltonian_sets(this%hm, subs%hm, cnfg)
    end do
    call base_model_end(iter)
    nullify(cnfg, subs)

    POP_SUB(base_model__copy__finish)
  end subroutine base_model__copy__finish

  ! ---------------------------------------------------------
  recursive subroutine base_model_copy_type(this, that)
    type(base_model_t), intent(inout) :: this
    type(base_model_t), intent(in)    :: that

    type(base_model_iterator_t)  :: iter
    type(base_model_t),  pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_model_copy_type)

    nullify(cnfg, osub, isub)
    call base_model_end(this)
    call base_model__copy__(this, that)
    call base_model_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_model_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_new(this, osub)
      call base_model_copy(osub, isub)
      call base_model_sets(this, osub, cnfg)
    end do
    call base_model_end(iter)
    call base_model__copy__(this)
    nullify(cnfg, osub, isub)

    POP_SUB(base_model_copy_type)
  end subroutine base_model_copy_type

  ! ---------------------------------------------------------
  subroutine base_model__iend__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__iend__)

    call base_model__inull__(this)
    call config_dict_end(this%dict)
    call base_model_hash_end(this%hash)
    call base_model_list_end(this%list)

    POP_SUB(base_model__iend__)
  end subroutine base_model__iend__

  ! ---------------------------------------------------------
  subroutine base_model__end__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__end__)

    call base_model__iend__(this)
    call base_hamiltonian__end__(this%hm)
    call base_system__end__(this%sys)
    call simulation_end(this%sim)

    POP_SUB(base_model__end__)
  end subroutine base_model__end__

  ! ---------------------------------------------------------
  recursive subroutine base_model_end_type(this)
    type(base_model_t), intent(inout) :: this

    type(base_model_t), pointer :: subs

    PUSH_SUB(base_model_end_type)

    do
      nullify(subs)
      call base_model_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_model_end(subs)
      call base_model__idel__(subs)
    end do
    nullify(subs)
    call base_model__end__(this)

    POP_SUB(base_model_end_type)
  end subroutine base_model_end_type

#define TEMPLATE_PREFIX base_model
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_model_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
