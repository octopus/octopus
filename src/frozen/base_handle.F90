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

#define HASH_TEMPLATE_NAME base_handle
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_handle

module base_handle_oct_m

  use base_model_oct_m
  use config_dict_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m

#define LIST_TEMPLATE_NAME base_handle
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_handle
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::         &
    HNDL_TYPE_NONE

  public ::                  &
    BASE_HANDLE_OK,          &
    BASE_HANDLE_KEY_ERROR,   &
    BASE_HANDLE_EMPTY_ERROR

  public ::        &
    base_handle_t

  public ::                &
    base_handle__init__,   &
    base_handle__start__,  &
    base_handle__update__, &
    base_handle__stop__,   &
    base_handle__reset__,  &
    base_handle__copy__,   &
    base_handle__end__

  public ::             &
    base_handle_new,    &
    base_handle_del,    &
    base_handle_init,   &
    base_handle_start,  &
    base_handle_update, &
    base_handle_stop,   &
    base_handle_sets,   &
    base_handle_gets,   &
    base_handle_get,    &
    base_handle_copy,   &
    base_handle_end

#define LIST_TEMPLATE_NAME base_handle
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: HNDL_TYPE_NONE = 0

  integer, parameter :: BASE_HANDLE_OK          = BASE_HANDLE_HASH_OK
  integer, parameter :: BASE_HANDLE_KEY_ERROR   = BASE_HANDLE_HASH_KEY_ERROR
  integer, parameter :: BASE_HANDLE_EMPTY_ERROR = BASE_HANDLE_HASH_EMPTY_ERROR

  type :: base_handle_t
    private
    type(json_object_t), pointer :: config =>null()
    type(base_handle_t), pointer :: prnt   =>null()
    integer                      :: type   = HNDL_TYPE_NONE
    type(simulation_t)           :: sim
    type(base_model_t)           :: model
    type(config_dict_t)          :: dict
    type(base_handle_hash_t)     :: hash
    type(base_handle_list_t)     :: list
  end type base_handle_t

  interface base_handle__init__
    module procedure base_handle__init__begin
    module procedure base_handle__init__finish
    module procedure base_handle__init__copy
  end interface base_handle__init__

  interface base_handle__copy__
    module procedure base_handle__copy__begin
    module procedure base_handle__copy__finish
  end interface base_handle__copy__

  interface base_handle_init
    module procedure base_handle_init_type
    module procedure base_handle_init_pass
    module procedure base_handle_init_copy
  end interface base_handle_init

  interface base_handle_gets
    module procedure base_handle_gets_config
    module procedure base_handle_gets_name
  end interface base_handle_gets

  interface base_handle_get
    module procedure base_handle_get_info
    module procedure base_handle_get_config
    module procedure base_handle_get_simulation
    module procedure base_handle_get_model
  end interface base_handle_get

  interface base_handle_copy
    module procedure base_handle_copy_type
  end interface base_handle_copy

  interface base_handle_end
    module procedure base_handle_end_type
  end interface base_handle_end

#define TEMPLATE_PREFIX base_handle
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_handle
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_handle__new__(this)
    type(base_handle_t), pointer :: this

    PUSH_SUB(base_handle__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_handle__new__)
  end subroutine base_handle__new__

  ! ---------------------------------------------------------
  subroutine base_handle__del__(this)
    type(base_handle_t), pointer :: this

    PUSH_SUB(base_handle__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_handle__del__)
  end subroutine base_handle__del__

  ! ---------------------------------------------------------
  subroutine base_handle_new(this, that)
    type(base_handle_t),  intent(inout) :: this
    type(base_handle_t), pointer        :: that

    PUSH_SUB(base_handle_new)

    nullify(that)
    call base_handle__new__(that)
    call base_handle__set__(that, this)
    call base_handle_list_push(this%list, that)

    POP_SUB(base_handle_new)
  end subroutine base_handle_new

  ! ---------------------------------------------------------
  subroutine base_handle_del(this)
    type(base_handle_t), pointer :: this

    PUSH_SUB(base_handle_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_handle_list_del(this%prnt%list, this)
        call base_handle_end(this)
        call base_handle__del__(this)
      end if
    end if

    POP_SUB(base_handle_del)
  end subroutine base_handle_del

  ! ---------------------------------------------------------
  subroutine base_handle__init__begin(this, config)
    type(base_handle_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_handle__init__begin)

    nullify(cnfg)
    this%config => config
    call json_get(this%config, "type", this%type, ierr)
    if(ierr/=JSON_OK) this%type = HNDL_TYPE_NONE
    call json_get(this%config, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_model__init__(this%model, cnfg)
    nullify(cnfg)
    call config_dict_init(this%dict)
    call base_handle_hash_init(this%hash)
    call base_handle_list_init(this%list)
    
    POP_SUB(base_handle__init__begin)
  end subroutine base_handle__init__begin

  ! ---------------------------------------------------------
  subroutine base_handle__init__build(this, list)
    type(base_handle_t), intent(inout) :: this
    type(json_array_t),  intent(in)    :: list

    type(json_array_iterator_t)  :: iter
    type(json_object_t), pointer :: cnfg
    type(base_handle_t), pointer :: hndl
    integer                      :: ierr

    PUSH_SUB(base_handle__init__build)

    call json_init(iter, list)
    do
      nullify(cnfg, hndl)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call base_handle_new(this, hndl)
      call base_handle_init(hndl, cnfg)
      call base_handle__set__(hndl, this)
      call base_handle_sets(this, hndl, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, hndl)

    POP_SUB(base_handle__init__build)
  end subroutine base_handle__init__build

  ! ---------------------------------------------------------
  subroutine base_handle__init__finish(this)
    type(base_handle_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(space_t),       pointer :: space
    type(geometry_t),    pointer :: geo
    integer                      :: ierr

    PUSH_SUB(base_handle__init__finish)

    nullify(cnfg, space, geo)
    call base_model__init__(this%model)
    call base_model_get(this%model, space)
    ASSERT(associated(space))
    call base_model_get(this%model, geo)
    ASSERT(associated(geo))
    call json_get(this%config, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call simulation_init(this%sim, geo, space, cnfg)
    nullify(cnfg, space, geo)

    POP_SUB(base_handle__init__finish)
  end subroutine base_handle__init__finish

  ! ---------------------------------------------------------
  subroutine base_handle__init__copy(this, that)
    type(base_handle_t), intent(out) :: this
    type(base_handle_t), intent(in)  :: that

    PUSH_SUB(base_handle__init__copy)

    ASSERT(associated(that%config))
    call base_handle__init__(this, that%config)
    if(simulation_assoc(that%sim))then
      call simulation_copy(this%sim, that%sim)
      call base_model__start__(this%model, this%sim)
    end if

    POP_SUB(base_handle__init__copy)
  end subroutine base_handle__init__copy

  ! ---------------------------------------------------------
  recursive subroutine base_handle_init_type(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_handle_init_type)

    call base_handle_init_pass(this, config, base_handle_init_type)

    POP_SUB(base_handle_init_type)
  end subroutine base_handle_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_handle_init_pass(this, config, handle_init)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    interface
      subroutine handle_init(this, config)
        use json_oct_m
        import :: base_handle_t
        type(base_handle_t), intent(out) :: this
        type(json_object_t), intent(in)  :: config
      end subroutine handle_init
    end interface

    type(json_array_iterator_t)   :: iter
    type(json_object_t),  pointer :: cnfg
    type(json_array_t),   pointer :: list
    type(base_handle_t),  pointer :: hndl
    integer                       :: ierr

    PUSH_SUB(base_handle_init_pass)

    nullify(cnfg, list, hndl)
    call base_handle__init__(this, config)
    call json_get(config, "systems", list, ierr)
    if(ierr==JSON_OK)then
      call json_init(iter, list)
      do
        nullify(cnfg, hndl)
        call json_next(iter, cnfg, ierr)
        if(ierr/=JSON_OK)exit
        call base_handle_new(this, hndl)
        call handle_init(hndl, cnfg)
        call base_handle__set__(hndl, this)
        call base_handle_sets(this, hndl, cnfg)
      end do
      call json_end(iter)
      nullify(cnfg, hndl)
    end if
    call base_handle__init__(this)
    nullify(list)

    POP_SUB(base_handle_init_pass)
  end subroutine base_handle_init_pass

  ! ---------------------------------------------------------
  recursive subroutine base_handle_init_copy(this, that)
    type(base_handle_t), intent(out) :: this
    type(base_handle_t), intent(in)  :: that

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_handle_init_copy)

    call base_handle__init__(this, that)
    call base_handle_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_handle_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_new(this, osub)
      call base_handle_init(osub, isub)
      call base_handle_sets(this, osub, cnfg)
    end do
    call base_handle_end(iter)
    call base_handle__init__(this)
    nullify(cnfg, osub, isub)

    POP_SUB(base_handle_init_copy)
  end subroutine base_handle_init_copy

  ! ---------------------------------------------------------
  subroutine base_handle__build__(this)
    type(base_handle_t), intent(inout) :: this

    type(base_handle_iterator_t) :: iter
    type(json_object_t), pointer :: cnfg
    type(base_handle_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_handle__build__)

    call base_handle_init(iter, this)
    do
      nullify(cnfg, subs)
      call base_handle_next(iter, cnfg, subs, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call simulation_extend(this%sim, subs%sim, cnfg)
    end do
    call base_handle_end(iter)
    nullify(cnfg, subs)

    POP_SUB(base_handle__build__)
  end subroutine base_handle__build__

  ! ---------------------------------------------------------
  subroutine base_handle__start__(this, grid)
    type(base_handle_t),    intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid

    PUSH_SUB(base_handle__start__)

    ASSERT(associated(this%config))
    call base_handle__build__(this)
    call simulation_start(this%sim, grid)
    ASSERT(simulation_assoc(this%sim))
    call base_model__start__(this%model, this%sim)

    POP_SUB(base_handle__start__)
  end subroutine base_handle__start__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_start(this, grid)
    type(base_handle_t),    intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_handle_start)

    call base_handle_init(iter, this)
    do
      nullify(subs)
      call base_handle_next(iter, subs, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_start(subs, grid)
    end do
    call base_handle_end(iter)
    nullify(subs)
    call base_handle__start__(this, grid)

    POP_SUB(base_handle_start)
  end subroutine base_handle_start

  ! ---------------------------------------------------------
  subroutine base_handle__update__(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__update__)

    ASSERT(associated(this%config))
    call base_model__update__(this%model)

    POP_SUB(base_handle__update__)
  end subroutine base_handle__update__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_update(this)
    type(base_handle_t), intent(inout) :: this

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_handle_update)

    call base_handle_init(iter, this)
    do
      nullify(subs)
      call base_handle_next(iter, subs, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_update(subs)
    end do
    call base_handle_end(iter)
    nullify(subs)
    call base_handle__update__(this)

    POP_SUB(base_handle_update)
  end subroutine base_handle_update

  ! ---------------------------------------------------------
  subroutine base_handle__stop__(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__stop__)

    call base_model__stop__(this%model)

    POP_SUB(base_handle__stop__)
  end subroutine base_handle__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_handle_stop)

    call base_handle_init(iter, this)
    do
      nullify(subs)
      call base_handle_next(iter, subs, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_stop(subs)
    end do
    call base_handle_end(iter)
    nullify(subs)
    call base_handle__stop__(this)

    POP_SUB(base_handle_stop)
  end subroutine base_handle_stop

  ! ---------------------------------------------------------
  subroutine base_handle__reset__(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__reset__)

    call base_model__reset__(this%model)

    POP_SUB(base_handle__reset__)
  end subroutine base_handle__reset__

  ! ---------------------------------------------------------
  subroutine base_handle__set__(this, that)
    type(base_handle_t),         intent(inout) :: this
    type(base_handle_t), target, intent(in)    :: that

    PUSH_SUB(base_handle__set__)

    this%prnt => that

    POP_SUB(base_handle__set__)
  end subroutine base_handle__set__

  ! ---------------------------------------------------------
  subroutine base_handle__sets__(this, that, config)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    PUSH_SUB(base_handle__sets__)

    call base_model_sets(this%model, that%model, config)

    POP_SUB(base_handle__sets__)
  end subroutine base_handle__sets__

  ! ---------------------------------------------------------
  subroutine base_handle_sets(this, that, config)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    character(len=BASE_HANDLE_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_handle_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_handle_hash_set(this%hash, config, that)
    call base_handle__sets__(this, that, config)

    POP_SUB(base_handle_sets)
  end subroutine base_handle_sets

  ! ---------------------------------------------------------
  subroutine base_handle_gets_config(this, config, that)
    type(base_handle_t),  intent(in) :: this
    type(json_object_t),  intent(in) :: config
    type(base_handle_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_handle_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_handle_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_HANDLE_OK) nullify(that)

    POP_SUB(base_handle_gets_config)
  end subroutine base_handle_gets_config

  ! ---------------------------------------------------------
  subroutine base_handle_gets_name(this, name, that)
    type(base_handle_t),  intent(in) :: this
    character(len=*),     intent(in) :: name
    type(base_handle_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_handle_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_handle_gets(this, config, that)

    POP_SUB(base_handle_gets_name)
  end subroutine base_handle_gets_name

  ! ---------------------------------------------------------
  subroutine base_handle_get_info(this, type)
    type(base_handle_t), intent(in)  :: this
    integer,  optional,  intent(out) :: type

    PUSH_SUB(base_handle_get_info)

    if(present(type)) type = this%type

    POP_SUB(base_handle_get_info)
  end subroutine base_handle_get_info

  ! ---------------------------------------------------------
  subroutine base_handle_get_config(this, that)
    type(base_handle_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_handle_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_handle_get_config)
  end subroutine base_handle_get_config

  ! ---------------------------------------------------------
  subroutine base_handle_get_simulation(this, that)
    type(base_handle_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_handle_get_simulation)

    that => this%sim

    POP_SUB(base_handle_get_simulation)
  end subroutine base_handle_get_simulation

  ! ---------------------------------------------------------
  subroutine base_handle_get_model(this, that)
    type(base_handle_t), target, intent(in) :: this
    type(base_model_t), pointer             :: that

    PUSH_SUB(base_handle_get_model)

    that => this%model

    POP_SUB(base_handle_get_model)
  end subroutine base_handle_get_model

  ! ---------------------------------------------------------
  subroutine base_handle__copy__begin(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(base_handle__copy__begin)

    call base_handle__end__(this)
    if(associated(that%config))then
      call base_handle__init__(this, that)
      if(simulation_assoc(that%sim)) call base_model__copy__(this%model, that%model)
    end if

    POP_SUB(base_handle__copy__begin)
  end subroutine base_handle__copy__begin

  ! ---------------------------------------------------------
  subroutine base_handle__copy__finish(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__copy__finish)

    call base_model__copy__(this%model)

    POP_SUB(base_handle__copy__finish)
  end subroutine base_handle__copy__finish

  ! ---------------------------------------------------------
  recursive subroutine base_handle_copy_type(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_handle_copy_type)

    nullify(cnfg, osub, isub)
    call base_handle_end(this)
    call base_handle__copy__(this, that)
    call base_handle_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_handle_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_new(this, osub)
      call base_handle_copy(osub, isub)
      call base_handle_sets(this, osub, cnfg)
    end do
    call base_handle_end(iter)
    call base_handle__copy__(this)
    nullify(cnfg, osub, isub)

    POP_SUB(base_handle_copy_type)
  end subroutine base_handle_copy_type

  ! ---------------------------------------------------------
  subroutine base_handle__end__(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(base_handle__end__)

    nullify(this%config, this%prnt)
    this%type = HNDL_TYPE_NONE
    call simulation_end(this%sim)
    call base_model__end__(this%model)
    call config_dict_end(this%dict)
    call base_handle_hash_end(this%hash)
    call base_handle_list_end(this%list)

    POP_SUB(base_handle__end__)
  end subroutine base_handle__end__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_end_type(this)
    type(base_handle_t), intent(inout) :: this

    type(base_handle_t), pointer :: subs

    PUSH_SUB(base_handle_end_type)

    do
      nullify(subs)
      call base_handle_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_handle_end(subs)
      call base_handle__del__(subs)
    end do
    nullify(subs)
    call base_handle__end__(this)

    POP_SUB(base_handle_end_type)
  end subroutine base_handle_end_type

#define TEMPLATE_PREFIX base_handle
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_handle_oct_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
