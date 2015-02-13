#include "global.h"

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

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_handle_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m, only: JSON_OK, json_object_t, json_init, json_get, json_end
  use json_m, only: json_array_t, json_array_iterator_t, json_next

  use config_dict_m, only: &
    CONFIG_DICT_OK,        &
    CONFIG_DICT_NAME_LEN

  use config_dict_m, only: &
    config_dict_t,         &
    config_dict_len,       &
    config_dict_init,      &
    config_dict_set,       &
    config_dict_get,       &
    config_dict_end

  use igrid_m, only: &
    grid_t

  use base_model_m, only: &
    base_model__init__,   &
    base_model__start__,  &
    base_model__update__, &
    base_model__stop__,   &
    base_model__add__

  use base_model_m, only: &
    base_model_t,         &
    base_model_get,       &
    base_model_copy,      &
    base_model_end

  implicit none

  private
  public ::                &
    base_handle__init__,   &
    base_handle__start__,  &
    base_handle__update__, &
    base_handle__stop__,   &
    base_handle__add__,    &
    base_handle__get__,    &
    base_handle__copy__,   &
    base_handle__end__

  public ::             &
    base_handle_init,   &
    base_handle_start,  &
    base_handle_update, &
    base_handle_stop,   &
    base_handle_next,   &
    base_handle_get,    &
    base_handle_copy,   &
    base_handle_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_handle_t
    private
    type(json_object_t), pointer :: config =>null()
    integer                      :: type
    type(base_model_t)           :: model
    type(config_dict_t)          :: dict
    type(base_handle_hash_t)     :: hash
  end type base_handle_t

  type, public :: base_handle_iterator_t
    private
    type(base_handle_t),      pointer :: self =>null()
    type(base_handle_hash_iterator_t) :: iter
  end type base_handle_iterator_t

  interface base_handle__init__
    module procedure base_handle__init__begin
    module procedure base_handle__init__finish
  end interface base_handle__init__

  interface base_handle_init
    module procedure base_handle_init_handle
    module procedure base_handle_init_pass
    module procedure base_handle_iterator_init
  end interface base_handle_init

  interface base_handle_next
    module procedure base_handle_iterator_next_config_handle
    module procedure base_handle_iterator_next_config
    module procedure base_handle_iterator_next_handle
  end interface base_handle_next

  interface base_handle_get
    module procedure base_handle_get_type
    module procedure base_handle_get_config
    module procedure base_handle_get_model
  end interface base_handle_get

  interface base_handle_copy
    module procedure base_handle_copy_handle
    module procedure base_handle_iterator_copy
  end interface base_handle_copy

  interface base_handle_end
    module procedure base_handle_end_handle
    module procedure base_handle_iterator_end
  end interface base_handle_end

  integer, public, parameter :: HNDL_TYPE_NONE = 0

  integer, public, parameter :: HNDL_NAME_LEN = CONFIG_DICT_NAME_LEN

  integer, public, parameter :: BASE_HANDLE_OK          = BASE_HANDLE_HASH_OK
  integer, public, parameter :: BASE_HANDLE_KEY_ERROR   = BASE_HANDLE_HASH_KEY_ERROR
  integer, public, parameter :: BASE_HANDLE_EMPTY_ERROR = BASE_HANDLE_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_handle__init__begin(this, config)
    type(base_handle_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_handle__init__begin)
    this%config=>config
    nullify(cnfg)
    call json_get(this%config, "type", this%type, ierr)
    if(ierr/=JSON_OK)this%type=HNDL_TYPE_NONE
    call json_get(this%config, "model", cnfg, ierr)
    if(ierr==JSON_OK)call base_model__init__(this%model, cnfg)
    nullify(cnfg)
    call config_dict_init(this%dict)
    call base_handle_hash_init(this%hash)
    POP_SUB(base_handle__init__begin)
    return
  end subroutine base_handle__init__begin

  ! ---------------------------------------------------------
  subroutine base_handle__init__finish(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__init__finish)
    call base_model__init__(this%model)
    POP_SUB(base_handle__init__finish)
    return
  end subroutine base_handle__init__finish

  ! ---------------------------------------------------------
  recursive subroutine base_handle_init_handle(this, config)
    type(base_handle_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_array_iterator_t)   :: iter
    type(json_object_t),  pointer :: cnfg
    type(json_array_t),   pointer :: list
    type(base_handle_t),  pointer :: hndl
    integer                       :: ierr
    !
    PUSH_SUB(base_handle_init_handle)
    nullify(cnfg, list, hndl)
    call base_handle__init__(this, config)
    call json_get(config, "systems", list, ierr)
    if(ierr==JSON_OK)then
      call json_init(iter, list)
      do
        nullify(cnfg, hndl)
        call json_next(iter, cnfg, ierr)
        if(ierr/=JSON_OK)exit
        SAFE_ALLOCATE(hndl)
        call base_handle_init(hndl, cnfg)
        call base_handle__add__(this, hndl, cnfg)
      end do
      call json_end(iter)
    end if
    call base_handle__init__(this)
    nullify(cnfg, list, hndl)
    POP_SUB(base_handle_init_handle)
    return
  end subroutine base_handle_init_handle

  ! ---------------------------------------------------------
  subroutine base_handle_init_pass(this, handle_init, config)
    type(base_handle_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_array_iterator_t)   :: iter
    type(json_object_t),  pointer :: cnfg
    type(json_array_t),   pointer :: list
    type(base_handle_t),  pointer :: hndl
    integer                       :: ierr
    !
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
        SAFE_ALLOCATE(hndl)
        call handle_init(hndl, cnfg)
        call base_handle__add__(this, hndl, cnfg)
      end do
      call json_end(iter)
    end if
    call base_handle__init__(this)
    nullify(cnfg, list, hndl)
    POP_SUB(base_handle_init_pass)
    return
  end subroutine base_handle_init_pass

  ! ---------------------------------------------------------
  subroutine base_handle__start__(this, grid)
    type(base_handle_t),    intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid
    !
    PUSH_SUB(base_handle_start)
    call base_model__start__(this%model, grid)
    POP_SUB(base_handle_start)
    return
  end subroutine base_handle__start__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_start(this, grid)
    type(base_handle_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid
    !
    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: subs
    integer                      :: ierr
    !
    PUSH_SUB(base_handle_start)
    nullify(subs)
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
    return
  end subroutine base_handle_start

  ! ---------------------------------------------------------
  subroutine base_handle__update__(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__update__)
    ASSERT(associated(this%config))
    call base_model__update__(this%model)
    POP_SUB(base_handle__update__)
    return
  end subroutine base_handle__update__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_update(this)
    type(base_handle_t), intent(inout) :: this
    !
    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: subs
    integer                      :: ierr
    !
    PUSH_SUB(base_handle_update)
    nullify(subs)
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
    return
  end subroutine base_handle_update

  ! ---------------------------------------------------------
  subroutine base_handle__stop__(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__stop__)
    call base_model__stop__(this%model)
    POP_SUB(base_handle__stop__)
    return
  end subroutine base_handle__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_stop(this)
    type(base_handle_t), intent(inout) :: this
    !
    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: subs
    integer                      :: ierr
    !
    PUSH_SUB(base_handle_stop)
    nullify(subs)
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
    return
  end subroutine base_handle_stop

  ! ---------------------------------------------------------
  subroutine base_handle__add__(this, that, config)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    character(len=HNDL_NAME_LEN) :: name
    integer                      :: ierr
    !
    PUSH_SUB(base_handle__add__)
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_handle_hash_set(this%hash, config, that)
    call base_model__add__(this%model, that%model, config)
    POP_SUB(base_handle__add__)
    return
  end subroutine base_handle__add__

  ! ---------------------------------------------------------
  subroutine base_handle__get__(this, name, that)
    type(base_handle_t),  intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(base_handle_t), pointer        :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_handle__get__)
    nullify(that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_handle_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_HANDLE_OK)nullify(that)
    end if
    POP_SUB(base_handle__get__)
    return
  end subroutine base_handle__get__

  ! ---------------------------------------------------------
  subroutine base_handle_get_type(this, that)
    type(base_handle_t), intent(in)  :: this
    integer,             intent(out) :: that
    !
    PUSH_SUB(base_handle_get_type)
    that=this%type
    POP_SUB(base_handle_get_type)
    return
  end subroutine base_handle_get_type

  ! ---------------------------------------------------------
  subroutine base_handle_get_config(this, that)
    type(base_handle_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(base_handle_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_handle_get_config)
    return
  end subroutine base_handle_get_config

  ! ---------------------------------------------------------
  subroutine base_handle_get_model(this, that)
    type(base_handle_t), target, intent(in) :: this
    type(base_model_t), pointer             :: that
    !
    PUSH_SUB(base_handle_get_model)
    that=>this%model
    POP_SUB(base_handle_get_model)
    return
  end subroutine base_handle_get_model

  ! ---------------------------------------------------------
  subroutine base_handle__copy__(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that
    !
    PUSH_SUB(base_handle__copy__)
    this%config=>that%config
    this%type=that%type
    call base_model_copy(this%model, that%model)
    call config_dict_init(this%dict, config_dict_len(that%dict))
    call base_handle_hash_init(this%hash, base_handle_hash_len(that%hash))
    POP_SUB(base_handle__copy__)
    return
  end subroutine base_handle__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_copy_handle(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that
    !
    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_handle_copy_handle)
    nullify(cnfg, osub, isub)
    call base_handle_end(this)
    call base_handle__copy__(this, that)
    call base_handle_init(iter, this)
    do
      nullify(cnfg, osub, isub)
      call base_handle_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      SAFE_ALLOCATE(osub)
      call base_handle_copy(osub, isub)
      call base_handle__add__(this, osub, cnfg)
    end do
    call base_handle_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(base_handle_copy_handle)
    return
  end subroutine base_handle_copy_handle

  ! ---------------------------------------------------------
  subroutine base_handle__end__(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__end__)
    nullify(this%config)
    this%type=HNDL_TYPE_NONE
    call base_model_end(this%model)
    call config_dict_end(this%dict)
    call base_handle_hash_end(this%hash)
    POP_SUB(base_handle__end__)
    return
  end subroutine base_handle__end__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_end_handle(this)
    type(base_handle_t), intent(inout) :: this
    !
    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: subs
    integer                      :: ierr
    !
    PUSH_SUB(base_handle_end_handle)
    nullify(subs)
    call base_handle_init(iter, this)
    do
      nullify(subs)
      call base_handle_next(iter, subs, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_end(subs)
      SAFE_DEALLOCATE_P(subs)
    end do
    call base_handle_end(iter)
    nullify(subs)
    call base_handle__end__(this)
    POP_SUB(base_handle_end_handle)
    return
  end subroutine base_handle_end_handle

  ! ---------------------------------------------------------
  subroutine base_handle_iterator_init(this, that)
    type(base_handle_iterator_t), intent(out) :: this
    type(base_handle_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_handle_iterator_init)
    this%self=>that
    call base_handle_hash_init(this%iter, that%hash)
    POP_SUB(base_handle_iterator_init)
    return
  end subroutine base_handle_iterator_init

  ! ---------------------------------------------------------
  subroutine base_handle_iterator_next_config_handle(this, config, that, ierr)
    type(base_handle_iterator_t), intent(inout) :: this
    type(json_object_t),         pointer        :: config
    type(base_handle_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_handle_iterator_next_config_handle)
    call base_handle_hash_next(this%iter, config, that, ierr)
    POP_SUB(base_handle_iterator_next_config_handle)
    return
  end subroutine base_handle_iterator_next_config_handle

  ! ---------------------------------------------------------
  subroutine base_handle_iterator_next_config(this, that, ierr)
    type(base_handle_iterator_t), intent(inout) :: this
    type(json_object_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_handle_iterator_next_config)
    call base_handle_hash_next(this%iter, that, ierr)
    POP_SUB(base_handle_iterator_next_config)
    return
  end subroutine base_handle_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_handle_iterator_next_handle(this, that, ierr)
    type(base_handle_iterator_t), intent(inout) :: this
    type(base_handle_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_handle_iterator_next_handle)
    call base_handle_hash_next(this%iter, that, ierr)
    POP_SUB(base_handle_iterator_next_handle)
    return
  end subroutine base_handle_iterator_next_handle

  ! ---------------------------------------------------------
  subroutine base_handle_iterator_copy(this, that)
    type(base_handle_iterator_t), intent(inout) :: this
    type(base_handle_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_handle_iterator_copy)
    this%self=>that%self
    call base_handle_hash_copy(this%iter, that%iter)
    POP_SUB(base_handle_iterator_copy)
    return
  end subroutine base_handle_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_handle_iterator_end(this)
    type(base_handle_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle_iterator_end)
    nullify(this%self)
    call base_handle_hash_end(this%iter)
    POP_SUB(base_handle_iterator_end)
    return
  end subroutine base_handle_iterator_end

end module base_handle_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
