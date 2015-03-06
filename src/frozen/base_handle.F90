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

#define LIST_TEMPLATE_NAME base_handle
#define LIST_INCLUDE_PREFIX
#include "tlist.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

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
    base_model__add__,    &
    base_model__copy__,   &
    base_model__end__

  use base_model_m, only: &
    base_model_t,         &
    base_model_get

  implicit none

  private
  public ::                &
    base_handle__init__,   &
    base_handle__start__,  &
    base_handle__update__, &
    base_handle__stop__,   &
    base_handle__add__,    &
    base_handle__copy__,   &
    base_handle__end__

  public ::             &
    base_handle_new,    &
    base_handle_del,    &
    base_handle_init,   &
    base_handle_start,  &
    base_handle_update, &
    base_handle_stop,   &
    base_handle_next,   &
    base_handle_get,    &
    base_handle_copy,   &
    base_handle_end

#define LIST_TEMPLATE_NAME base_handle
#define LIST_INCLUDE_HEADER
#include "tlist.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  integer, public, parameter :: HNDL_TYPE_NONE = 0

  type, public :: base_handle_t
    private
    type(json_object_t), pointer :: config =>null()
    type(base_handle_t), pointer :: prnt   =>null()
    integer                      :: type   = HNDL_TYPE_NONE
    type(base_model_t)           :: model
    type(config_dict_t)          :: dict
    type(base_handle_hash_t)     :: hash
    type(base_handle_list_t)     :: list
  end type base_handle_t

  type, public :: base_handle_iterator_t
    private
    type(base_handle_t),      pointer :: self =>null()
    type(base_handle_hash_iterator_t) :: iter
  end type base_handle_iterator_t

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
    module procedure base_handle_init_handle
    module procedure base_handle_init_pass
    module procedure base_handle_init_copy
    module procedure base_handle_iterator_init
  end interface base_handle_init

  interface base_handle_next
    module procedure base_handle_iterator_next_config_handle
    module procedure base_handle_iterator_next_config
    module procedure base_handle_iterator_next_handle
  end interface base_handle_next

  interface base_handle_get
    module procedure base_handle_get_handle
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

  integer, public, parameter :: HNDL_NAME_LEN = CONFIG_DICT_NAME_LEN

  integer, public, parameter :: BASE_HANDLE_OK          = BASE_HANDLE_HASH_OK
  integer, public, parameter :: BASE_HANDLE_KEY_ERROR   = BASE_HANDLE_HASH_KEY_ERROR
  integer, public, parameter :: BASE_HANDLE_EMPTY_ERROR = BASE_HANDLE_HASH_EMPTY_ERROR

contains

#define LIST_TEMPLATE_NAME base_handle
#define LIST_INCLUDE_BODY
#include "tlist.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  subroutine base_handle_new(this, that)
    type(base_handle_t),  intent(inout) :: this
    type(base_handle_t), pointer        :: that
    !
    PUSH_SUB(base_handle_new)
    nullify(that)
    SAFE_ALLOCATE(that)
    call base_handle__set__(that, this)
    call base_handle_list_push(this%list, that)
    POP_SUB(base_handle_new)
    return
  end subroutine base_handle_new

  ! ---------------------------------------------------------
  subroutine base_handle__idel__(this)
    type(base_handle_t), pointer :: this
    !
    PUSH_SUB(base_handle__idel__)
    SAFE_DEALLOCATE_P(this)
    nullify(this)
    POP_SUB(base_handle__idel__)
    return
  end subroutine base_handle__idel__

  ! ---------------------------------------------------------
  subroutine base_handle_del(this)
    type(base_handle_t), pointer :: this
    !
    PUSH_SUB(base_handle_del)
    if(associated(this))then
      if(associated(this%prnt))then
        call base_handle_list_del(this%prnt%list, this)
        call base_handle_end(this)
        call base_handle__idel__(this)
      end if
    end if
    POP_SUB(base_handle_del)
    return
  end subroutine base_handle_del

  ! ---------------------------------------------------------
  subroutine base_handle__inull__(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__inull__)
    nullify(this%config, this%prnt)
    this%type=HNDL_TYPE_NONE
    POP_SUB(base_handle__inull__)
    return
  end subroutine base_handle__inull__

  ! ---------------------------------------------------------
  subroutine base_handle__iinit__(this, config)
    type(base_handle_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: ierr
    !
    PUSH_SUB(base_handle__iinit__)
    call base_handle__inull__(this)
    this%config=>config
    call json_get(this%config, "type", this%type, ierr)
    if(ierr/=JSON_OK)this%type=HNDL_TYPE_NONE
    call config_dict_init(this%dict)
    call base_handle_hash_init(this%hash)
    call base_handle_list_init(this%list)
    POP_SUB(base_handle__iinit__)
    return
  end subroutine base_handle__iinit__

  ! ---------------------------------------------------------
  subroutine base_handle__init__begin(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_handle__init__begin)
    nullify(cnfg)
    call base_handle__iinit__(this, config)
    call json_get(this%config, "model", cnfg, ierr)
    if(ierr==JSON_OK)call base_model__init__(this%model, cnfg)
    nullify(cnfg)
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
  subroutine base_handle__init__copy(this, that)
    type(base_handle_t), intent(out) :: this
    type(base_handle_t), intent(in)  :: that
    !
    PUSH_SUB(base_handle__init__copy)
    if(associated(that%config))then
      call base_handle__iinit__(this, that%config)
      call base_model__init__(this%model, that%model)
    end if
    POP_SUB(base_handle__init__copy)
    return
  end subroutine base_handle__init__copy

  ! ---------------------------------------------------------
  recursive subroutine base_handle_init_handle(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config
    !
    PUSH_SUB(base_handle_init_handle)
    call base_handle_init_pass(this, config, base_handle_init_handle)
    POP_SUB(base_handle_init_handle)
    return
  end subroutine base_handle_init_handle

  ! ---------------------------------------------------------
  subroutine base_handle_init_pass(this, config, handle_init)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config
    interface
      subroutine handle_init(this, config)
        use json_m, only: json_object_t
        import :: base_handle_t
        type(base_handle_t), intent(out) :: this
        type(json_object_t), intent(in)  :: config
      end subroutine handle_init
    end interface
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
        call base_handle_new(this, hndl)
        call handle_init(hndl, cnfg)
        call base_handle__set__(hndl, this)
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
  recursive subroutine base_handle_init_copy(this, that)
    type(base_handle_t), intent(out) :: this
    type(base_handle_t), intent(in)  :: that
    !
    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_handle_init_copy)
    nullify(cnfg, osub, isub)
    call base_handle__init__(this, that)
    call base_handle_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_handle_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_new(this, osub)
      call base_handle_init(osub, isub)
      call base_handle__add__(this, osub, cnfg)
    end do
    call base_handle_end(iter)
    call base_handle__init__(this)
    nullify(cnfg, osub, isub)
    POP_SUB(base_handle_init_copy)
    return
  end subroutine base_handle_init_copy

  ! ---------------------------------------------------------
  subroutine base_handle__istart__(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__istart__)
    ASSERT(associated(this%config))
    POP_SUB(base_handle__istart__)
    return
  end subroutine base_handle__istart__

  ! ---------------------------------------------------------
  subroutine base_handle__start__(this, grid)
    type(base_handle_t),    intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid
    !
    PUSH_SUB(base_handle__start__)
    call base_handle__istart__(this)
    call base_model__start__(this%model, grid)
    POP_SUB(base_handle__start__)
    return
  end subroutine base_handle__start__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_start(this, grid)
    type(base_handle_t),    intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid
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
    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_handle_hash_set(this%hash, config, that)
    call base_model__add__(this%model, that%model, config)
    POP_SUB(base_handle__add__)
    return
  end subroutine base_handle__add__

  ! ---------------------------------------------------------
  subroutine base_handle__set__(this, that)
    type(base_handle_t),         intent(inout) :: this
    type(base_handle_t), target, intent(in)    :: that
    !
    PUSH_SUB(base_handle__set__)
    this%prnt=>that
    POP_SUB(base_handle__set__)
    return
  end subroutine base_handle__set__

  ! ---------------------------------------------------------
  subroutine base_handle_get_handle(this, name, that)
    type(base_handle_t),  intent(in) :: this
    character(len=*),     intent(in) :: name
    type(base_handle_t), pointer     :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_handle_get_handle)
    nullify(that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_handle_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_HANDLE_OK)nullify(that)
    end if
    POP_SUB(base_handle_get_handle)
    return
  end subroutine base_handle_get_handle

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
  subroutine base_handle__icopy__(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that
    !
    PUSH_SUB(base_handle__icopy__)
    call base_handle__iend__(this)
    if(associated(that%config))then
      call base_handle__iinit__(this, that%config)
      call base_handle__istart__(this)
    end if
    POP_SUB(base_handle__icopy__)
    return
  end subroutine base_handle__icopy__

  ! ---------------------------------------------------------
  subroutine base_handle__copy__begin(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that
    !
    PUSH_SUB(base_handle__copy__begin)
    call base_handle__icopy__(this, that)
    call base_model__copy__(this%model, that%model)
    POP_SUB(base_handle__copy__begin)
    return
  end subroutine base_handle__copy__begin

  ! ---------------------------------------------------------
  subroutine base_handle__copy__finish(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__copy__finish)
    call base_model__copy__(this%model)
    POP_SUB(base_handle__copy__finish)
    return
  end subroutine base_handle__copy__finish

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
    call base_handle_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_handle_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_new(this, osub)
      call base_handle_copy(osub, isub)
      call base_handle__add__(this, osub, cnfg)
    end do
    call base_handle_end(iter)
    call base_handle__copy__(this)
    nullify(cnfg, osub, isub)
    POP_SUB(base_handle_copy_handle)
    return
  end subroutine base_handle_copy_handle

  ! ---------------------------------------------------------
  subroutine base_handle__iend__(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__iend__)
    call base_handle__inull__(this)
    call config_dict_end(this%dict)
    call base_handle_hash_end(this%hash)
    call base_handle_list_end(this%list)
    POP_SUB(base_handle__iend__)
    return
  end subroutine base_handle__iend__

  ! ---------------------------------------------------------
  subroutine base_handle__end__(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle__end__)
    call base_handle__iend__(this)
    call base_model__end__(this%model)
    POP_SUB(base_handle__end__)
    return
  end subroutine base_handle__end__

  ! ---------------------------------------------------------
  recursive subroutine base_handle_end_handle(this)
    type(base_handle_t), intent(inout) :: this
    !
    type(base_handle_t), pointer :: subs
    !
    PUSH_SUB(base_handle_end_handle)
    do
      nullify(subs)
      call base_handle_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_handle_end(subs)
      call base_handle__idel__(subs)
    end do
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
