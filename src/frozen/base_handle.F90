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

  use json_m, only: JSON_OK, json_object_t, json_get

  use grid_m, only: &
    grid_t

  use base_model_m, only:              &
    model_t      => base_model_t,      &
    model_init   => base_model_init,   &
    model_start  => base_model_start,  &
    model_update => base_model_update, &
    model_stop   => base_model_stop,   &
    model_get    => base_model_get,    &
    model_copy   => base_model_copy,   &
    model_end    => base_model_end

  implicit none

  private

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
    type(model_t)                :: model
    type(base_handle_hash_t)     :: hash
  end type base_handle_t

  type, public :: base_handle_iterator_t
    private
    type(base_handle_t),      pointer :: self =>null()
    type(base_handle_hash_iterator_t) :: iter
  end type base_handle_iterator_t

  interface base_handle_init
    module procedure base_handle_init_begin
    module procedure base_handle_init_build
    module procedure base_handle_init_finish
    module procedure base_handle_iterator_init
  end interface base_handle_init

  interface base_handle_next
    module procedure base_handle_iterator_next_config_handle
    module procedure base_handle_iterator_next_config
    module procedure base_handle_iterator_next_handle
  end interface base_handle_next

  interface base_handle_get
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

  integer, public, parameter :: BASE_HANDLE_OK          = BASE_HANDLE_HASH_OK
  integer, public, parameter :: BASE_HANDLE_KEY_ERROR   = BASE_HANDLE_HASH_KEY_ERROR
  integer, public, parameter :: BASE_HANDLE_EMPTY_ERROR = BASE_HANDLE_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_handle_init_begin(this, config)
    type(base_handle_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_handle_init_begin)
    this%config=>config
    nullify(cnfg)
    call json_get(this%config, "model", cnfg, ierr)
    if(ierr==JSON_OK)call model_init(this%model, cnfg)
    nullify(cnfg)
    call base_handle_hash_init(this%hash)
    POP_SUB(base_handle_init_begin)
    return
  end subroutine base_handle_init_begin

  ! ---------------------------------------------------------
  subroutine base_handle_init_build(this, that, config)
    type(base_handle_t),         intent(inout) :: this
    type(base_handle_t),         intent(in)    :: that
    type(json_object_t), target, intent(in)    :: config
    !
    PUSH_SUB(base_handle_init_build)
    call base_handle_hash_set(this%hash, config, that)
    call model_init(this%model, that%model, config)
    POP_SUB(base_handle_init_build)
    return
  end subroutine base_handle_init_build

  ! ---------------------------------------------------------
  subroutine base_handle_init_finish(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle_init_finish)
    call model_init(this%model)
    POP_SUB(base_handle_init_finish)
    return
  end subroutine base_handle_init_finish

  ! ---------------------------------------------------------
  subroutine base_handle_start(this, grid)
    type(base_handle_t),    intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid
    !
    PUSH_SUB(base_handle_start)
    call model_start(this%model, grid)
    POP_SUB(base_handle_start)
    return
  end subroutine base_handle_start

  ! ---------------------------------------------------------
  subroutine base_handle_update(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle_update)
    call model_update(this%model)
    POP_SUB(base_handle_update)
    return
  end subroutine base_handle_update

  ! ---------------------------------------------------------
  subroutine base_handle_stop(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle_stop)
    call model_stop(this%model)
    POP_SUB(base_handle_stop)
    return
  end subroutine base_handle_stop

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
    type(model_t),      pointer             :: that
    !
    PUSH_SUB(base_handle_get_model)
    that=>this%model
    POP_SUB(base_handle_get_model)
    return
  end subroutine base_handle_get_model

  ! ---------------------------------------------------------
  subroutine base_handle_copy_handle(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that
    !
    PUSH_SUB(base_handle_copy_handle)
    this%config=>that%config
    call base_handle_hash_copy(this%hash, that%hash)
    call model_copy(this%model, that%model)
    POP_SUB(base_handle_copy_handle)
    return
  end subroutine base_handle_copy_handle

  ! ---------------------------------------------------------
  subroutine base_handle_end_handle(this)
    type(base_handle_t), intent(inout) :: this
    !
    PUSH_SUB(base_handle_end_handle)
    call base_handle_hash_end(this%hash)
    call model_end(this%model)
    nullify(this%config)
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
