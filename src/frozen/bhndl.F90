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

#define HASH_TEMPLATE_NAME bhndl
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bhndl

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module bhndl_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m, only: JSON_OK, json_object_t, json_array_t, json_get
  use json_m, only: json_array_iterator_t, json_init, json_next, json_end

  use grid_m, only: &
    grid_t

  use bmodl_m, only:              &
    model_t      => bmodl_t,      &
    model_init   => bmodl_init,   &
    model_start  => bmodl_start,  &
    model_update => bmodl_update, &
    model_get    => bmodl_get,    &
    model_copy   => bmodl_copy,   &
    model_end    => bmodl_end

  implicit none

  private

  public ::       &
    bhndl_init,   &
    bhndl_start,  &
    bhndl_update, &
    bhndl_next,   &
    bhndl_get,    &
    bhndl_copy,   &
    bhndl_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: bhndl_t
    private
    type(json_object_t), pointer :: config =>null()
    type(model_t)                :: model
    type(bhndl_hash_t)           :: hash
  end type bhndl_t

  type, public :: bhndl_iterator_t
    private
    type(bhndl_t),      pointer :: self =>null()
    type(bhndl_hash_iterator_t) :: iter
  end type bhndl_iterator_t

  interface bhndl_init
    module procedure bhndl_init_bhndl
    module procedure bhndl_iterator_init
  end interface bhndl_init

  interface bhndl_next
    module procedure bhndl_iterator_next_config_bhndl
    module procedure bhndl_iterator_next_config
    module procedure bhndl_iterator_next_bhndl
  end interface bhndl_next

  interface bhndl_get
    module procedure bhndl_get_config
    module procedure bhndl_get_model
  end interface bhndl_get

  interface bhndl_copy
    module procedure bhndl_copy_bhndl
    module procedure bhndl_iterator_copy
  end interface bhndl_copy

  interface bhndl_end
    module procedure bhndl_end_bhndl
    module procedure bhndl_iterator_end
  end interface bhndl_end

  integer, public, parameter :: BHNDL_OK          = BHNDL_HASH_OK
  integer, public, parameter :: BHNDL_KEY_ERROR   = BHNDL_HASH_KEY_ERROR
  integer, public, parameter :: BHNDL_EMPTY_ERROR = BHNDL_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  recursive subroutine bhndl_init_bhndl(this, config)
    type(bhndl_t),               intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    integer                      :: ierr
    !
    PUSH_SUB(bhndl_init_bhndl)
    this%config=>config
    nullify(cnfg, list)
    call json_get(this%config, "model", cnfg, ierr)
    if(ierr==JSON_OK)call model_init(this%model, cnfg)
    nullify(cnfg)
    call bhndl_hash_init(this%hash)
    call json_get(this%config, "systems", list, ierr)
    if(ierr==JSON_OK)call bhndl_init_hash(this%hash, this%model, list)
    nullify(list)
    call model_init(this%model)
    POP_SUB(bhndl_init_bhndl)
    return
  end subroutine bhndl_init_bhndl

  ! ---------------------------------------------------------
  recursive subroutine bhndl_init_hash(this, model, list)
    use json_m
    type(bhndl_hash_t), intent(inout) :: this
    type(model_t),      intent(inout) :: model
    type(json_array_t), intent(in)    :: list
    !
    type(bhndl_t),       pointer :: hndl
    type(json_object_t), pointer :: cnfg
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr
    !
    PUSH_SUB(bhndl_init_hash)
    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      SAFE_ALLOCATE(hndl)
      call bhndl_init(hndl, cnfg)
      call model_init(model, hndl%model, cnfg)
      call bhndl_hash_set(this, cnfg, hndl)
    end do
    call json_end(iter)
    nullify(cnfg)
    POP_SUB(bhndl_init_hash)
    return
  end subroutine bhndl_init_hash

  ! ---------------------------------------------------------
  recursive subroutine bhndl_start(this, grid)
    type(bhndl_t),          intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid
    !
    type(bhndl_iterator_t) :: iter
    type(bhndl_t), pointer :: hndl
    integer                :: ierr
    !
    PUSH_SUB(bhndl_start)
    call bhndl_init(iter, this)
    do
      nullify(hndl)
      call bhndl_next(iter, hndl, ierr)
      if(ierr/=BHNDL_OK)exit
      call bhndl_start(hndl, grid)
    end do
    call bhndl_end(iter)
    nullify(hndl)
    call model_start(this%model, grid)
    POP_SUB(bhndl_start)
    return
  end subroutine bhndl_start

  ! ---------------------------------------------------------
  recursive subroutine bhndl_update(this)
    type(bhndl_t), intent(inout) :: this
    !
    type(bhndl_iterator_t) :: iter
    type(bhndl_t), pointer :: hndl
    integer                :: ierr
    !
    PUSH_SUB(bhndl_update)
    call bhndl_init(iter, this)
    do
      nullify(hndl)
      call bhndl_next(iter, hndl, ierr)
      if(ierr/=BHNDL_HASH_OK)exit
      call bhndl_update(hndl)
    end do
    call bhndl_end(iter)
    nullify(hndl)
    call model_update(this%model)
    POP_SUB(bhndl_update_hash)
    return
  end subroutine bhndl_update

  ! ---------------------------------------------------------
  subroutine bhndl_get_config(this, that)
    type(bhndl_t),        target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(bhndl_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(bhndl_get_config)
    return
  end subroutine bhndl_get_config

  ! ---------------------------------------------------------
  subroutine bhndl_get_model(this, that)
    type(bhndl_t), target, intent(in) :: this
    type(model_t), pointer             :: that
    !
    PUSH_SUB(bhndl_get_model)
    that=>this%model
    POP_SUB(bhndl_get_model)
    return
  end subroutine bhndl_get_model

  ! ---------------------------------------------------------
  recursive subroutine bhndl_copy_bhndl(this, that)
    type(bhndl_t), intent(inout) :: this
    type(bhndl_t), intent(in)    :: that
    !
    PUSH_SUB(bhndl_copy_bhndl)
    call bhndl_end(this)
    this%config=>that%config
    call model_copy(this%model, that%model)
    call bhndl_copy_hash(this, that)
    POP_SUB(bhndl_copy_bhndl)
    return
  end subroutine bhndl_copy_bhndl

  ! ---------------------------------------------------------
  recursive subroutine bhndl_copy_hash(this, that)
    type(bhndl_t), intent(inout) :: this
    type(bhndl_t), intent(in)    :: that
    !
    type(bhndl_iterator_t)       :: iter
    type(json_object_t), pointer :: cnfg
    type(bhndl_t),       pointer :: ihnd, ohnd
    integer                      :: ierr
    !
    PUSH_SUB(bhndl_copy_hash)
    call bhndl_init(iter, this)
    do
      nullify(cnfg, ihnd, ohnd)
      call bhndl_next(iter, cnfg, ihnd, ierr)
      if(ierr/=BHNDL_OK)exit
      SAFE_ALLOCATE(ohnd)
      call bhndl_copy(ohnd, ihnd)
      call bhndl_hash_set(this%hash, cnfg, ohnd)
    end do
    call bhndl_end(iter)
    nullify(cnfg, ihnd, ohnd)
    POP_SUB(bhndl_copy_hash)
    return
  end subroutine bhndl_copy_hash

  ! ---------------------------------------------------------
  recursive subroutine bhndl_end_bhndl(this)
    type(bhndl_t), intent(inout) :: this
    !
    PUSH_SUB(bhndl_end_bhndl)
    call bhndl_end_hash(this%hash)
    call model_end(this%model)
    nullify(this%config)
    POP_SUB(bhndl_end_bhndl)
    return
  end subroutine bhndl_end_bhndl

  ! ---------------------------------------------------------
  recursive subroutine bhndl_end_hash(this)
    type(bhndl_hash_t), intent(inout) :: this
    !
    type(bhndl_t), pointer :: hndl
    !
    PUSH_SUB(bhndl_end_hash)
    do
      nullify(hndl)
      call bhndl_hash_pop(this, hndl)
      if(.not.associated(hndl))exit
      call bhndl_end(hndl)
      SAFE_DEALLOCATE_P(hndl)
    end do
    call bhndl_hash_end(this)
    POP_SUB(bhndl_end_hash)
    return
  end subroutine bhndl_end_hash

  ! ---------------------------------------------------------
  subroutine bhndl_iterator_init(this, that)
    type(bhndl_iterator_t), intent(out) :: this
    type(bhndl_t),  target, intent(in)  :: that
    !
    PUSH_SUB(bhndl_iterator_init)
    this%self=>that
    call bhndl_hash_init(this%iter, that%hash)
    POP_SUB(bhndl_iterator_init)
    return
  end subroutine bhndl_iterator_init

  ! ---------------------------------------------------------
  subroutine bhndl_iterator_next_config_bhndl(this, config, handle, ierr)
    type(bhndl_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: config
    type(bhndl_t),         pointer        :: handle
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bhndl_iterator_next_config_bhndl)
    call bhndl_hash_next(this%iter, config, handle, ierr)
    POP_SUB(bhndl_iterator_next_config_bhndl)
    return
  end subroutine bhndl_iterator_next_config_bhndl

  ! ---------------------------------------------------------
  subroutine bhndl_iterator_next_config(this, that, ierr)
    type(bhndl_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bhndl_iterator_next_config)
    call bhndl_hash_next(this%iter, that, ierr)
    POP_SUB(bhndl_iterator_next_config)
    return
  end subroutine bhndl_iterator_next_config

  ! ---------------------------------------------------------
  subroutine bhndl_iterator_next_bhndl(this, that, ierr)
    type(bhndl_iterator_t), intent(inout) :: this
    type(bhndl_t),         pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bhndl_iterator_next_bhndl)
    call bhndl_hash_next(this%iter, that, ierr)
    POP_SUB(bhndl_iterator_next_bhndl)
    return
  end subroutine bhndl_iterator_next_bhndl

  ! ---------------------------------------------------------
  subroutine bhndl_iterator_copy(this, that)
    type(bhndl_iterator_t), intent(inout) :: this
    type(bhndl_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(bhndl_iterator_copy)
    this%self=>that%self
    call bhndl_hash_copy(this%iter, that%iter)
    POP_SUB(bhndl_iterator_copy)
    return
  end subroutine bhndl_iterator_copy

  ! ---------------------------------------------------------
  subroutine bhndl_iterator_end(this)
    type(bhndl_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(bhndl_iterator_end)
    nullify(this%self)
    call bhndl_hash_end(this%iter)
    POP_SUB(bhndl_iterator_end)
    return
  end subroutine bhndl_iterator_end

end module bhndl_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
